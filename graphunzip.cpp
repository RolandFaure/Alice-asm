#include "graphunzip.h"
#include "basic_graph_manipulation.h"

#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::mutex;
using std::pair;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::stringstream;
using robin_hood::unordered_map;
using std::set;

/*
 * Read paths (GAF) storage and per-segment-end consensus computation.
 *
 * Semantics reproduced exactly from the previous implementation (Segment::add_neighbor / compute_consensuses /
 * get_strong_neighbors_*), for one end of a segment:
 *  - N = the list of "neighbor paths" of that end, in GAF order: for every read sub-path (a read path is cut wherever two
 *    consecutive segments are not linked) and every occurrence of the segment in it, the rest of the sub-path on that side
 *    (oriented away from the segment). Empty neighbor paths are ignored.
 *  - N is std::sort-ed by decreasing length (unstable: equal-length paths are ordered by whatever permutation std::sort
 *    produces; sorting references with the same comparator reproduces that permutation since it only depends on the
 *    comparison outcomes).
 *  - maximal paths = the distinct paths of N that are not a prefix of another path, ordered by their first occurrence in
 *    the sorted N.
 *  - for a maximal path P and a position c: first = #paths of N that have P[0..c] as prefix,
 *    second = #paths of N that have P[0..c-1] as prefix, are longer than c and differ from P at c.
 *  - consensus: the first maximal path (in order) for which every position satisfies 0.2*(first-1) > second-1 && second < 5;
 *    the consensus is its prefix made of the positions with first > 1. The end is "haploid" if such a path exists or if N is
 *    empty (then the consensus is empty).
 *  - strong neighbors: for each maximal path in order, its longest prefix whose positions all satisfy
 *    first >= min_coverage && (second == 0 || first > 4); kept if non-empty.
 *
 * Instead of materialising every neighbor path (O(L^2) memory per read of L steps) and comparing all of them pairwise, the
 * GAF is parsed once into a binary file D of int32 steps (segment ID << 1 | orientation), which is memory-mapped. D holds
 * the forward sub-paths X followed by the reverse-complement R of X, so that every neighbor path is a contiguous range of
 * D (PathRef). For each segment end, a compressed trie of its neighbor paths is built from PathRefs (edges are ranges of D,
 * edge comparisons are accelerated with prefix hashes stored in the same file), giving first/second for all positions of
 * all maximal paths in time and memory linear in the number of neighbor paths.
 */
namespace {

struct PathRef{
    uint64_t start = 0; //index in D of the first step
    uint32_t len = 0;   //number of steps
};

const uint64_t HASH_MOD = (1ULL << 61) - 1;
const uint64_t HASH_BASE = 0x1f3a5c7e9b2d4f61ULL % HASH_MOD;

inline uint64_t mulmod(uint64_t a, uint64_t b){
    __uint128_t p = (__uint128_t) a * b;
    uint64_t r = (uint64_t) (p & HASH_MOD) + (uint64_t) (p >> 61);
    r = (r & HASH_MOD) + (r >> 61);
    if (r >= HASH_MOD){
        r -= HASH_MOD;
    }
    return r;
}

/**
 * @brief Memory-mapped binary representation of all the read sub-paths of the GAF.
 * File layout: D (2N uint32: X then R) | H (2N+1 uint64 prefix hashes of D) | offsets (S+1 uint64, start of each sub-path in X)
 */
class PathStore{
public:
    ~PathStore(){ release(); }

    bool build(const string& gaf_file, const string& bin_file, vector<Segment> &segments, unordered_map<string, int> &segments_IDs, vector<uint64_t>& occurrences);
    void release(){
        if (map != MAP_FAILED){
            munmap(map, map_len);
            map = MAP_FAILED;
        }
        D = nullptr; H = nullptr; off = nullptr;
        pw.clear(); pw.shrink_to_fit();
    }

    uint64_t nb_steps() const {return N;}
    uint64_t nb_subpaths() const {return S;}
    uint64_t subpath_begin(uint64_t sp) const {return off[sp];}
    uint32_t step(uint64_t i) const {return D[i];}

    //neighbor path starting just after position g of X, going forward
    PathRef forward_ref(uint64_t g, uint32_t len) const {return PathRef{g + 1, len};}
    //neighbor path made of X[g-1], X[g-2], ... reverse-complemented, i.e. R[N-g ...]
    PathRef reverse_ref(uint64_t g, uint32_t len) const {return PathRef{2 * N - g, len};}

    //length of the longest common prefix of D[a..a+maxlen) and D[b..b+maxlen)
    uint32_t lcp(uint64_t a, uint64_t b, uint32_t maxlen) const {
        if (a == b){
            return maxlen;
        }
        uint32_t k = 0;
        uint32_t lim = std::min<uint32_t>(maxlen, 8);
        for ( ; k < lim ; k++){
            if (D[a+k] != D[b+k]){
                return k;
            }
        }
        if (k == maxlen || range_hash(a, maxlen) == range_hash(b, maxlen)){
            return maxlen;
        }
        uint32_t lo = k, hi = maxlen; //prefix of length lo is equal, prefix of length hi is not
        while (hi - lo > 1){
            uint32_t mid = lo + (hi - lo) / 2;
            if (range_hash(a, mid) == range_hash(b, mid)){
                lo = mid;
            }
            else{
                hi = mid;
            }
        }
        while (lo < maxlen && D[a+lo] == D[b+lo]){ //only possible after a hash collision
            lo++;
        }
        return lo;
    }

    void materialize(const PathRef& r, vector<pair<int,bool>>& out) const {
        out.clear();
        out.reserve(r.len);
        for (uint64_t i = r.start ; i < r.start + r.len ; i++){
            out.push_back({(int) (D[i] >> 1), (bool) (D[i] & 1)});
        }
    }

private:
    uint64_t range_hash(uint64_t a, uint32_t len) const {
        uint64_t sub = mulmod(H[a], pw[len]);
        uint64_t h = H[a + len] + HASH_MOD - sub;
        if (h >= HASH_MOD){
            h -= HASH_MOD;
        }
        return h;
    }

    const uint32_t* D = nullptr;
    const uint64_t* H = nullptr;
    const uint64_t* off = nullptr;
    uint64_t N = 0; //number of steps in X
    uint64_t S = 0; //number of sub-paths
    void* map = MAP_FAILED;
    size_t map_len = 0;
    vector<uint64_t> pw; //powers of HASH_BASE
};

static bool write_all(int fd, const void* data, size_t bytes){
    const char* p = (const char*) data;
    while (bytes > 0){
        ssize_t w = ::write(fd, p, bytes);
        if (w < 0){
            if (errno == EINTR){
                continue;
            }
            return false;
        }
        p += w;
        bytes -= w;
    }
    return true;
}

/**
 * @brief Parse the GAF file once, cut the read paths where consecutive segments are not linked (exactly as before), and
 * write the sub-paths of length >= 2 (shorter ones have no neighbors) to bin_file, which is then memory-mapped read-only
 * and unlinked (the mapping stays valid until release()).
 * occurrences[s] is incremented for each occurrence of segment s in the stored sub-paths.
 */
bool PathStore::build(const string& gaf_file, const string& bin_file, vector<Segment> &segments, unordered_map<string, int> &segments_IDs, vector<uint64_t>& occurrences){

    int fd = ::open(bin_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd < 0){
        std::cerr << "ERROR: graphunzip cannot create temporary file " << bin_file << endl;
        return false;
    }

    vector<uint32_t> buffer;
    const size_t buffer_size = 1 << 16;
    buffer.reserve(buffer_size);
    bool write_ok = true;
    auto flush = [&](){
        if (!buffer.empty()){
            write_ok = write_ok && write_all(fd, buffer.data(), buffer.size() * sizeof(uint32_t));
            buffer.clear();
        }
    };

    vector<uint64_t> offsets = {0};
    uint64_t n = 0;
    uint32_t max_len = 0;

    ifstream gaf(gaf_file);
    string line;
    vector<pair<int,bool>> segments_now;
    auto store_subpath = [&](const vector<pair<int,bool>>& sub){
        if (sub.size() < 2){
            return;
        }
        for (const pair<int,bool>& st : sub){
            buffer.push_back(((uint32_t) st.first << 1) | (st.second ? 1u : 0u));
            occurrences[st.first]++;
            if (buffer.size() >= buffer_size){
                flush();
            }
        }
        n += sub.size();
        offsets.push_back(n);
        max_len = std::max<uint32_t>(max_len, sub.size());
    };

    while (getline(gaf, line)){

        //get the path as the 6th tab-delimited column (index 5), without going through a stringstream
        string path;
        size_t start_of_path = 0;
        bool enough_fields = true;
        for (int field = 0 ; field < 5 ; field++){
            start_of_path = line.find('\t', start_of_path);
            if (start_of_path == string::npos){
                enough_fields = false;
                break;
            }
            start_of_path++;
        }
        if (enough_fields){
            size_t end_of_path = line.find('\t', start_of_path);
            path = line.substr(start_of_path, (end_of_path == string::npos) ? string::npos : end_of_path - start_of_path);
        }

        bool orientation_now = true; //a path step without '<' or '>' (e.g. a stable-ID path with a single segment) is taken as forward
        string str_now;
        segments_now.clear();

        for (int c = 0 ; c < path.size(); c++){
            if (path[c] == '>' || path[c] == '<' || c == path.size()-1){

                if (c == path.size()-1){
                    str_now += path[c];
                }

                if (c != 0){
                    auto it_ID = segments_IDs.find(str_now);
                    if (it_ID == segments_IDs.end()){
                        cout << "Error: the segment " << str_now << " is in the gaf but cannot be found in the GFA" << endl;
                        cout << line << endl;
                        ::close(fd);
                        unlink(bin_file.c_str());
                        exit(1);
                    }
                    int new_contig_ID = it_ID->second;

                    //check that the segment can indeed be found next to the previous segment
                    bool found = false;
                    if (segments_now.size() > 0){
                        for (const pair<int,int>& link : segments[segments_now.back().first].links[segments_now.back().second].first){
                            if (link.first == new_contig_ID && link.second == !orientation_now){
                                found = true;
                            }
                        }
                    }
                    if (segments_now.size() == 0 || found){
                        segments_now.push_back({new_contig_ID, orientation_now});
                    }
                    else{ //then the neighbor is not the neighbor of a previous segment, cut
                        store_subpath(segments_now);
                        segments_now = {{new_contig_ID, orientation_now}};
                    }
                }

                if (path[c] == '>'){
                    orientation_now = true;
                }
                else{
                    orientation_now = false;
                }
                str_now.clear();
            }
            else{
                str_now += path[c];
            }
        }
        store_subpath(segments_now);
    }
    gaf.close();
    flush();

    N = n;
    S = offsets.size() - 1;
    map_len = (2 * N) * sizeof(uint32_t) + (2 * N + 1) * sizeof(uint64_t) + (S + 1) * sizeof(uint64_t);
    if (!write_ok || ftruncate(fd, map_len) != 0){
        std::cerr << "ERROR: graphunzip could not write temporary file " << bin_file << endl;
        ::close(fd);
        unlink(bin_file.c_str());
        return false;
    }

    //fill R, the hashes and the offsets in place
    void* wmap = mmap(nullptr, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (wmap == MAP_FAILED){
        std::cerr << "ERROR: graphunzip could not memory-map temporary file " << bin_file << endl;
        ::close(fd);
        unlink(bin_file.c_str());
        return false;
    }
    uint32_t* Dw = (uint32_t*) wmap;
    for (uint64_t j = 0 ; j < N ; j++){
        Dw[N + j] = Dw[N - 1 - j] ^ 1u;
    }
    uint64_t* Hw = (uint64_t*) ((char*) wmap + (2 * N) * sizeof(uint32_t));
    Hw[0] = 0;
    for (uint64_t i = 0 ; i < 2 * N ; i++){
        uint64_t h = mulmod(Hw[i], HASH_BASE) + (uint64_t) Dw[i] + 1;
        if (h >= HASH_MOD){
            h -= HASH_MOD;
        }
        Hw[i+1] = h;
    }
    uint64_t* offw = Hw + (2 * N + 1);
    std::copy(offsets.begin(), offsets.end(), offw);
    offsets.clear();
    offsets.shrink_to_fit();
    munmap(wmap, map_len);

    map = mmap(nullptr, map_len, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd);
    unlink(bin_file.c_str()); //the file stays readable through the mapping and disappears when it is unmapped
    if (map == MAP_FAILED){
        std::cerr << "ERROR: graphunzip could not memory-map temporary file " << bin_file << endl;
        return false;
    }
    D = (const uint32_t*) map;
    H = (const uint64_t*) ((const char*) map + (2 * N) * sizeof(uint32_t));
    off = H + (2 * N + 1);

    pw.resize((size_t) max_len + 1);
    pw[0] = 1;
    for (size_t i = 1 ; i < pw.size() ; i++){
        pw[i] = mulmod(pw[i-1], HASH_BASE);
    }
    return true;
}

const uint32_t NO_NODE = UINT32_MAX;

struct TrieNode{
    uint64_t lab;   //the edge from the parent is D[lab .. lab+len)
    uint32_t len;
    uint32_t cnt;   //number of neighbor paths having the path root->this node as prefix
    uint32_t end;   //number of neighbor paths ending exactly here
    uint32_t depth; //length of the path root->this node
    uint32_t child;
    uint32_t sib;
};

struct EndScratch{
    vector<PathRef> refs;
    vector<TrieNode> nodes;
    vector<pair<uint64_t,uint32_t>> maximal; //start of the path in D, leaf node
    vector<uint32_t> stack;
    vector<char> cons_ok;
    vector<uint32_t> cons_len;
    vector<uint32_t> strong_len;
    vector<uint32_t> strong_end; //where the strong prefix ends in the trie: 2*node if exactly at node, 2*node+1 if one step into the edge leading to node
    vector<char> strong_seen;
};

inline bool consensual(int first, int second){
    return 0.2*(first-1) > second-1 && second < 5;
}

inline bool strong(int first, int second, int min_coverage){
    return first >= min_coverage && (second == 0 || first > 4);
}

/**
 * @brief Compute the consensus and the strong neighbors of one segment end from its neighbor paths (sc.refs, in GAF order).
 * @return whether this end is haploid
 */
bool process_end(const PathStore& ps, EndScratch& sc, int min_coverage, PathRef& consensus, vector<PathRef>& strong_neighbors){

    consensus = PathRef();
    strong_neighbors.clear();
    vector<PathRef>& refs = sc.refs;
    if (refs.empty()){
        return true;
    }

    //same comparator as before on the same sequence of lengths => same permutation
    std::sort(refs.begin(), refs.end(), [](const PathRef& a, const PathRef& b){return a.len > b.len;});

    vector<TrieNode>& nodes = sc.nodes;
    nodes.clear();
    nodes.push_back(TrieNode{0, 0, 0, 0, 0, NO_NODE, NO_NODE});
    sc.maximal.clear();

    for (const PathRef& r : refs){
        uint32_t node = 0;
        uint32_t pos = 0;
        nodes[0].cnt++;
        while (true){
            if (pos == r.len){ //r is a prefix of a path already inserted
                nodes[node].end++;
                break;
            }
            uint32_t first_step = ps.step(r.start + pos);
            uint32_t prev = NO_NODE;
            uint32_t c = nodes[node].child;
            while (c != NO_NODE && ps.step(nodes[c].lab) != first_step){
                prev = c;
                c = nodes[c].sib;
            }
            if (c == NO_NODE){ //new branch: r is a new maximal path
                uint32_t leaf = nodes.size();
                nodes.push_back(TrieNode{r.start + pos, r.len - pos, 1, 1, nodes[node].depth + (r.len - pos), NO_NODE, nodes[node].child});
                nodes[node].child = leaf;
                sc.maximal.push_back({r.start, leaf});
                break;
            }
            uint32_t edge_len = nodes[c].len;
            uint32_t m = ps.lcp(nodes[c].lab, r.start + pos, std::min(edge_len, r.len - pos));
            if (m == edge_len){
                nodes[c].cnt++;
                node = c;
                pos += m;
                continue;
            }
            //split the edge after m steps
            uint32_t mid = nodes.size();
            nodes.push_back(TrieNode{nodes[c].lab, m, nodes[c].cnt + 1, 0, nodes[node].depth + m, c, nodes[c].sib});
            if (prev == NO_NODE){
                nodes[node].child = mid;
            }
            else{
                nodes[prev].sib = mid;
            }
            nodes[c].lab += m;
            nodes[c].len -= m;
            nodes[c].sib = NO_NODE;
            pos += m;
            if (pos == r.len){
                nodes[mid].end = 1;
                break;
            }
            uint32_t leaf = nodes.size();
            nodes.push_back(TrieNode{r.start + pos, r.len - pos, 1, 1, nodes[mid].depth + (r.len - pos), NO_NODE, c});
            nodes[mid].child = leaf;
            sc.maximal.push_back({r.start, leaf});
            break;
        }
    }

    //propagate the per-position strengths from the root. On an edge u->w of length e, the first position has
    //(first, second) = (cnt(w), cnt(u)-end(u)-cnt(w)), the other positions (cnt(w), 0)
    sc.cons_ok.assign(nodes.size(), 0);
    sc.cons_len.assign(nodes.size(), 0);
    sc.strong_len.assign(nodes.size(), 0);
    sc.strong_end.assign(nodes.size(), 0);
    sc.cons_ok[0] = 1;
    sc.stack.clear();
    sc.stack.push_back(0);
    while (!sc.stack.empty()){
        uint32_t u = sc.stack.back();
        sc.stack.pop_back();
        for (uint32_t w = nodes[u].child ; w != NO_NODE ; w = nodes[w].sib){
            int C = nodes[w].cnt;
            int X = (int) nodes[u].cnt - (int) nodes[u].end - C;
            uint32_t e = nodes[w].len;
            sc.cons_ok[w] = sc.cons_ok[u] && consensual(C, X) && (e == 1 || consensual(C, 0));
            sc.cons_len[w] = sc.cons_len[u] + (C > 1 ? e : 0);
            if (sc.strong_len[u] < nodes[u].depth){
                sc.strong_len[w] = sc.strong_len[u];
                sc.strong_end[w] = sc.strong_end[u];
            }
            else if (!strong(C, X, min_coverage)){
                sc.strong_len[w] = nodes[u].depth;
                sc.strong_end[w] = 2 * u;
            }
            else if (e == 1 || strong(C, 0, min_coverage)){
                sc.strong_len[w] = nodes[w].depth;
                sc.strong_end[w] = 2 * w;
            }
            else{
                sc.strong_len[w] = nodes[u].depth + 1;
                sc.strong_end[w] = 2 * w + 1;
            }
            sc.stack.push_back(w);
        }
    }

    bool haploid = false;
    for (const pair<uint64_t,uint32_t>& p : sc.maximal){
        if (sc.cons_ok[p.second]){
            haploid = true;
            consensus = PathRef{p.first, sc.cons_len[p.second]};
            break;
        }
    }
    //strong neighbors, in the order of the maximal paths. Several maximal paths often share the same strong prefix: only the
    //first occurrence is kept, since list_non_represented_paths / UnrepresentedPathsAdder give the same result on a repeated path
    sc.strong_seen.assign(2 * nodes.size(), 0);
    for (const pair<uint64_t,uint32_t>& p : sc.maximal){
        if (sc.strong_len[p.second] > 0 && !sc.strong_seen[sc.strong_end[p.second]]){
            sc.strong_seen[sc.strong_end[p.second]] = 1;
            strong_neighbors.push_back(PathRef{p.first, sc.strong_len[p.second]});
        }
    }
    return haploid;
}

/**
 * @brief Everything graphunzip needs to know about the reads: which segments are haploid, their consensus left/right and
 * their strong neighbor paths, all stored as references into the memory-mapped PathStore.
 */
struct ReadPathsInfo{
    PathStore store;
    vector<char> haploid;
    vector<PathRef> consensus; //index 2*segment + end (0 = left, 1 = right)
    vector<vector<PathRef>> strong_neighbors; //index 2*segment + end

    bool is_haploid(size_t s) const {return haploid[s];}

    pair<int,bool> step_at(uint64_t i) const {
        uint32_t st = store.step(i);
        return {(int) (st >> 1), (bool) (st & 1)};
    }

    const vector<PathRef>& strong_neighbor_refs(size_t s, int end) const {return strong_neighbors[2*s+end];}

    vector<pair<int,bool>> get_consensus(size_t s, int end) const {
        vector<pair<int,bool>> res;
        store.materialize(consensus[2*s+end], res);
        return res;
    }

    vector<vector<pair<int,bool>>> get_strong_neighbors(size_t s, int end) const {
        vector<vector<pair<int,bool>>> res(strong_neighbors[2*s+end].size());
        for (size_t i = 0 ; i < res.size() ; i++){
            store.materialize(strong_neighbors[2*s+end][i], res[i]);
        }
        return res;
    }

    //free everything except the haploid flags
    void release_paths(){
        store.release();
        consensus.clear(); consensus.shrink_to_fit();
        strong_neighbors.clear(); strong_neighbors.shrink_to_fit();
    }

    void compute(size_t nb_segments, int min_coverage, int threads, const vector<uint64_t>& occurrences, uint64_t memory_budget);
};

/**
 * @brief Compute haploidy, consensuses and strong neighbors of all segments, by batches of segments whose occurrence index
 * fits in memory_budget bytes. Each batch scans the memory-mapped sub-paths once.
 */
void ReadPathsInfo::compute(size_t nb_segments, int min_coverage, int threads, const vector<uint64_t>& occurrences, uint64_t memory_budget){

    struct Occurrence{
        uint64_t g;          //position in X
        uint32_t left_len;   //number of steps of the sub-path before g
        uint32_t right_len;  //number of steps of the sub-path after g
    };

    haploid.assign(nb_segments, 0);
    consensus.assign(2 * nb_segments, PathRef());
    strong_neighbors.assign(2 * nb_segments, {});

    const uint64_t S = store.nb_subpaths();

    size_t lo = 0;
    while (lo < nb_segments){
        //choose the batch [lo, hi)
        size_t hi = lo;
        uint64_t nb_occ = 0;
        while (hi < nb_segments && (hi == lo || (nb_occ + occurrences[hi]) * sizeof(Occurrence) <= memory_budget)){
            nb_occ += occurrences[hi];
            hi++;
        }

        //index the occurrences of the segments of the batch, in GAF order
        vector<uint64_t> first_occ(hi - lo + 1, 0);
        for (size_t s = lo ; s < hi ; s++){
            first_occ[s - lo + 1] = first_occ[s - lo] + occurrences[s];
        }
        vector<Occurrence> occs(nb_occ);
        {
            vector<uint64_t> fill(first_occ.begin(), first_occ.end() - 1);
            for (uint64_t sp = 0 ; sp < S ; sp++){
                uint64_t a = store.subpath_begin(sp);
                uint64_t b = store.subpath_begin(sp + 1);
                for (uint64_t g = a ; g < b ; g++){
                    uint32_t id = store.step(g) >> 1;
                    if (id >= lo && id < hi){
                        occs[fill[id - lo]++] = Occurrence{g, (uint32_t) (g - a), (uint32_t) (b - g - 1)};
                    }
                }
            }
        }

        #pragma omp parallel num_threads(std::max(1, threads))
        {
            EndScratch sc;
            #pragma omp for schedule(dynamic, 16)
            for (size_t s = lo ; s < hi ; s++){
                bool haploid_ends[2];
                for (int end = 0 ; end < 2 ; end++){
                    sc.refs.clear();
                    for (uint64_t o = first_occ[s - lo] ; o < first_occ[s - lo + 1] ; o++){
                        const Occurrence& occ = occs[o];
                        bool forward = store.step(occ.g) & 1;
                        //right of a forward occurrence / left of a reverse occurrence: the rest of the sub-path
                        PathRef r = (forward == (end == 1)) ? store.forward_ref(occ.g, occ.right_len) : store.reverse_ref(occ.g, occ.left_len);
                        if (r.len > 0){
                            sc.refs.push_back(r);
                        }
                    }
                    haploid_ends[end] = process_end(store, sc, min_coverage, consensus[2*s+end], strong_neighbors[2*s+end]);
                }
                haploid[s] = haploid_ends[0] && haploid_ends[1];
            }
        }
        lo = hi;
    }
}

} //namespace


/**
 * @brief Unzip the graph and create a new list of segments. 
 * 
 * @param old_segments 
 * @param new_segments 
 * @param min_coverage 
 * @param unordered_map<int, vector<set<int>>> already_built_bridges; //associates to each haploid contig two sets, one of the contigs it is connected to on the left, one on the right

 */
void create_haploid_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, unordered_map<int, vector<set<int>>> &already_built_bridges, int min_coverage, bool contiguity, const ReadPathsInfo& paths){

    new_segments = {};

    //associates all old contigs ids to their new IDs (can be multiple if contig is duplicated)
    int new_ID = 0;


    //begin by building bridges between haploid contigs
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){

        
        if (paths.is_haploid(old_segment)){


            //first create the contig if needed
            if (old_ID_to_new_IDs.find(old_segment) == old_ID_to_new_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_segment].name + "_0" , new_ID, old_segments[old_segment].get_pos_in_file(), old_segments[old_segment].get_length(), old_segments[old_segment].get_coverage()));
                old_ID_to_new_IDs[old_segment] = {new_ID};
                new_ID++;
            }
            int ID_of_new_contig = old_ID_to_new_IDs[old_segment][0];

            //now see if it already has neighbors left and right and create them until the next haploid contig if not

            //right
            vector<pair<int,bool>> cons_right = paths.get_consensus(old_segment, 1);
            bool there_is_a_bridge = false;
            double coverage_bridge = old_segments[old_segment].get_coverage();
            int otherEndOfBridge = 0;
            int endOfOtherEndOfBridge = 0;
            for (auto contig_and_orientation : cons_right){
                if (paths.is_haploid(contig_and_orientation.first)){
                    //in case contiguity mode is on, check that the bridge is reciprocal
                    bool reciprocal = false;
                    if (contig_and_orientation.second == 1 && contiguity){
                        for (auto neighbor : paths.get_consensus(contig_and_orientation.first, 0)){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }
                    if (contig_and_orientation.second == 0 && contiguity){
                        for (auto neighbor : paths.get_consensus(contig_and_orientation.first, 1)){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }

                    if (!contiguity || reciprocal){
                        there_is_a_bridge = true;
                        coverage_bridge = std::min(old_segments[contig_and_orientation.first].get_original_coverage() , old_segments[old_segment].get_original_coverage());
                        otherEndOfBridge = contig_and_orientation.first;
                        endOfOtherEndOfBridge = contig_and_orientation.second;
                    }
                    break;
                }
            }
            if(there_is_a_bridge && already_built_bridges[old_segment][1].find(otherEndOfBridge) == already_built_bridges[old_segment][1].end()){

                int contig_in_cons = 0;
                int previous_ID = ID_of_new_contig;
                int previous_old_ID = old_segment;
                int previous_end = 1;
                while (true){

                    int end_of_contig_right = 0;
                    if (!cons_right[contig_in_cons].second){
                        end_of_contig_right = 1;
                    }

                    //create the contig
                    if (!paths.is_haploid(cons_right[contig_in_cons].first)){
                        double coverage_of_this_contig = std::min(coverage_bridge, old_segments[cons_right[contig_in_cons].first].get_coverage());
                        new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_right[contig_in_cons].first].size()) , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), old_segments[cons_right[contig_in_cons].first].get_length(), coverage_of_this_contig)); //each copy keeps its own length (not the length of the whole bridge)
                        if (old_ID_to_new_IDs.find(cons_right[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                            old_ID_to_new_IDs[cons_right[contig_in_cons].first] = {};
                        }
                        old_ID_to_new_IDs[cons_right[contig_in_cons].first].push_back(new_ID);
                        old_segments[cons_right[contig_in_cons].first].decrease_coverage(coverage_of_this_contig);
                        new_ID++;
                    }
                    else{
                        //create the haploid contig only if not already created
                        if (old_ID_to_new_IDs.find(cons_right[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                            new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), old_segments[cons_right[contig_in_cons].first].get_length(), old_segments[cons_right[contig_in_cons].first].get_coverage()));
                            old_ID_to_new_IDs[cons_right[contig_in_cons].first] = {new_ID};
                            new_ID++;
                        }
                    }

                    int ID_of_new_contig_right = old_ID_to_new_IDs[cons_right[contig_in_cons].first][old_ID_to_new_IDs[cons_right[contig_in_cons].first].size() - 1];

                    //add the link
                    //find the CIGAR in the old segment
                    string cigar = "*";
                    int idx = 0;
                    for (auto link : old_segments[cons_right[contig_in_cons].first].links[end_of_contig_right].first){
                        if (link.first == previous_old_ID && link.second == previous_end){
                            cigar = old_segments[cons_right[contig_in_cons].first].links[end_of_contig_right].second[idx];
                        }
                        idx++;
                    }
                    new_segments[previous_ID].links[previous_end].first.push_back({ID_of_new_contig_right, end_of_contig_right});
                    new_segments[previous_ID].links[previous_end].second.push_back(cigar);
                    new_segments[ID_of_new_contig_right].links[end_of_contig_right].first.push_back({previous_ID, previous_end});
                    new_segments[ID_of_new_contig_right].links[end_of_contig_right].second.push_back(cigar);

                    //go to the next contig
                    if (paths.is_haploid(cons_right[contig_in_cons].first)){
                        already_built_bridges[cons_right[contig_in_cons].first][end_of_contig_right].insert(old_segment);
                        already_built_bridges[old_segment][1].insert(cons_right[contig_in_cons].first);
                        break;
                    }
                    previous_end = 1-end_of_contig_right;
                    previous_ID = ID_of_new_contig_right;
                    previous_old_ID = cons_right[contig_in_cons].first;
                    contig_in_cons++;
                }
            }   

            //left
            vector<pair<int,bool>> cons_left = paths.get_consensus(old_segment, 0);
            there_is_a_bridge = false;
            coverage_bridge = old_segments[old_segment].get_coverage();
            for (auto contig_and_orientation : cons_left){
                if (paths.is_haploid(contig_and_orientation.first)){
                    bool reciprocal = false;

                    if (contig_and_orientation.second == 1 && contiguity){
                        for (auto neighbor : paths.get_consensus(contig_and_orientation.first, 0)){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }
                    if (contig_and_orientation.second == 0 && contiguity){
                        for (auto neighbor : paths.get_consensus(contig_and_orientation.first, 1)){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }

                    if (!contiguity || reciprocal){
                        there_is_a_bridge = true;
                        coverage_bridge = std::min(old_segments[contig_and_orientation.first].get_original_coverage() , old_segments[old_segment].get_original_coverage());
                        otherEndOfBridge = contig_and_orientation.first;
                        endOfOtherEndOfBridge = contig_and_orientation.second;
                    }
                    break;
                }
            }

            if(there_is_a_bridge && already_built_bridges[old_segment][0].find(otherEndOfBridge) == already_built_bridges[old_segment][0].end()){

                int contig_in_cons = 0;
                int previous_ID = ID_of_new_contig;
                int previous_old_ID = old_segment;
                int previous_end = 0;
                while (true){

                    int end_of_contig_left = 0;
                    if (!cons_left[contig_in_cons].second){
                        end_of_contig_left = 1;
                    }

                    //create the contig
                    if (!paths.is_haploid(cons_left[contig_in_cons].first)){
                        double coverage_of_this_contig = std::min(coverage_bridge, old_segments[cons_left[contig_in_cons].first].get_coverage());
                        new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_left[contig_in_cons].first].size()) , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), old_segments[cons_left[contig_in_cons].first].get_length(), coverage_of_this_contig)); //each copy keeps its own length (not the length of the whole bridge)
                        if (old_ID_to_new_IDs.find(cons_left[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                            old_ID_to_new_IDs[cons_left[contig_in_cons].first] = {};
                        }
                        old_ID_to_new_IDs[cons_left[contig_in_cons].first].push_back(new_ID);
                        old_segments[cons_left[contig_in_cons].first].decrease_coverage(coverage_of_this_contig);
                        new_ID++;
                    }
                    else{
                        //create the haploid contig only if not already created
                        if (old_ID_to_new_IDs.find(cons_left[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                            new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), old_segments[cons_left[contig_in_cons].first].get_length(), old_segments[cons_left[contig_in_cons].first].get_coverage()));
                            old_ID_to_new_IDs[cons_left[contig_in_cons].first] = {new_ID};
                            new_ID++;
                        }
                    }

                    int ID_of_new_contig_left = old_ID_to_new_IDs[cons_left[contig_in_cons].first][old_ID_to_new_IDs[cons_left[contig_in_cons].first].size() - 1];

                    //add the link
                    //find the CIGAR in the old segment
                    string cigar = "*";
                    int idx = 0;
                    for (auto link : old_segments[cons_left[contig_in_cons].first].links[end_of_contig_left].first){
                        if (link.first == previous_old_ID && link.second == previous_end){
                            cigar = old_segments[cons_left[contig_in_cons].first].links[end_of_contig_left].second[idx];
                        }
                        idx++;
                    }

                    new_segments[previous_ID].links[previous_end].first.push_back({ID_of_new_contig_left, end_of_contig_left});
                    new_segments[previous_ID].links[previous_end].second.push_back(cigar);
                    new_segments[ID_of_new_contig_left].links[end_of_contig_left].first.push_back({previous_ID, previous_end});
                    new_segments[ID_of_new_contig_left].links[end_of_contig_left].second.push_back(cigar);

                    //go to the next contig
                    if (paths.is_haploid(cons_left[contig_in_cons].first)){
                        already_built_bridges[cons_left[contig_in_cons].first][end_of_contig_left].insert(old_segment);
                        already_built_bridges[old_segment][0].insert(cons_left[contig_in_cons].first);
                        break;
                    }
                    previous_end = 1-end_of_contig_left;
                    previous_ID = ID_of_new_contig_left;
                    previous_old_ID = cons_left[contig_in_cons].first;
                    contig_in_cons++;
                }
            }        
        }
    }
}

/**
 * @brief Take the non represented paths of the graph one by one, in the order in which list_non_represented_paths finds them,
 * and add contigs and links until all paths are represented. The paths are consumed as soon as they are found instead of
 * being all stored first (same result, since list_non_represented_paths does not read what this modifies).
 */
class UnrepresentedPathsAdder{
public:
    UnrepresentedPathsAdder(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, int min_coverage, const ReadPathsInfo& paths)
        : old_segments(old_segments), new_segments(new_segments), old_ID_to_new_IDs(old_ID_to_new_IDs), min_coverage(min_coverage){
        for (size_t s_idx = 0 ; s_idx < old_segments.size() ; s_idx++){ //we are not going to create new versions of haploid contigs
            Segment& s = old_segments[s_idx];
            if (paths.is_haploid(s_idx)){
                old_IDs_to_new_non_haploid_IDs[s.ID] = old_ID_to_new_IDs[s.ID][0];
            }
        }
    }

    //convert an unrepresented path in a list of links that must be there in the final graph
    void add(const vector<pair<int,bool>>& path){
        //compute the coverage of the path
        double coverage = old_segments[path[0].first].get_coverage();
        for (const pair<int,bool>& contig : path){
            coverage = std::min(coverage, old_segments[contig.first].get_coverage());
        }
        //if the coverage is too low, we don't add the path
        if (coverage < min_coverage){
            return;
        }

        for (int contig = 0 ; contig < path.size() - 1 ; contig++){

            int old_ID1 = path[contig].first;
            int old_ID2 = path[contig+1].first;

            int end1 = 1;
            int end2 = 0;
            if (!path[contig].second){
                end1 = 0;
            }
            if (!path[contig+1].second){
                end2 = 1;
            }

            //everything below only depends on (old_ID1, end1, old_ID2, end2) and the state it leaves, and doing it again is a
            //no-op (contigs are created once, links_to_add is a set): skip the links already processed, which are seen again and again
            uint64_t link_key = ((uint64_t) (uint32_t) old_ID1 << 33) | ((uint64_t) end1 << 32) | ((uint64_t) (uint32_t) old_ID2 << 1) | (uint64_t) end2;
            if (!processed_links.insert(link_key).second){
                continue;
            }

            //create the contigs if not already done
            if (old_IDs_to_new_non_haploid_IDs.find(old_ID1) == old_IDs_to_new_non_haploid_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_ID1].name + "_" + std::to_string(old_ID_to_new_IDs[old_ID1].size()) , new_segments.size(), old_segments[old_ID1].get_pos_in_file(), old_segments[old_ID1].get_length(), old_segments[old_ID1].get_coverage()));
                old_IDs_to_new_non_haploid_IDs[old_ID1] = new_segments.size() - 1;
                if (old_ID_to_new_IDs.find(old_ID1) == old_ID_to_new_IDs.end()){
                    old_ID_to_new_IDs[old_ID1] = {(int) new_segments.size() - 1};
                }
                else{
                    old_ID_to_new_IDs[old_ID1].push_back(new_segments.size() - 1);
                }
            }

            if (old_IDs_to_new_non_haploid_IDs.find(old_ID2) == old_IDs_to_new_non_haploid_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_ID2].name + "_" + std::to_string(old_ID_to_new_IDs[old_ID2].size()) , new_segments.size(), old_segments[old_ID2].get_pos_in_file(), old_segments[old_ID2].get_length(), old_segments[old_ID2].get_coverage()));
                old_IDs_to_new_non_haploid_IDs[old_ID2] = new_segments.size() - 1;
                if (old_ID_to_new_IDs.find(old_ID2) == old_ID_to_new_IDs.end()){
                    old_ID_to_new_IDs[old_ID2] = {(int) new_segments.size() - 1};
                }
                else{
                    old_ID_to_new_IDs[old_ID2].push_back(new_segments.size() - 1);
                }
            }

            int new_ID1 = old_ID_to_new_IDs[old_ID1][old_ID_to_new_IDs[old_ID1].size()-1];
            int new_ID2 = old_ID_to_new_IDs[old_ID2][old_ID_to_new_IDs[old_ID2].size()-1];
            string cigar = "*";
            int idx = 0;
            for (const pair<int,int>& link : old_segments[old_ID2].links[end2].first){
                if (link.first == old_ID1 && link.second == end1){
                    cigar = old_segments[old_ID2].links[end2].second[idx];
                }
                idx ++;
            }


            if (links_to_add.find({{pair<int,int>(new_ID2, end2), pair<int,int>(new_ID1, end1)}, cigar}) == links_to_add.end()){
                links_to_add.insert({{pair<int,int>(new_ID1, end1), pair<int,int>(new_ID2, end2)}, cigar});
            }
        }
    }

    //add the links
    void finish(){
        for (const pair<pair<pair<int,int>, pair<int,int>>, string>& link : links_to_add){
            new_segments[link.first.first.first].links[link.first.first.second].first.push_back({link.first.second.first, link.first.second.second});
            new_segments[link.first.first.first].links[link.first.first.second].second.push_back(link.second);
            new_segments[link.first.second.first].links[link.first.second.second].first.push_back({link.first.first.first, link.first.first.second});
            new_segments[link.first.second.first].links[link.first.second.second].second.push_back(link.second);
        }
    }

private:
    vector<Segment> &old_segments;
    vector<Segment> &new_segments;
    unordered_map<int, std::vector<int>>& old_ID_to_new_IDs;
    int min_coverage;
    unordered_map<int, int> old_IDs_to_new_non_haploid_IDs; //associates old IDs to new IDs for the contig we are going to create
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add;
    robin_hood::unordered_flat_set<uint64_t> processed_links;
};

/**
 * @brief Check what is not seen yet in new_segments
 * 
 * @param old_segments 
 * @param new_segments 
 * @param old_ID_to_new_IDs 
 * @param min_coverage 
 * @param all_paths
 * @return vector<vector<pair<int,bool>>> all the non represented paths
 */
void list_non_represented_paths(vector<Segment> &old_segments, unordered_map<int, vector<set<int>>> &already_built_bridges, int min_coverage, const ReadPathsInfo& paths, UnrepresentedPathsAdder& unrepresented_paths){

    //first index represented paths
    vector<vector<pair<int,bool>>> represented_paths;
    unordered_map<int, vector<pair<int,int>>> where_is_this_contig_represented; //associates to an ID (indices of the path, position in the path)

    //go through all the haploids segments and list the represented paths left and right in new_segments
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){
        if (paths.is_haploid(old_segment)){
            
            std::vector<std::pair<int,bool>> consensus_left = {{old_segments[old_segment].ID, false}};
            std::vector<std::pair<int,bool>> new_elements = paths.get_consensus(old_segment, 0);
            consensus_left.insert(consensus_left.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            bool there_is_a_bridge = false;
            int idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_left){
                if (paths.is_haploid(contig_and_orientation.first) && idx_of_last_haploid_contig != 0){
                    if (already_built_bridges[old_segment][0].find(contig_and_orientation.first) != already_built_bridges[old_segment][0].end()){
                        there_is_a_bridge = true;
                    }
                    break;
                }
                idx_of_last_haploid_contig++;
            }
            if (there_is_a_bridge){
                auto represented_path = vector<pair<int,bool>>(consensus_left.begin(), consensus_left.begin() + idx_of_last_haploid_contig + 1);
                represented_paths.push_back(represented_path);
                std::reverse(represented_path.begin(), represented_path.end());
                //reverse the orientations
                for (auto &contig : represented_path){
                    contig.second = !contig.second;
                }
                represented_paths.push_back(represented_path);

                //fill where_is_this_contig_represented
                int number_of_haploid_contigs = 0;
                int idx = 0;
                for (auto contig_and_orientation : consensus_left){

                    if (where_is_this_contig_represented.find(contig_and_orientation.first) == where_is_this_contig_represented.end()){
                        where_is_this_contig_represented[contig_and_orientation.first] = {};
                    }
                    int symmetrical_idx = represented_paths[represented_paths.size() - 1].size() - 1 - idx;
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 1, symmetrical_idx});
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 2, idx});

                    if (paths.is_haploid(contig_and_orientation.first)){
                        number_of_haploid_contigs++;
                        if (number_of_haploid_contigs > 1){
                            break;
                        }
                    }
                    idx++;
                }
            }

            vector<pair<int,bool>> consensus_right = {{old_segments[old_segment].ID, true}};
            new_elements = paths.get_consensus(old_segment, 1);
            consensus_right.insert(consensus_right.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            there_is_a_bridge = false;
            idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_right){
                if (paths.is_haploid(contig_and_orientation.first) && idx_of_last_haploid_contig != 0){
                    if (already_built_bridges[old_segment][1].find(contig_and_orientation.first) != already_built_bridges[old_segment][1].end()){
                        there_is_a_bridge = true;
                    }
                    break;
                }
                idx_of_last_haploid_contig++;
            }

            if (there_is_a_bridge){

                auto represented_path = vector<pair<int,bool>>(consensus_right.begin(), consensus_right.begin() + idx_of_last_haploid_contig + 1);
                represented_paths.push_back(represented_path);
                std::reverse(represented_path.begin(), represented_path.end());
                //reverse the orientations
                for (auto &contig : represented_path){
                    contig.second = !contig.second;
                }
                represented_paths.push_back(represented_path);

                //fill where_is_this_contig_represented
                int number_of_haploid_contigs = 0;
                int idx = 0;
                for (auto contig_and_orientation : consensus_right){

                    if (where_is_this_contig_represented.find(contig_and_orientation.first) == where_is_this_contig_represented.end()){
                        where_is_this_contig_represented[contig_and_orientation.first] = {};
                    }
                    int symmetrical_idx = represented_paths[represented_paths.size() - 1].size() - 1 - idx;
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 2, idx});
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 1, symmetrical_idx});

                    if (paths.is_haploid(contig_and_orientation.first)){
                        number_of_haploid_contigs++;
                        if (number_of_haploid_contigs > 1){
                            break;
                        }
                    }
                    idx++;
                }
            }
        }
    }


    //now go through all the paths and see if they are represented
    //A path of a segment s is cut into pieces ending at each haploid contig: the first piece starts with s, the next ones start
    //with the haploid contig that ended the previous piece. Pieces are kept as (head, range of the strong neighbor path) and only
    //materialised when they are handed to unrepresented_paths
    struct Piece{
        pair<int,bool> head;
        uint64_t begin; //range [begin, end) of D following the head
        uint64_t end;
    };
    vector<Piece> unrepresented_paths_tmp; //inventory missing paths, but wait to have a path left and right before concatenating this to the unrepresented paths
    vector<pair<int,bool>> piece_buffer;

    //is head + D[begin..end) equal to the represented path starting at it
    auto piece_equal = [&paths](const Piece& piece, vector<pair<int,bool>>::const_iterator it){
        if (*it != piece.head){
            return false;
        }
        ++it;
        for (uint64_t i = piece.begin ; i < piece.end ; i++, ++it){
            if (*it != paths.step_at(i)){
                return false;
            }
        }
        return true;
    };

    for (size_t s_idx = 0 ; s_idx < old_segments.size() ; s_idx++){
        Segment& s = old_segments[s_idx];

        unrepresented_paths_tmp.clear();
        bool unrepresented_left = false;
        bool unrepresented_right = false;
        for (int end_of_s = 0 ; end_of_s < 2 ; end_of_s++){ //left paths first, then right paths
            for (const PathRef& path : paths.strong_neighbor_refs(s_idx, end_of_s)){

                //which side of s this path is on. Must be tracked explicitly: after the first haploid contig H, the pieces
                //start with H, whose orientation says nothing about the side of s
                bool path_is_on_the_right = (end_of_s == 1);

                //check if the path is represented
                Piece piece{{s.ID, path_is_on_the_right}, path.start, path.start};
                for (uint64_t i = path.start ; i < path.start + path.len ; i++){
                    pair<int,bool> contig = paths.step_at(i);
                    piece.end = i + 1;
                    if (paths.is_haploid(contig.first)){
                        size_t piece_size = 1 + (piece.end - piece.begin);
                        bool found = false;
                        auto it_where = where_is_this_contig_represented.find(contig.first);
                        if (it_where != where_is_this_contig_represented.end()){
                            for (const pair<int,int>& path_and_pos : it_where->second){
                                //compare the piece with the beginning (if pos==0) or the end of the represented path
                                const vector<pair<int,bool>>& represented = represented_paths[path_and_pos.first];
                                if (piece_size > represented.size()){
                                    continue;
                                }
                                if (path_and_pos.second == 0){
                                    found = piece_equal(piece, represented.begin());
                                }
                                else {
                                    found = piece_equal(piece, represented.end() - piece_size);
                                }
                                if (found){
                                    break;
                                }
                            }
                        }
                        if (!found){
                            unrepresented_paths_tmp.push_back(piece);
                            if (path_is_on_the_right){
                                unrepresented_right = true;
                            }
                            else {
                                unrepresented_left = true;
                            }
                        }
                        piece = Piece{contig, i + 1, i + 1};
                    }
                }
                //check the last piece
                {
                    size_t piece_size = 1 + (piece.end - piece.begin);
                    bool found = false;

                    auto it_where = where_is_this_contig_represented.find(piece.head.first);
                    if (it_where != where_is_this_contig_represented.end()){
                        for (const pair<int,int>& path_and_pos : it_where->second){
                            const vector<pair<int,bool>>& represented = represented_paths[path_and_pos.first];
                            if (represented.size() >= path_and_pos.second + piece_size
                                && piece_equal(piece, represented.begin() + path_and_pos.second)){
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found){
                        unrepresented_paths_tmp.push_back(piece);
                        if (path_is_on_the_right){
                            unrepresented_right = true;
                        }
                        else {
                            unrepresented_left = true;
                        }
                    }
                }

                if (unrepresented_left && unrepresented_right) { //if there is no path either left or right, it looks like this contig is not represented in the haploid contigs right. If only left or right, this kind of looks like an erroneous path
                    for (const Piece& p : unrepresented_paths_tmp){
                        piece_buffer.clear();
                        piece_buffer.push_back(p.head);
                        for (uint64_t i = p.begin ; i < p.end ; i++){
                            piece_buffer.push_back(paths.step_at(i));
                        }
                        unrepresented_paths.add(piece_buffer);
                    }
                    unrepresented_paths_tmp.clear();
                }
            }
        }
    }
}

/**
 * @brief Check that a graph file written by output_graph exists (and is not empty if it should contain something).
 * output_graph does not report errors, so this is the only way to detect e.g. an unwritable path.
 */
static bool graph_file_written(const string& path, const vector<Segment>& segments){
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f){
        return false;
    }
    bool expect_content = false;
    for (const Segment& s : segments){
        if (s.name != "delete_me"){
            expect_content = true;
            break;
        }
    }
    return !expect_content || f.tellg() > 0;
}

int main(int argc, char *argv[])
{

    if (argc != 11){
        //if -h or --help is passed as an argument, print the help
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)){
            std::cout << "Usage: graphunzip <gfa_input> <gaf_file> <min_coverage> <threads> <rename> <gfa_output> <contiguity> <single_genome> <last_k_used> <logfile>" << std::endl;
            return 0;
        }
        std::cout << "Usage: graphunzip <gfa_input> <gaf_file> <min_coverage> <threads> <rename> <gfa_output> <contiguity> <single_genome> <last_k_used> <logfile>" << std::endl;
        return 1;
    }

    std::string gfa_input = argv[1];
    std::string gaf_file = argv[2];
    int min_coverage = std::stoi(argv[3]);
    int threads = std::stoi(argv[4]);
    bool rename = std::stoi(argv[5]);
    std::string gfa_output = argv[6];
    bool contiguity = std::stoi(argv[7]);
    bool single_genome = std::stoi(argv[8]);
    int last_k = std::stoi(argv[9]);
    std::string logfile = argv[10];

    ofstream log(logfile);
    if (!log){
        std::cerr << "ERROR: graphunzip could not open log file " << logfile << endl;
        return 1;
    }

    // Export the command line used to the log file
    log << "Command line used: ";
    for (int i = 0; i < argc; ++i) {
        log << argv[i] << " ";
    }
    log << endl;

    //fail early if the inputs cannot be read or the output cannot be written
    if (!std::ifstream(gfa_input) || !std::ifstream(gaf_file)){
        std::cerr << "ERROR: graphunzip cannot read " << gfa_input << " or " << gaf_file << endl;
        log << "ERROR: cannot read " << gfa_input << " or " << gaf_file << endl;
        return 1;
    }
    {
        ofstream test_output(gfa_output);
        if (!test_output){
            std::cerr << "ERROR: graphunzip cannot write to " << gfa_output << endl;
            log << "ERROR: cannot write to " << gfa_output << endl;
            return 1;
        }
    }

    //load the segments from the GFA file
    unordered_map<string, int> segment_IDs;
    vector<Segment> segments; //segments is a dict of pairs
    load_GFA(gfa_input, segments, segment_IDs, false); //false is to not load the seqs in memory
    log << "Segments loaded" << endl;

    //load the paths from the GAF file
    //parse the GAF once into a memory-mapped binary file, then compute the consensuses by batches of segments
    ReadPathsInfo paths;
    {
        vector<uint64_t> occurrences(segments.size(), 0);
        if (!paths.store.build(gaf_file, gfa_output + ".paths.bin", segments, segment_IDs, occurrences)){
            log << "ERROR: could not build the binary representation of " << gaf_file << endl;
            return 1;
        }
        uint64_t memory_budget = 1ULL << 30; //bytes of occurrence index per batch of segments
        const char* budget_env = std::getenv("GRAPHUNZIP_BATCH_BYTES");
        if (budget_env != nullptr && std::strtoull(budget_env, nullptr, 10) > 0){
            memory_budget = std::strtoull(budget_env, nullptr, 10);
        }
        paths.compute(segments.size(), min_coverage, threads, occurrences, memory_budget);
    }

    log << "Paths loaded and haploid contigs determined" << endl;

    vector<Segment> unzipped_segments;
    unordered_map<int, std::vector<int>> old_ID_to_new_IDs;
    unordered_map<int, vector<set<int>>> already_built_bridges;
    already_built_bridges.reserve(segments.size());
    for (auto old_segment = 0 ; old_segment < segments.size() ; old_segment++){
       already_built_bridges[old_segment] = {set<int>(), set<int>()};
    }
    create_haploid_contigs(segments, unzipped_segments, old_ID_to_new_IDs, already_built_bridges, min_coverage, contiguity, paths);

    {
        UnrepresentedPathsAdder unrepresented_paths(segments, unzipped_segments, old_ID_to_new_IDs, min_coverage, paths);
        list_non_represented_paths(segments, already_built_bridges, min_coverage, paths, unrepresented_paths);
        unrepresented_paths.finish();
    }
    paths.release_paths(); //the read paths are not needed anymore

    vector<Segment> merged_segments;
    merge_adjacent_contigs(unzipped_segments, merged_segments, gfa_input, rename, threads);

    bool success = true;
    if (!contiguity){
        string gfa_output_tmp = gfa_output + "_tmp.gfa";
        output_graph(gfa_output_tmp, gfa_input, merged_segments);
        if (!graph_file_written(gfa_output_tmp, merged_segments)){
            std::cerr << "ERROR: graphunzip could not write " << gfa_output_tmp << endl;
            return 1;
        }

        //now merge the contigs
        segments.clear();
        segment_IDs.clear();
        load_GFA(gfa_output_tmp, segments, segment_IDs, false);
        vector<Segment> new_merged_segments;
        merge_adjacent_contigs(segments, new_merged_segments, gfa_output_tmp, rename, threads);
        output_graph(gfa_output, gfa_output_tmp, new_merged_segments);
        success = graph_file_written(gfa_output, new_merged_segments);

        //remove the temporary files
        remove(gfa_output_tmp.c_str());
    }
    else{
        string gfa_output_tmp = gfa_output + "_tmp.gfa";
        output_graph(gfa_output_tmp, gfa_input, merged_segments);
        if (!graph_file_written(gfa_output_tmp, merged_segments)){
            std::cerr << "ERROR: graphunzip could not write " << gfa_output_tmp << endl;
            return 1;
        }

        string gfa_output_tmp2 = gfa_output + "_tmp2.gfa";
        pop_and_shave_graph(gfa_output_tmp, min_coverage, 2*last_k+10, last_k, gfa_output_tmp2, 0, threads, single_genome);
        if (!std::ifstream(gfa_output_tmp2)){
            std::cerr << "ERROR: graphunzip could not write " << gfa_output_tmp2 << endl;
            remove(gfa_output_tmp.c_str());
            return 1;
        }
        //now merge the contigs
        segments.clear();
        segment_IDs.clear();
        load_GFA(gfa_output_tmp2, segments, segment_IDs, false);
        vector<Segment> merged_segments;
        merge_adjacent_contigs(segments, merged_segments, gfa_output_tmp2, rename, threads);
        output_graph(gfa_output, gfa_output_tmp2, merged_segments);
        success = graph_file_written(gfa_output, merged_segments);

        //remove the temporary files
        remove(gfa_output_tmp.c_str());
        remove(gfa_output_tmp2.c_str());
    }

    if (!success){
        std::cerr << "ERROR: graphunzip could not write the output graph " << gfa_output << endl;
        log << "ERROR: could not write the output graph " << gfa_output << endl;
        return 1;
    }
    log << "Output written to " << gfa_output << endl;
    log.close();
    if (!log){
        std::cerr << "ERROR: graphunzip could not write the log file " << logfile << endl;
        return 1;
    }
    return 0;
}
