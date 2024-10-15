#include "graphunzip.h"
#include "basic_graph_manipulation.h"

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

// string reverse_complement(string& seq){
//     string rc (seq.size(), 'N');
//     for (int i = seq.size() - 1 ; i >= 0; i--){
//         switch (seq[i]){
//             case 'A':
//                 rc[seq.size() - 1 - i] = 'T';
//                 break;
//             case 'T':
//                 rc[seq.size() - 1 - i] = 'A';
//                 break;
//             case 'C':
//                 rc[seq.size() - 1 - i] = 'G';
//                 break;
//             case 'G':
//                 rc[seq.size() - 1 - i] = 'C';
//                 break;
//             default:
//                 rc[seq.size() - 1 - i] = 'N';
//                 break;
//         }
//     }
//     return rc;
// }

void Segment::add_neighbor(vector<pair<int,bool>> new_neighbor, bool left){

    if (new_neighbor.size() == 0){
        return;
    }
    if (left){
        neighbors_left.push_back(new_neighbor);
    }
    else{
        neighbors_right.push_back(new_neighbor);
    }
}

/**
 * @brief Compute the consensual path of contigs left and right
 * 
 */
void Segment::compute_consensuses(){

    //right of the contig

    //first sort the neighbors by decreasing length. This way, a path seen after can never include a path seen before
    std::sort(neighbors_right.begin(), neighbors_right.end(), [](vector<pair<int,bool>> a, vector<pair<int,bool>> b){return a.size() > b.size();});

    vector<vector<pair<int,bool>>> all_paths_right; //inventoriate a subset of paths which contains all the other paths
    //iteratively go through the neighbors and check if they are subreads of inventoried paths. If not, add them to the list of all paths
    for (int n1 = 0 ; n1 < neighbors_right.size() ; n1++){
        int nb_of_supra_reads = 0;
        int idx_of_supra_read = -1;
        for (int n2 = 0 ; n2<all_paths_right.size() ; n2++){
            bool is_subread = true;
            for (int c = 0 ; c < neighbors_right[n1].size() ; c++){
                if (all_paths_right[n2][c] != neighbors_right[n1][c]){
                    is_subread = false;
                    break;
                }
            }
            if (is_subread){
                nb_of_supra_reads++;
                idx_of_supra_read = n2;
            }
        }
        if (nb_of_supra_reads == 0){
            all_paths_right.push_back(neighbors_right[n1]);
        }
    }

    //now inventoriate how many reads support each path
    for (auto new_path : all_paths_right){
        neighbors_right_with_strengths.push_back({});
        for (auto contig : new_path){
            neighbors_right_with_strengths[neighbors_right_with_strengths.size() - 1].push_back({contig, {0,0}});
        }
    }
    for (int n1 = 0 ; n1 < neighbors_right.size() ; n1++){
        for (int n2 = 0 ; n2<all_paths_right.size() ; n2++){
            for (int c = 0 ; c < neighbors_right[n1].size() ; c++){
                if (all_paths_right[n2][c] != neighbors_right[n1][c]){
                    neighbors_right_with_strengths[n2][c].second.second++;
                    break;
                }
                else{
                    neighbors_right_with_strengths[n2][c].second.first++;
                }
            }
        }
    }

    // if (this->name == "8997"){
    //     cout << "****here are the neighbors right with their strengths: " << endl;
    //     for (auto neighbor : neighbors_right_with_strengths){
    //         for (auto contig : neighbor){
    //             cout << contig.first.first << " ";
    //         }
    //         cout << endl;
    //         for (auto contig : neighbor){
    //             cout << contig.second.first << "/" << contig.second.second << " ";
    //         }
    //         cout << endl;
    //     }
    // }

    //go through the neighbors with strengths and see if there is a consensus
    bool haploid_right = false;
    for (auto neighbor : neighbors_right_with_strengths){
        bool consensus_here = true;
        for (pair<pair<int,bool>,pair<int,int>> contig : neighbor){
            if (0.2*(contig.second.first-1) > contig.second.second-1 && contig.second.second < 5){ //then it is relatively consensual
                if (contig.second.first > 1){
                    consensus_right.push_back(contig.first);
                }
            }
            else{
                consensus_right = {};
                consensus_here = false;
                break;
            }
        }
        if (consensus_here){
            haploid_right = true;
            break;
        }
    }

    //left of the contig

    //first sort the neighbors by decreasing length. This way, a path seen after can never include a path seen before
    std::sort(neighbors_left.begin(), neighbors_left.end(), [](vector<pair<int,bool>> a, vector<pair<int,bool>> b){return a.size() > b.size();});

    vector<vector<pair<int,bool>>> all_paths_left; //inventoriate a subset of paths which contains all the other paths
    //iteratively go through the neighbors and check if they are subreads of inventoried paths. If not, add them to the list of all paths
    for (int n1 = 0 ; n1 < neighbors_left.size() ; n1++){
        int nb_of_supra_reads = 0;
        int idx_of_supra_read = -1;
        for (int n2 = 0 ; n2<all_paths_left.size() ; n2++){
            bool is_subread = true;
            for (int c = 0 ; c < neighbors_left[n1].size() ; c++){
                if (all_paths_left[n2][c] != neighbors_left[n1][c]){
                    is_subread = false;
                    break;
                }
            }
            if (is_subread){
                nb_of_supra_reads++;
                idx_of_supra_read = n2;
            }
        }
        if (nb_of_supra_reads == 0){
            all_paths_left.push_back(neighbors_left[n1]);
        }
    }

    //now inventoriate how many reads support each path
    for (auto new_path : all_paths_left){
        neighbors_left_with_strengths.push_back({});
        for (auto contig : new_path){
            neighbors_left_with_strengths[neighbors_left_with_strengths.size() - 1].push_back({contig, {0,0}});
        }
    }

    for (int n1 = 0 ; n1 < neighbors_left.size() ; n1++){
        for (int n2 = 0 ; n2<all_paths_left.size() ; n2++){
            for (int c = 0 ; c < neighbors_left[n1].size() ; c++){
                if (all_paths_left[n2][c] != neighbors_left[n1][c]){
                    neighbors_left_with_strengths[n2][c].second.second++;
                    break;
                }
                else{
                    neighbors_left_with_strengths[n2][c].second.first++;
                }
            }
        }
    }

    //go through the neighbors with strengths and see if there is a consensus
    bool haploid_left = false;
    for (auto neighbor : neighbors_left_with_strengths){
        bool consensus_here = true;
        for (pair<pair<int,bool>,pair<int,int>> contig : neighbor){
            if (0.2*(contig.second.first-1) > contig.second.second-1 && contig.second.second < 5){ //then it is relatively consensual
                if (contig.second.first > 1){
                    consensus_left.push_back(contig.first);
                }
            }
            else{
                consensus_left = {};
                consensus_here = false;
                break;
            }
        }
        if (consensus_here){
            haploid_left = true;
            break;
        }
    }

    if (neighbors_left.size() == 0){
        haploid_left = true;
        consensus_left = {};
    }
    if (neighbors_right.size() == 0){
        haploid_right = true;
        consensus_right = {};
    }

    this->haploid = haploid_right && haploid_left;
}

vector<vector<pair<int,bool>>> Segment::get_strong_neighbors_left(int min_coverage){
    vector<vector<pair<int,bool>>> strong_neighbors;
    
    for (vector<pair<pair<int,bool>,pair<int,int>>> neighbor : neighbors_left_with_strengths){
        vector<pair<int,bool>> strong_neighbor;
        for (pair<pair<int,bool>,pair<int,int>> contig : neighbor){
            if (contig.second.first >= min_coverage && (contig.second.second == 0 || contig.second.first > 4)){
                strong_neighbor.push_back(contig.first);
            }
            else{
                break;
            }
        }
        if (strong_neighbor.size() > 0){
            strong_neighbors.push_back(strong_neighbor);
        }
    }
    
    return strong_neighbors;
}

vector<vector<pair<int,bool>>> Segment::get_strong_neighbors_right(int min_coverage){
    vector<vector<pair<int,bool>>> strong_neighbors;
    
    for (vector<pair<pair<int,bool>,pair<int,int>>> neighbor : neighbors_right_with_strengths){
        vector<pair<int,bool>> strong_neighbor;
        for (pair<pair<int,bool>,pair<int,int>> contig : neighbor){
            if (contig.second.first >= min_coverage && (contig.second.second == 0 || contig.second.first > 4)){
                strong_neighbor.push_back(contig.first);
            }
            else{
                break;
            }
        }
        if (strong_neighbor.size() > 0){
            strong_neighbors.push_back(strong_neighbor);
        }
    }
    
    return strong_neighbors;
}

/*
void load_GFA(string gfa_file, vector<Segment> &segments, unordered_map<string, int> &segment_IDs){
    //load the segments from the GFA file
    
    //in a first pass index all the segments by their name
    ifstream gfa(gfa_file);
    string line;
    while (getline(gfa, line)){
        if (line[0] == 'S'){
            long int pos_in_file = (long int) gfa.tellg() - line.size() - 1;
            stringstream ss(line);
            string nothing, name, seq;
            ss >> nothing >> name >> seq;

            double coverage = 0;
            //try to find a DP: tag
            string tag;
            while (ss >> tag){
                if (tag.substr(0, 3) == "DP:"){
                    coverage = std::stof(tag.substr(5, tag.size() - 5));
                }
            }

            Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, seq.size(), coverage);

            segment_IDs[name] = s.ID;
            segments.push_back(s);
        }
    }
    gfa.close();

    //in a second pass, load the links
    gfa.open(gfa_file);
    while (getline(gfa, line)){
        if (line[0] == 'L'){
            stringstream ss(line);
            string nothing, name1, name2;
            string orientation1, orientation2;
            string cigar;

            ss >> nothing >> name1 >> orientation1 >> name2 >> orientation2 >> cigar;

            int end1 = 1;
            int end2 = 0;

            if (orientation1 == "-"){
                end1 = 0;
            }
            if (orientation2 == "-"){
                end2 = 1;
            }

            int ID1 = segment_IDs[name1];
            int ID2 = segment_IDs[name2];

            //check that the link did not already exist
            bool already_exists = false;
            for (pair<int,int> link : segments[ID1].links[end1].first){
                if (link.first == ID2 && link.second == end2){
                    already_exists = true;
                }
            }
            if (!already_exists){

                segments[ID1].links[end1].first.push_back({ID2, end2});
                segments[ID1].links[end1].second.push_back(cigar);

                segments[ID2].links[end2].first.push_back({ID1, end1});
                segments[ID2].links[end2].second.push_back(cigar);
            }
        }
    }
    gfa.close();
}
*/

void load_GAF(string gaf_file, vector<Segment> &segments, unordered_map<string, int> &segments_IDs){
    //load the paths from the GAF file
    ifstream gaf(gaf_file);
    string line;
    int nb_lines = 0;
    while (getline(gaf, line)){

        //get the name of the read as the first space-delimited word, and the path as the 8th tab-delimited column
        string name_of_read;
        string path;
        string nothing;
        int field = 0;
        string field_content = "";
        stringstream line_stream(line);
        while (std::getline(line_stream, field_content, '\t')){
            if (field == 0){
                name_of_read = field_content.substr(0, field_content.find_first_of(' '));
            }
            else if (field == 5){
                path = field_content;
            }
            field++;
        }

        if (nb_lines % 100000 == 0){
            cout << nb_lines << " lines read in load_GAF\r" << std::flush;
        }
        // if (nb_lines > 1000000){
        //     cout << "NOT READING EVERYTHING FOR DEBGUGGING PUROSPS" << endl;
        //     break;
        // }
        nb_lines += 1;
        
        bool orientation_now;
        string str_now;
        vector<pair<int,bool>> segments_now;
        vector<vector<pair<int,bool>>> paths_in_this_read;

        for (int c = 0 ; c < path.size(); c++){
            if (path[c] == '>' || path[c] == '<' || c == path.size()-1){

                if (c == path.size()-1){
                    str_now += path[c];
                }

                if (c != 0){
                    if (segments_IDs.find(str_now) == segments_IDs.end()){
                        cout << "Error: the segment " << str_now << " is in the gaf but cannot be found in the GFA" << endl;
                        cout << line << endl;
                        exit(1);
                    }
                    int new_contig_ID = segments_IDs[str_now];

                    //check that the segment can indeed be found next to the previous segment
                    bool found = false;
                    if (segments_now.size() > 0){
                        for (pair<int,int> link : segments[segments_now.back().first].links[segments_now.back().second].first){
                            if (link.first == new_contig_ID && link.second == !orientation_now){
                                found = true;
                            }
                        }
                    }
                    if (segments_now.size() == 0 || found){
                        segments_now.push_back({new_contig_ID, orientation_now});
                    }
                    else{ //then the neighbor is not the neighbor of a previous segment, cut
                        paths_in_this_read.push_back(segments_now);
                        segments_now = {{new_contig_ID, orientation_now}};  
                    }
                }

                if (path[c] == '>'){
                    orientation_now = true;
                }
                else{
                    orientation_now = false;
                }
                str_now = "";
            }
            else{
                str_now += path[c];
            }
        }
        paths_in_this_read.push_back(segments_now);

        for (vector<pair<int,bool>> pair_path : paths_in_this_read){

            //compute the reverse path too
            vector<pair<int,bool>> reverse_path_contigs;
            for (int contig = pair_path.size() - 1 ; contig >= 0 ; contig--){
                reverse_path_contigs.push_back({pair_path[contig].first, !pair_path[contig].second});
            }


            for (int contig = 0 ; contig < pair_path.size()  ; contig++){
                if (pair_path[contig].second){
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(pair_path.begin() + contig + 1, pair_path.end()), false);
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(reverse_path_contigs.begin() + reverse_path_contigs.size() - contig, reverse_path_contigs.end()), true);
                } 
                else{
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(pair_path.begin() + contig + 1, pair_path.end()), true);
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(reverse_path_contigs.begin() + reverse_path_contigs.size() - contig, reverse_path_contigs.end()), false);
                }

            }
        }
    }
    gaf.close();
}

/**
 * @brief Determine which contigs are haploid.
 * 
 * @param segments 
 * @param min_coverage 
 */
void determine_haploid_contigs(vector<Segment> &segments, int min_coverage){
    //determine which contigs are haploid. For this, a simple metric: the paths should be consensual left and right
    for (auto s = 0 ; s < segments.size() ; s++){
    
        segments[s].compute_consensuses();
        if (segments[s].name == "4972" && false){

            cout << "****here are the neighbors right with their strengths: " << endl;
            auto neighbors_right_with_strengths = segments[s].get_neighbors_right_with_strengths();
            for (auto neighbor : neighbors_right_with_strengths){
                for (auto contig : neighbor){
                    cout << contig.first.first << " ";
                }
                cout << endl;
                for (auto contig : neighbor){
                    cout << contig.second.first << "/" << contig.second.second << " ";
                }
                cout << endl;
            }
            cout << "****here are the neighbors left with their strengths: " << endl;
            auto neighbors_left_with_strengths = segments[s].get_neighbors_left_with_strengths();
            for (auto neighbor : neighbors_left_with_strengths){
                for (auto contig : neighbor){
                    cout << contig.first.first << " ";
                }
                cout << endl;
                for (auto contig : neighbor){
                    cout << contig.second.first << "/" << contig.second.second << " ";
                }
                cout << endl;
            }
            
            cout << "haploid: " << segments[s].is_haploid() << endl;

            cout << "consensus left: ";
            for (auto contig : segments[s].get_consensus_left()){
                cout << contig.first << " ";
            }
            cout << endl;

            cout << "consensus right: ";
            for (auto contig : segments[s].get_consensus_right()){
                cout << contig.first << " ";
            }
            cout << endl;
        }
    }    
}

/**
 * @brief Unzip the graph and create a new list of segments. 
 * 
 * @param old_segments 
 * @param new_segments 
 * @param min_coverage 
 * @param unordered_map<int, vector<set<int>>> already_built_bridges; //associates to each haploid contig two sets, one of the contigs it is connected to on the left, one on the right

 */
void create_haploid_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, unordered_map<int, vector<set<int>>> &already_built_bridges, int min_coverage, bool contiguity){

    new_segments = {};

    //associates all old contigs ids to their new IDs (can be multiple if contig is duplicated)
    int new_ID = 0;


    //begin by building bridges between haploid contigs
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){

        // if (old_segments[old_segment].name != "2" && old_segments[old_segment].name != "666603"){
        //     continue;
        // }
        
        if (old_segments[old_segment].is_haploid()){

            // cout << "in haploid contig " << old_segments[old_segment].name << endl;

            //first create the contig if needed
            if (old_ID_to_new_IDs.find(old_segment) == old_ID_to_new_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_segment].name + "_0" , new_ID, old_segments[old_segment].get_pos_in_file(), old_segments[old_segment].get_length(), old_segments[old_segment].get_coverage()));
                old_ID_to_new_IDs[old_segment] = {new_ID};
                new_ID++;
            }
            int ID_of_new_contig = old_ID_to_new_IDs[old_segment][0];

            //now see if it already has neighbors left and right and create them until the next haploid contig if not

            //right
            vector<pair<int,bool>> cons_right = old_segments[old_segment].get_consensus_right();
            bool there_is_a_bridge = false;
            double coverage_bridge = old_segments[old_segment].get_coverage();
            int length_bridge = old_segments[old_segment].get_length();
            int otherEndOfBridge = 0;
            int endOfOtherEndOfBridge = 0;
            for (auto contig_and_orientation : cons_right){
                length_bridge += old_segments[contig_and_orientation.first].get_length();
                if (old_segments[contig_and_orientation.first].is_haploid()){
                    //in case contiguity mode is on, check that the bridge is reciprocal
                    bool reciprocal = false;
                    if (contig_and_orientation.second == 1 && contiguity){
                        // if (old_segments[old_segment].name == "4925"){
                        //     cout << "checked againts: ";
                        //     for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_left()){
                        //         cout << neighbor.first << " ";
                        //     }
                        // }
                        for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_left()){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }
                    if (contig_and_orientation.second == 0 && contiguity){
                        // if (old_segments[old_segment].name == "4925"){
                        //     cout << "checked2 againts: ";
                        //     for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_right()){
                        //         cout << neighbor.first << " ";
                        //     }
                        // }
                        for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_right()){
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
            // if (old_segments[old_segment].name == "4925"){
            //     cout << "there is a bridge right: " << there_is_a_bridge << " with " << otherEndOfBridge << " and " << endOfOtherEndOfBridge << endl;
            // }
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
                    if (!old_segments[cons_right[contig_in_cons].first].is_haploid()){
                        double coverage_of_this_contig = std::min(coverage_bridge, old_segments[cons_right[contig_in_cons].first].get_coverage());
                        new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_right[contig_in_cons].first].size()) , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), length_bridge, coverage_of_this_contig));
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
                            new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), length_bridge, coverage_bridge));
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
                    if (old_segments[cons_right[contig_in_cons].first].is_haploid()){
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
            vector<pair<int,bool>> cons_left = old_segments[old_segment].get_consensus_left();
            there_is_a_bridge = false;
            coverage_bridge = old_segments[old_segment].get_coverage();
            length_bridge = old_segments[old_segment].get_length();
            for (auto contig_and_orientation : cons_left){
                length_bridge += old_segments[contig_and_orientation.first].get_length();
                if (old_segments[contig_and_orientation.first].is_haploid()){
                    bool reciprocal = false;

                    if (contig_and_orientation.second == 1 && contiguity){
                        for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_left()){
                            if (neighbor.first == old_segment){
                                reciprocal = true;
                            }
                        }
                    }
                    if (contig_and_orientation.second == 0 && contiguity){
                        for (auto neighbor : old_segments[contig_and_orientation.first].get_consensus_right()){
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
                    if (!old_segments[cons_left[contig_in_cons].first].is_haploid()){
                        double coverage_of_this_contig = std::min(coverage_bridge, old_segments[cons_left[contig_in_cons].first].get_coverage());
                        new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_left[contig_in_cons].first].size()) , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), length_bridge, coverage_of_this_contig));
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
                            new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), length_bridge, coverage_bridge));
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
                    if (old_segments[cons_left[contig_in_cons].first].is_haploid()){
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
 * @brief Check what is not seen yet in new_segments
 * 
 * @param old_segments 
 * @param new_segments 
 * @param old_ID_to_new_IDs 
 * @param min_coverage 
 * @param all_paths
 * @return vector<vector<pair<int,bool>>> all the non represented paths
 */
vector<vector<pair<int,bool>>> list_non_represented_paths(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, unordered_map<int, vector<set<int>>> &already_built_bridges, int min_coverage){

    //first index represented paths
    vector<vector<pair<int,bool>>> represented_paths;
    unordered_map<int, vector<pair<int,int>>> where_is_this_contig_represented; //associates to an ID (indices of the path, position in the path)

    //go through all the haploids segments and list the represented paths left and right in new_segments
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){
        if (old_segments[old_segment].is_haploid()){
            
            std::vector<std::pair<int,bool>> consensus_left = {{old_segments[old_segment].ID, false}};
            std::vector<std::pair<int,bool>> new_elements = old_segments[old_segment].get_consensus_left();
            consensus_left.insert(consensus_left.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            bool there_is_a_bridge = false;
            int idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_left){
                if (old_segments[contig_and_orientation.first].is_haploid() && idx_of_last_haploid_contig != 0){
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

                    if (old_segments[contig_and_orientation.first].is_haploid()){
                        number_of_haploid_contigs++;
                        if (number_of_haploid_contigs > 1){
                            break;
                        }
                    }
                    idx++;
                }
            }

            vector<pair<int,bool>> consensus_right = {{old_segments[old_segment].ID, true}};
            new_elements = old_segments[old_segment].get_consensus_right();
            consensus_right.insert(consensus_right.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            there_is_a_bridge = false;
            idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_right){
                if (old_segments[contig_and_orientation.first].is_haploid() && idx_of_last_haploid_contig != 0){
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

                    if (old_segments[contig_and_orientation.first].is_haploid()){
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

    // cout << "finished indexing represented paths" << endl;
    // cout << "for example, here are all the paths where I find 1: " << endl;
    // for (pair<int,int> path_and_pos : where_is_this_contig_represented[1]){
    //     cout << "path " << path_and_pos.first << " at position " << path_and_pos.second << " : ";
    //     for (pair<int,bool> contig : represented_paths[path_and_pos.first]){
    //         cout << contig.first << " ";
    //     }
    //     cout << endl;
    // }

    //now go through all the paths and see if they are represented
    vector<vector<pair<int,bool>>> unrepresented_paths;
    int idx = 0;
    for (Segment s : old_segments){
        vector<vector<pair<int,bool>>> paths_segment = s.get_strong_neighbors_left(min_coverage);

        int number_of_left_paths = paths_segment.size();
        vector<vector<pair<int,bool>>> neighbors_right = s.get_strong_neighbors_right(min_coverage);

        paths_segment.insert(paths_segment.end(), neighbors_right.begin(), neighbors_right.end());

        int idx_here = 0;
        for (vector<pair<int,bool>> path : paths_segment){

            if (idx % 100000 == 0){
                cout << idx << " paths checked\r" << std::flush;
            }

            //check if the path is represented
            vector<pair<int,bool>> path_until_haploid_contig;
            if (idx_here < number_of_left_paths){
                path_until_haploid_contig = {{s.ID, false}};
            }
            else{
                path_until_haploid_contig = {{s.ID, true}};
            }
            for (pair<int,bool> contig : path){

                // if (s.name == "446"){
                //     cout << "checking if paeeth " << idx_here << " ";
                //     for (pair<int,bool> contig : path_until_haploid_contig){
                //         cout << old_segments[contig.first].name << " ";
                //     }
                //     cout << " is represented" << endl;
                // }

                if (!old_segments[contig.first].is_haploid()){
                    path_until_haploid_contig.push_back(contig);
                }
                else{
                    path_until_haploid_contig.push_back(contig);
                    bool found = false;
                    if (where_is_this_contig_represented.find(contig.first) != where_is_this_contig_represented.end()){
                        for (pair<int,int> path_and_pos : where_is_this_contig_represented[contig.first]){
                            vector<pair<int,bool>> path_to_check;
                            if (path_and_pos.second == 0){
                                path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].begin(), represented_paths[path_and_pos.first].begin()+std::min(path_until_haploid_contig.size(),represented_paths[path_and_pos.first].size()));
                            }
                            else {
                                path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].end()-std::min(path_until_haploid_contig.size(),represented_paths[path_and_pos.first].size()), represented_paths[path_and_pos.first].end());
                            }
                            // if (s.name == "446"){
                            //     cout << "represented path: ";
                            //     for (pair<int,bool> contig : represented_paths[path_and_pos.first]){
                            //         cout << old_segments[contig.first].name << " ";
                            //     }
                            //     cout << endl;
                            //     cout << "path to check: ";
                            //     for (pair<int,bool> contig : path_to_check){
                            //         cout << old_segments[contig.first].name << " ";
                            //     }
                            //     cout << endl;
                            // }
                            if (path_until_haploid_contig == path_to_check){
                                found = true;
                            }
                        }
                    }
                    if (!found){
                        unrepresented_paths.push_back(path_until_haploid_contig);
                        // for (pair<int,bool> contig : path_until_haploid_contig){
                            // if (old_segments[contig.first].name == "1530"){
                            //     cout << "Here is a nzzon represented path contianing 9661: ";
                            //     for (pair<int,bool> contig : path_until_haploid_contig){
                            //         cout << old_segments[contig.first].name << " ";
                            //     }
                            //     cout << " (computed from segment " << s.name << endl;
                            // }
                        // }
                    }
                    path_until_haploid_contig = {contig};
                }
            }
            //check the last path
            if (path_until_haploid_contig.size() > 0){
                bool found = false;

                for (pair<int,int> path_and_pos : where_is_this_contig_represented[path_until_haploid_contig[0].first]){
                    vector<pair<int,bool>> path_to_check;

                    // if (s.name == "446"){
                    //     cout << "represented path: ";
                    //     for (pair<int,bool> contig : represented_paths[path_and_pos.first]){
                    //         cout << old_segments[contig.first].name << " ";
                    //     }
                    //     cout << endl;
                    // }

                    if (represented_paths[path_and_pos.first].size() >= path_and_pos.second + path_until_haploid_contig.size()){
                        path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].begin()+path_and_pos.second, represented_paths[path_and_pos.first].begin() + path_and_pos.second + path_until_haploid_contig.size());
                    
                        // if (s.name == "13278"){
                        //     cout << "path to check: ";
                        //     for (pair<int,bool> contig : path_to_check){
                        //         cout << old_segments[contig.first].name << " ";
                        //     }
                        //     cout << endl;
                        // }
                        
                        if (path_until_haploid_contig == path_to_check){
                            found = true;
                        }
                    }
                }
                if (!found){
                    unrepresented_paths.push_back(path_until_haploid_contig);
                    //check if 2392 is in the path
                    // for (pair<int,bool> contig : path_until_haploid_contig){
                    //     if (old_segments[contig.first].name == "1530"){
                    //         cout << "Here is a nzzon represented path contianing 9661: ";
                    //         for (pair<int,bool> contig : path_until_haploid_contig){
                    //             cout << old_segments[contig.first].name << " ";
                    //         }
                    //         cout << " (computed from segment " << s.name << endl;
                    //     }
                    // }
                }
            }
            idx_here++;
            idx++;
        }
    }

    return unrepresented_paths;

}

/**
 * @brief Take a list of non represented path in the graph and add contigs and link until all paths are represented
 * 
 * @param old_segments 
 * @param new_segments 
 * @param old_ID_to_new_IDs 
 * @param min_coverage 
 * @param unrepresented_paths 
 */
void add_unrepresented_paths(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, int min_coverage, vector<vector<pair<int,bool>>>& unrepresented_paths){

    unordered_map<int, int> old_IDs_to_new_non_haploid_IDs; //associates old IDs to new IDs for the contig we are going to create
    for (Segment s: old_segments){ //we are not going to create new versions of haploid contigs
        if (s.is_haploid()){
            old_IDs_to_new_non_haploid_IDs[s.ID] = old_ID_to_new_IDs[s.ID][0];
        }
    }

    //convert unrepresented path in a list of links that must be there in the final graph
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add;
    int index = 0;

    for (vector<pair<int,bool>> path : unrepresented_paths){

        // cout << "adding the contigs of path " << endl;
        // for (pair<int,bool> contig : path){
        //     cout << contig.first << " ";
        // }
        // cout << endl;

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
            for (pair<int,int> link : old_segments[old_ID2].links[end2].first){
                if (link.first == old_ID1 && link.second == end1){
                    cigar = old_segments[old_ID2].links[end2].second[idx];
                }
                idx ++;
            }

            // if ("10336_0" == new_segments[new_ID1].name){
            //     cout << "adding link between " << new_segments[new_ID1].name << " and " << new_segments[new_ID2].name << " with cigar " << cigar << endl;
            // }
            // else if ("10336_0" == new_segments[new_ID2].name){
            //     cout << "adding2 link between " << new_segments[new_ID1].name << " and " << new_segments[new_ID2].name << " with cigar " << cigar << endl;
            // }

            if (links_to_add.find({{pair<int,int>(new_ID2, end2), pair<int,int>(new_ID1, end1)}, cigar}) == links_to_add.end()){
                links_to_add.insert({{pair<int,int>(new_ID1, end1), pair<int,int>(new_ID2, end2)}, cigar});
            }
        }
        // if (index > 3){
        //     cout << "breaking" << endl;
        //     break;
        // }
        index++;
    }

    //add the links
    for (pair<pair<pair<int,int>, pair<int,int>>, string> link : links_to_add){
        new_segments[link.first.first.first].links[link.first.first.second].first.push_back({link.first.second.first, link.first.second.second});
        new_segments[link.first.first.first].links[link.first.first.second].second.push_back(link.second);
        new_segments[link.first.second.first].links[link.first.second.second].first.push_back({link.first.first.first, link.first.first.second});
        new_segments[link.first.second.first].links[link.first.second.second].second.push_back(link.second);
    }
}

/**
 * @brief Takes in a graph and duplicates contigs that can be duplicated to improve completeness
 * 
 * @param segments The graph
 * @param min_relative_coverage The minimum relative coverage to consider a neighbor as sufficiently solid to entail duplication
 */
void duplicate_contigs(vector<Segment> &segments, float min_relative_coverage, float min_absolute_coverage, float relative_lengths_difference){

    //the basic idea of this function is that if a contig is the unique neighbor of several contigs, it can be duplicated to be next to all of them

    int number_of_changes = 1;

    while (number_of_changes > 0){
        number_of_changes = 0;
        int index_segment = 0;
        int original_segments_size = segments.size();
        while (index_segment < original_segments_size){

            Segment s = segments[index_segment];
            for (int end : {0,1}){

                if (s.links[end].first.size() <= 1){
                    continue;
                }

                //std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>  links;
                bool duplicable = true;
                float highest_coverage = 0;
                float lowest_coverage = 10000;
                float smallest_length = 1000000;
                // cout << "checking2 segment " << s.name << " at end " << end << " " << s.links.size() << endl;
                for (int n = 0 ; n < s.links[end].first.size() ; n++){
                    if (segments[s.links[end].first[n].first].links[s.links[end].first[n].second].first.size() != 1){
                        duplicable = false;
                        break;
                    }
                    if (segments[s.links[end].first[n].first].get_coverage() < min_absolute_coverage){
                        duplicable = false;
                        break;
                    }
                    if (s.links[end].first[n].first == s.ID){ //self loop, do NOT duplicate
                        duplicable = false;
                        break;
                    }
                    if (segments[s.links[end].first[n].first].get_coverage() > highest_coverage){
                        highest_coverage = segments[s.links[end].first[n].first].get_coverage();
                    }
                    if (segments[s.links[end].first[n].first].get_coverage() < lowest_coverage || lowest_coverage == 0){
                        lowest_coverage = segments[s.links[end].first[n].first].get_coverage();
                    }
                    if (segments[s.links[end].first[n].first].get_length() < smallest_length){
                        smallest_length = segments[s.links[end].first[n].first].get_length();
                    }
                }
                //don't duplicate because of a very small bubble
                if (highest_coverage*min_relative_coverage > lowest_coverage || smallest_length < relative_lengths_difference*s.get_length()){ 
                    duplicable = false;
                }
                //check that the coverage are coherent with a duplication
                double total_coverage = 0;
                double max_coverage = 0;
                if (duplicable){
                    for (int n = 0 ; n < s.links[end].first.size() ; n++){
                        total_coverage += segments[s.links[end].first[n].first].get_coverage();
                        if (segments[s.links[end].first[n].first].get_coverage() > max_coverage){
                            max_coverage = segments[s.links[end].first[n].first].get_coverage();
                        }
                    }
                    if (s.get_coverage() <= 0.8*total_coverage || s.get_coverage() <= max_coverage){
                        duplicable = false;
                    }
                }
                

                if (duplicable){
                    // number_of_changes++;
                    // cout << "I should duplicate " << s.name << " at end " << end << endl;

                    //create all the new contigs
                    for (int n = 0 ; n < s.links[end].first.size() ; n++){
                        string new_name = s.name + "_dup" + std::to_string(n);
                        double coverage_of_this_neighbor = s.get_coverage() * segments[s.links[end].first[n].first].get_coverage() / total_coverage;
                        segments.push_back(Segment(new_name, segments.size(), s.get_pos_in_file(), s.get_length(), coverage_of_this_neighbor));
                        segments[segments.size()-1].seq = s.seq;
                        //now add the link to the contig at the end
                        segments[segments.size()-1].links[end].first.push_back({s.links[end].first[n].first, s.links[end].first[n].second});
                        segments[segments.size()-1].links[end].second.push_back(s.links[end].second[n]);
                        segments[s.links[end].first[n].first].links[s.links[end].first[n].second].first.push_back({segments.size()-1, end});
                        segments[s.links[end].first[n].first].links[s.links[end].first[n].second].second.push_back(s.links[end].second[n]);
                        //add all the links on the other side
                        for (int m = 0 ; m < s.links[1-end].first.size() ; m++){
                            segments[segments.size()-1].links[1-end].first.push_back({s.links[1-end].first[m].first, s.links[1-end].first[m].second});
                            segments[segments.size()-1].links[1-end].second.push_back(s.links[1-end].second[m]);
                            segments[s.links[1-end].first[m].first].links[s.links[1-end].first[m].second].first.push_back({segments.size()-1, 1-end});
                            segments[s.links[1-end].first[m].first].links[s.links[1-end].first[m].second].second.push_back(s.links[1-end].second[m]);
                        }
                    }

                    //delete all the links of the contig
                    for (int delete_end : {0,1}){
                        for (int n = 0 ; n < s.links[delete_end].first.size() ; n++){
                            int other_ID = s.links[delete_end].first[n].first;
                            int other_end = s.links[delete_end].first[n].second;
                            int idx = 0;
                        
                            for (pair<int,int> link : segments[other_ID].links[other_end].first){
                                if (link.first == s.ID && link.second == delete_end){
                                    segments[other_ID].links[other_end].first.erase(segments[other_ID].links[other_end].first.begin() + idx);
                                    segments[other_ID].links[other_end].second.erase(segments[other_ID].links[other_end].second.begin() + idx);
                                    // cout << "deleting link between " << s.name << " and " << segments[other_ID].name << endl;
                                    break;
                                }
                                idx++;
                            }
                        }
                    }
                    segments[index_segment].links = {{{},{}},{{},{}}};
                    // cout << "deleting semgent " << s.name << endl;
                    segments[index_segment].name = "delete_me";
                } 
            }
            index_segment++;
        }
    }

}

/**
 * @brief Merge all old_segments into new_segments
 * 
 * @param old_segments 
 * @param new_segments 
 * @param original_gfa_file original gfa file to retrieve the sequences
 * @param rename rename the contigs in short names or keep the old names with underscores in between
 * @return * void 
 */
/*
void merge_adjacent_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, string original_gfa_file, bool rename){

    set<int> already_looked_at_segments; //old IDs of segments that have already been looked at and merged (don't want to merge them twice)
    int seg_idx = 0;
    unordered_map<pair<int,int>,pair<int,int>> old_ID_to_new_ID; //associates (old_id, old end) with (new_id_new_end)
    int number_of_merged_contigs = 0;
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add; //list of links to add, all in old IDs
    for (Segment old_seg : old_segments){

        // if (old_seg.name != "1665420_0"){
        //     seg_idx++;
        //     continue;
        // }

        if (already_looked_at_segments.find(old_seg.ID) != already_looked_at_segments.end()){
            seg_idx++;
            continue;
        }
        //check if it has either at least two neighbors left or that its neighbor left has at least two neighbors right
        // cout << "in merge, looking at segment " << seg_idx << " out of " << old_segments.size() << "\r" << std::flush;
        bool dead_end_left = false;
        if (old_seg.links[0].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].links[old_seg.links[0].first[0].second].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].ID == old_seg.ID){
            dead_end_left = true;
        }

        bool dead_end_right = false;
        if (old_seg.links[1].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].links[old_seg.links[1].first[0].second].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].ID == old_seg.ID){
            dead_end_right = true;
        }

        if (!dead_end_left && !dead_end_right){ //means this contig is in the middle of a long haploid contig, no need to merge
            seg_idx++;
            continue;
        }


        //create a new contig
        if (dead_end_left && dead_end_right){
            already_looked_at_segments.insert(old_seg.ID);
            string name = old_seg.name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), old_seg.get_length(), old_seg.get_coverage()));
            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 1};
            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[0].first){
                links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_seg.links[1].first){
                links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                idx_link++;
            }
            // cout << "yupee" << endl;
        }
        else if (dead_end_left){
            
            //let's see how far we can go right
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            vector<int> all_lengths = {old_seg.get_length()};
            int current_ID = old_seg.ID;
            int current_end = 1;

            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 1
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
                all_lengths.push_back(old_segments[current_ID].get_length());
            }
            already_looked_at_segments.insert(current_ID);

            //create the new contig
            string new_name = "";
            for (string name : all_names){
                new_name += name + "_";
            }
            new_name = new_name.substr(0, new_name.size()-1);
            string new_seq = "";
            for (string seq : all_seqs){
                new_seq += seq;
            }
            double new_coverage = 0;
            int new_length = 0;
            int idx = 0;
            for (double coverage : all_coverages){
                new_coverage += all_coverages[idx]*all_lengths[idx];
                new_length += all_lengths[idx];
                idx++;
            }
            new_coverage = new_coverage/new_length;
            string name = new_name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[0].first){
                links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "hurray" << endl;
        }
        else{
            
            //let's see how far we can go left
            vector<string> all_names = {old_seg.name};
            string seq = old_seg.get_seq(original_gfa_file);
            vector<string> all_seqs = {reverse_complement(seq)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            vector<int> all_lengths = {old_seg.get_length()};
            int current_ID = old_seg.ID;
            int current_end = 0;

            // cout << "exploring all the contigs left" << endl;
            // cout << "first exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name << endl;
            
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 0
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
                all_lengths.push_back(old_segments[current_ID].get_length());
                // cout << "exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name  << " " << old_segments[current_ID].links[current_end].first.size() << endl;
                // cout << "here are all the links : ";
                // for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                //     cout << link.first << " " << old_segments[link.first].name << " " << link.second << " ; ";
                // }
                // cout << endl;
            }
            already_looked_at_segments.insert(current_ID);

            //create the new contig
            string new_name = "";
            for (string name : all_names){
                new_name += name + "_";
            }
            new_name = new_name.substr(0, new_name.size()-1);
            string new_seq = "";
            for (string seq : all_seqs){
                new_seq += seq;
            }
            double new_coverage = 0;
            int new_length = 0;
            int idx = 0;
            for (double coverage : all_coverages){
                new_coverage += all_coverages[idx]*all_lengths[idx];
                new_length += all_lengths[idx];
                idx++;
            }
            new_coverage = new_coverage/new_length;
            string name = new_name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[1].first){
                links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "yay" << endl;
        }

        seg_idx++;
    }

    //some contigs are left: the ones that were in circular rings... go through them and add them
    for (Segment old_seg : old_segments){
        if (already_looked_at_segments.find(old_seg.ID) == already_looked_at_segments.end()){
            int current_ID = old_seg.ID;
            int current_end = 1;
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            vector<int> all_lengths = {old_seg.get_length()};
            bool circular_as_expected = true;
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                if (old_segments[current_ID].links[current_end].first[0].first == old_seg.ID){
                    break;
                }
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 0
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
                all_lengths.push_back(old_segments[current_ID].get_length());
            }
            if (old_segments[current_ID].links[current_end].first[0].first != old_seg.ID){
                circular_as_expected = false;
            }
            if (circular_as_expected){
                already_looked_at_segments.insert(current_ID);
                string new_name = "";
                for (string name : all_names){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (string seq : all_seqs){
                    new_seq += seq;
                }
                double new_coverage = 0;
                int new_length = 0;
                int idx = 0;
                for (double coverage : all_coverages){
                    new_coverage += all_coverages[idx]*all_lengths[idx];
                    new_length += all_lengths[idx];
                    idx++;
                }
                new_coverage = new_coverage/new_length;
                string name = new_name;
                if (rename){
                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }
                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
                new_segments[new_segments.size()-1].seq = new_seq;

                //add the links
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[0].first){
                    links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                    idx_link++;
                }

                old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};
            }
            else{
                cout << "ERROR: contig " << old_seg.name << " was discarded while merging the reads in graphunzip.cpp" << endl;
            }
        }
    }

    //now add the links in the new segments
    for (pair<pair<pair<int,int>, pair<int,int>>, string> link : links_to_add){
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].first.push_back(old_ID_to_new_ID[link.first.second]);
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].second.push_back(link.second);
    }
}
*/

/*
void output_graph(string gfa_output, string gfa_input, vector<Segment> &segments){
    ofstream gfa(gfa_output);
    for (Segment s : segments){
        gfa << "S\t" << s.name << "\t" << s.get_seq(gfa_input) << "\tDP:f:" << s.get_coverage() <<  "\n";
    }
    for (Segment s : segments){
        for (int end = 0 ; end < 2 ; end++){
            for (int neigh = 0 ; neigh < s.links[end].first.size() ; neigh++){

                //to make sure the link is not outputted twice
                if (s.ID > s.links[end].first[neigh].first || (s.ID == s.links[end].first[neigh].first && end > s.links[end].first[neigh].second) ){
                    continue;
                }

                string orientation = "+";
                if (end == 0){
                    orientation = "-";
                }
                gfa << "L\t" << s.name << "\t" << orientation << "\t" << segments[s.links[end].first[neigh].first].name << "\t";
                if (s.links[end].first[neigh].second == 0){
                    gfa << "+\t";
                }
                else{
                    gfa << "-\t";
                }
                gfa << s.links[end].second[neigh] << "\n";
            }
        }
    }
    gfa.close();
}
*/

int main(int argc, char *argv[])
{


    //HS_GraphUnzip <gfa_input> <gaf_file> <threads> <gfa_output> <exhaustive>
    if (argc != 9){
        //if -h or --help is passed as an argument, print the help
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)){
            std::cout << "Usage: graphunzip <gfa_input> <gaf_file> <min_coverage> <threads> <rename> <gfa_output> <contigutiy> <logfile>" << std::endl;
            return 0;
        }
        std::cout << "Usage: graphunzip <gfa_input> <gaf_file> <min_coverage> <threads> <rename> <gfa_output> <contiguity> <logfile>" << std::endl;
        return 1;
    }

    std::string gfa_input = argv[1];
    std::string gaf_file = argv[2];
    int min_coverage = std::stoi(argv[3]);
    int threads = std::stoi(argv[4]);
    bool rename = std::stoi(argv[5]);
    std::string gfa_output = argv[6];
    bool contiguity = std::stoi(argv[7]);
    std::string logfile = argv[8];

    ofstream log(logfile);


    //load the segments from the GFA file
    unordered_map<string, int> segment_IDs;
    vector<Segment> segments; //segments is a dict of pairs
    load_GFA(gfa_input, segments, segment_IDs);
    // cout << "Segments loaded" << endl;
    log << "Segments loaded" << endl;

    //load the paths from the GAF file
    load_GAF(gaf_file, segments, segment_IDs);
    // cout << "Paths loaded" << endl;
    log << "Paths loaded" << endl;

    //unzip the graph
    // cout << "Now unzipping the graph" << endl;
    vector<Segment> unzipped_segments;
    unordered_map<int, std::vector<int>> old_ID_to_new_IDs;
    unordered_map<int, vector<set<int>>> already_built_bridges;
    for (auto old_segment = 0 ; old_segment < segments.size() ; old_segment++){
       already_built_bridges[old_segment] = {set<int>(), set<int>()};
    }
    determine_haploid_contigs(segments, min_coverage);
    // cout << "Haploid contigs determined" << endl;
    // cout << "Creating haploid contigs" << endl;
    create_haploid_contigs(segments, unzipped_segments, old_ID_to_new_IDs, already_built_bridges, min_coverage, contiguity);
    // cout << "Haploid contigs created" << endl;

    // //output the graph
    // output_graph("out_alice/tmp/haploid.gfa", gfa_input, unzipped_segments);
    // exit(0);

    // cout << "Listing non represented paths" << endl;
    vector<vector<pair<int,bool>>> non_represented_paths = list_non_represented_paths(segments, unzipped_segments, old_ID_to_new_IDs, already_built_bridges, min_coverage);
    // cout << "Non represented paths listed" << endl;

    // cout << "Adding non represented paths" << endl;
    add_unrepresented_paths(segments, unzipped_segments, old_ID_to_new_IDs, min_coverage, non_represented_paths);
    // cout << "Non represented paths added" << endl;

    //output the graph
    // output_graph("out_alice/tmp/before_merge.gfa", gfa_input, unzipped_segments);
    // exit(0);


    // cout << "Graph unzipped" << endl;

    // cout << "Merging adjacent contigs" << endl;
    vector<Segment> merged_segments;
    merge_adjacent_contigs(unzipped_segments, merged_segments, gfa_input, rename);
    // cout << "Adjacent contigs merged" << endl;

    // output_graph("out_alice/tmp/before_dup.gfa", gfa_input, merged_segments);

    // cout << "Duplicating contigs" << endl;
    duplicate_contigs(merged_segments, 0.1, 5, 0.2);
    // cout << "Contigs duplicated" << endl;

    // output_graph("out_alice/tmp/after_dup.gfa", gfa_input, merged_segments);

    if (!contiguity){
        // cout << "Outputting the graph" << endl;
        output_graph(gfa_output, gfa_input, merged_segments);
        // cout << "Graph outputted" << endl;
    }
    else{
        string gfa_output_tmp = gfa_output + "_tmp.gfa";
        output_graph(gfa_output_tmp, gfa_input, merged_segments);

        string gfa_output_tmp2 = gfa_output + "_tmp2.gfa";
        pop_and_shave_graph(gfa_output_tmp, 5, 100, true, 31, gfa_output_tmp2, 0);
        //now merge the contigs
        segments.clear();
        segment_IDs.clear();
        load_GFA(gfa_output_tmp2, segments, segment_IDs);
        vector<Segment> merged_segments;
        merge_adjacent_contigs(segments, merged_segments, gfa_output_tmp2, rename);
        output_graph(gfa_output, gfa_output_tmp2, merged_segments);

        //remove the temporary files
        remove(gfa_output_tmp.c_str());
        remove(gfa_output_tmp2.c_str());
    }

    
}