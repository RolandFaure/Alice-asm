import re
import sys
import segment as sg
from input_output import read_GAF
from copy import deepcopy

import numpy as np

class Path :

    def __init__(self, contigs, orientations, read_name) :
        if len(contigs) != len(orientations) :
            raise ValueError("ERROR in simple_unzip.py: iiox")
        self.__contigs = contigs
        self.__orientations = []
        self.__read_name = read_name
        for i in orientations :
            if i == ">":
                self.__orientations.append(1)
            else :
                self.__orientations.append(0)

    def __len__(self):
        return len(self.__contigs)
    
    def __str__(self):
        string = ""
        for c in range(len(self.__contigs)) :
            string += "<>"[self.__orientations[c]] + str(self.__contigs[c].names)+ ", "
        return string

    def name(self):
        return self.__read_name
    
    def get_contigs(self):
        return self.__contigs
    
    def get_orientations(self):
        return self.__orientations
    
    def replace(self, contig_before, contig_after, pos_of_contig) :
        if self.__contigs[pos_of_contig] != contig_before :
            print("ERROR in simple_unzip.py: 9u8u")
            sys.exit()
        else:     
            self.__contigs[pos_of_contig] = contig_after
    
    def cancel(self, contig): #empty the path if it contains this contig
        co = 0
        for c in self.__contigs :
            if c.ID == contig.ID :
                # print("Cancelling ffry ", self.__read_name, " ", contig.names)
                # if co > 0 and co < len(self.__contigs)-1 :
                #     print("ERROR: cancelddling a path in the middle of it: ", self)
                #     sys.exit()
                self.__contigs = []
                self.__orientations = []
                return
            co += 1

    #check if the contigs are actually linked
    def split_if_invalid(self) :
        all_coherent_subpaths = []
        last_index = 0
        c = 0
        while c < len(self.__contigs)-1:
            index = sg.find_this_link(self.__contigs[c+1], 1-self.__orientations[c+1] , self.__contigs[c].links[self.__orientations[c]], self.__contigs[c].otherEndOfLinks[self.__orientations[c]])
            if index == -1 :
                # print("No link between ", self.__contigs[c+1].names , " and ", self.__contigs[c].names)
                # self.__contigs = []
                # self.__orientations = []
                all_coherent_subpaths += [Path(self.__contigs[last_index:c+1], ["<>"[i] for i in self.__orientations[last_index:c+1]], self.__read_name)]
                last_index = c+1
            c+= 1

        all_coherent_subpaths += [Path(self.__contigs[last_index:c+1], ["<>"[i] for i in self.__orientations[last_index:c+1]], self.__read_name)]
        return all_coherent_subpaths

    #get rid of the ends of the path that are just straight lines
    def trim(self):
        trim_beginning = 0
        contigs = self.__contigs
        orientations = self.__orientations
        while trim_beginning < len(contigs)-1 and len(contigs[trim_beginning].links[orientations[trim_beginning]]) == 1 and len(contigs[trim_beginning+1].links[1-orientations[trim_beginning+1]]) == 1 :
            trim_beginning += 1

        trim_end = 0
        while trim_end < len(contigs)-1-trim_beginning and len(contigs[-1-trim_end].links[1-orientations[-1-trim_end]]) == 1 and len(contigs[-1-trim_end-1].links[orientations[-1-trim_end-1]]) == 1 :
            trim_end += 1

        # if trim_beginning > 0 or True :
        #     print("Now trimming coiJ xco ", self)

        self.__contigs = self.__contigs[trim_beginning:len(contigs)-trim_end]
        self.__orientations = self.__orientations[trim_beginning:len(contigs)-trim_end]

        # if trim_beginning > 0 :
        #     print("Now trimming coiJ DS ", self)

#function to unzip the graph without making any assumptions
#input: the graph and the gaf file
#output: an unzipped graph
def simple_unzip(segments, names, gafFile) :

    lines = []
    print("Reading the gaf file...")
    read_GAF(gafFile, 0, 0, lines)
    print("Finished going through the gaf file.")

    #count the number of dead ends in the graph
    nbOfDeadEnds = 0
    for segment in segments :   
        for end in range(2) :
            if len(segment.links[end]) == 0 :
                nbOfDeadEnds += 1
    old_segments = segments.copy()

    #get rid of the links that are not in the gaf file
    segments = remove_unsupported_links(segments, names, lines)

    #count the number of dead ends in the graph now
    nbOfDeadEndsNow = 0
    for segment in segments :
        for end in range(2) :
            if len(segment.links[end]) == 0 :
                nbOfDeadEndsNow += 1

    #if too many dead ends where created, it means that the graph cannot be untangled properly
    # if nbOfDeadEndsNow > len(segments)/2 and nbOfDeadEndsNow > nbOfDeadEnds * 2 :
    #     print("WARNING: the graph cannot be untangled properly. That is probably because the reads are too short. The result remains valid, albeit less contiguous.")
    #     return old_segments

    on_which_paths_is_this_contig = {}
    for s in segments :
        on_which_paths_is_this_contig[s] = []

    # translate the paths in the GAF as paths in the graph
    paths = []
    for line in lines :
        # if line[0] != "SRR13128013.131014 131014 length=8539" :
        #     continue
        #split the line on '>' and '<' to get the path
        cont = re.split('[><]' , line[1].rstrip())
        orientations = "".join(re.findall("[<>]", line[1]))
        del cont[0] #because the first element is always ''
        try : #THIS SHOULD BE REMOVED IN FUTURE VERSIONS
            contigs = [segments[names[i]] for i in cont] 
        except KeyError :
            # print("fqjkldqjk pouet pouet ERROR ")
            pass

        paths.append(Path(contigs, orientations, line[0]))
        # if '53110' in cont :
        #     print("Here is a path: ", line[0], " ", paths[-1])

    #make sure all paths are coherent
    new_paths = []
    for p in paths:
        # if "53110" in str(p) :
        #     print("Here is sd ddisi Path: ", p.split_if_invalid())
        new_paths += [i for i in p.split_if_invalid()]

    paths = new_paths

    # #trim all the paths :
    # for p in paths : 
    #     p.trim()

    # print("All the paths are indexed: hccue")

    pa = 0
    for p in paths :
        co = 0
        for c in p.get_contigs() :
            on_which_paths_is_this_contig[c].append((pa, co))
            # if c.names == ['1630--1'] :
            #     print("my litt segment is here! ", p)
            co += 1
        pa += 1

    toDelete = set()
    go_on = True
    round = 0
    potentially_interesting_segments = set()
    for segment in segments :
        if len(segment.links[0]) > 1 or len(segment.links[1]) > 1 :
            potentially_interesting_segments.add(segment)
    while go_on :
        round += 1
        go_on = False
        next_potentially_interesting_segments = set()
        for segment in segments :

            if segment not in potentially_interesting_segments :
                continue
            # print("Looking icizzcce at segment : ", segment.names, " ", len(segment.links[0]), " ", len(segment.links[1]))
            segment_to_duplicate = False
            #see if it should be duplicated
            if len(segment.links[0]) > 1 or len(segment.links[1]) > 1 : 

                # print("Looking at segment ", segment.names, " ", len(segment.links[0]), " ", len(segment.links[1]))

                #list the paths going through the contig
                pairs = {}
                pair_to_paths = {}
                for p in on_which_paths_is_this_contig[segment] :
                    # if segment.names == ['119--1'] :
                    #     print("ciciuddou ", paths[p[0]], " ", p[1], " ", len(paths[p[0]])-1)
                    if len(paths[p[0]]) == 0 :
                        continue
                    path = paths[p[0]]
                    deadendleft = False
                    deadendright = False
                    contigs = path.get_contigs()
                    orientations = path.get_orientations()
                    # print("path ", p, " ", path)
                    # print(paths[5])
                    # print(paths[6])
                    if p[1] == 0 and (len(segment.links[0]) == 0 and orientations[p[1]] == 1) or (len(segment.links[1]) == 0 and orientations[p[1]] == 0) :
                        deadendleft = True
                    if p[1] == len(path)-1 and (len(segment.links[0]) == 0 and orientations[p[1]] == 0) or (len(segment.links[1]) == 0 and orientations[p[1]] == 1) :
                        deadendright = True
                    if (p[1] > 0 or deadendleft) and (p[1] < len(path)-1 or deadendright): #we want paths that go through the contigs, except if the contig is a dead end
                        #find the two links that are supported by this path
                        index_left = -2
                        if p[1] > 0 :
                            contig_left = contigs[p[1]-1]
                            end_left = orientations[p[1]-1]
                            index_left = sg.find_this_link(contig_left, end_left, segment.links[1-orientations[p[1]]], segment.otherEndOfLinks[1-orientations[p[1]]])
                            if index_left == -1 : #could happen if it was an esoteric link
                                continue
                                # print("ERROR: From semgent ", segment.names, " ", 1-orientations[p[1]] , ", did not find ", contig_left.names, " ", end_left, " among ", \
                                #       [(segment.links[1-orientations[p[1]]][i].names, segment.otherEndOfLinks[1-orientations[p[1]]][i]) for i in range(len(segment.links[1-orientations[p[1]]]))], \
                                #       " ", path)
                                # sys.exit()
                        index_right = -2
                        if p[1] < len(path)-1 :
                            contig_right = contigs[p[1]+1]
                            end_right = 1-orientations[p[1]+1]
                            index_right = sg.find_this_link(contig_right, end_right, segment.links[orientations[p[1]]], segment.otherEndOfLinks[orientations[p[1]]])
                            # print("looking for ", contig_right.names, " ", end_right, " ", contig_right.ID, " among ", \
                            #       [(segment.links[orientations[p[1]]][i].names, segment.otherEndOfLinks[orientations[p[1]]][i], segment.ID) for i in range(len(segment.links[orientations[p[1]]]))], \
                            #         " ", path, " ", index_right)
                            if index_right == -1 : #could happen if it was an esoteric link
                                continue

                        pair = (index_left, index_right)
                        # if segment.names == ['1657--1'] :
                        #     print("On contig ", segment.names, ", path ", path.name(), " supports the links towards ", segment.links[1-orientations[p[1]]][index_left].names, \
                        #           "and ", segment.links[orientations[p[1]]][index_right].names)
                        if orientations[p[1]] == 0 :
                            pair = (index_right, index_left)
                        if pair not in pairs :
                            pairs[pair] = 0
                            pair_to_paths[pair] = []
                        pairs[pair] += 1
                        pair_to_paths[pair].append(p)


                # if segment.names == ['53110'] :
                #     print("Looking at segment ", segment.names, " ", pairs, " ", [(segment.links[0][i[0]].names, segment.links[1][i[1]].names) for i in pairs])
                #     print("qldjl present on ",[paths[i[0]] for i in on_which_paths_is_this_contig[segment]])
                #     for i in on_which_paths_is_this_contig[segment] :
                #         if i[1] > 0 and i[1] < len(paths[i[0]].get_contigs()) -1 :
                #             print("idc : ", paths[i[0]])

                # if segment.names == ['130'] :
                #     print("Pairs off ccncnc : ", segment.names, " ", pairs)

                #test if the pairs are enough to support a duplication
                new_pairs = {}
                best_pair_for_each_left_link = [(-1,(-1,-1)) for i in segment.links[0]]
                best_pair_for_each_right_link = [(-1,(-1,-1)) for i in segment.links[1]]

                for p in pairs.keys() :
                    if p[0] > -1 and best_pair_for_each_left_link[p[0]][0] < pairs[p] :
                        best_pair_for_each_left_link[p[0]] = (pairs[p], p)
                    if p[1] > -1 and best_pair_for_each_right_link[p[1]][0] < pairs[p] :
                        best_pair_for_each_right_link[p[1]] = (pairs[p], p)

                #keep the pair if it is strong, or if it is the only option to avoid a dead end
                for p in pairs.keys() :
                    if pairs[p] >= 3 : 
                        new_pairs[p] = pairs[p]
                for pair in best_pair_for_each_left_link :
                    if pair[0] > 0 and (len(new_pairs)== 0 or len(segment.links[0][pair[1][0]].links[segment.otherEndOfLinks[0][pair[1][0]]]) == 1):
                        new_pairs[pair[1]] = pair[0]
                for pair in best_pair_for_each_right_link :
                    if pair[0] > 0 and (len(new_pairs)== 0 or len(segment.links[1][pair[1][1]].links[segment.otherEndOfLinks[1][pair[1][1]]]) == 1):
                        new_pairs[pair[1]] = pair[0]
                pairs = new_pairs

                # if segment.names == ['3151--1'] :
                #     print("newpais ", segment.names, " ", pairs, " ", [(segment.links[0][i[0]].names, segment.links[1][i[1]].names) for i in pairs])

                all_links_left = [False for i in range(len(segment.links[0]))]
                all_links_right = [False for i in range(len(segment.links[1]))]

                for pair in pairs.keys() :
                    if pair[0] > -1 :
                        all_links_left[pair[0]] = True
                    if pair[1] > -1 :
                        all_links_right[pair[1]] = True

                # segment_to_duplicate = all([i for i in all_links_right+all_links_left]) #to be sure not to add dead ends
                # segment_to_duplicate = True #this means that a link that is not supported by a path will be deleted
                segment_to_duplicate = any([pairs[p]>=3 for p in pairs.keys()]) #at least one pair is strong enough to support a duplication

                if segment_to_duplicate and len(pairs) > 0 :#and segment.names == ['NC_038882.1_9000_0'] :

                    #mark all the neighbors of the segment as potentially interesting
                    for neighbor in segment.links[0] :
                        next_potentially_interesting_segments.add(neighbor)
                    for neighbor in segment.links[1] :
                        next_potentially_interesting_segments.add(neighbor)

                    totalCoverage = np.sum(p for p in pairs.values())
                    for pair in pairs.keys() :
                        # if heeere :
                        #     print("On contig ", segment.names, " a pair is ", segment.links[0][pair[0]].names, " ", segment.links[1][pair[1]].names, " ", pairs[pair], " ", pair )


                        #create a new segment
                        new_coverages = [pairs[pair]/totalCoverage*segment.depth for i in range(len(segment.names)) ]
                        new_segment = sg.Segment(segment.names, segment.orientations, segment.lengths, segInsideCIGARs=segment.insideCIGARs, readCoverage=new_coverages)
                        if pair[0] >= 0 :
                            sg.add_link(segment.links[0][pair[0]], segment.otherEndOfLinks[0][pair[0]], new_segment, 0, CIGAR = segment.CIGARs[0][pair[0]])
                        if pair[1] >= 0 :
                            sg.add_link(segment.links[1][pair[1]], segment.otherEndOfLinks[1][pair[1]], new_segment, 1, CIGAR = segment.CIGARs[1][pair[1]])

                        # print("Pair to duplicate : ", pair, " ", pair_to_paths[pair])
                        for pa in pair_to_paths[pair] :
                            # if pa[0] == 6 :
                            #     print("1055 is heer e ", paths[pa[0]], " ", segment.ID, " ", new_segment.ID, " ", [i.ID for i in paths[pa[0]].get_contigs()])
                            paths[pa[0]].replace(segment, new_segment, pa[1])

                        on_which_paths_is_this_contig[new_segment] = pair_to_paths[pair]

                        segments.append(new_segment)
                        next_potentially_interesting_segments.add(new_segment)
                        potentially_interesting_segments.add(new_segment)

                    pa = 0
                    # for p in paths :
                    for p in on_which_paths_is_this_contig[segment] :
                        path = paths[p[0]]
                        if (p[1] == 0 and deadendleft) or (p[1] == len(path)-1 and deadendright) :
                            continue
                        # if segment.names == ['edge_389@0_0_2'] or True:
                        # # if pa == 31624 :
                        #     print("On contig ", segment.names, " the pairs are ", [(segment.links[0][pair[0]].names, segment.links[1][pair[1]].names) for pair in pairs.keys()] )
                        # #     print("path ", path, " ", segment.names, " ", pairs, " ", pa)
                        path.cancel(segment)
                        pa += 1

                    segment.cut_all_links()
                    toDelete.add(segment)
                    on_which_paths_is_this_contig[segment] = []

                    go_on = True
            
        potentially_interesting_segments = next_potentially_interesting_segments

    delIdx = 0
    while delIdx < len(segments) :
        if segments[delIdx] in toDelete :
            del segments[delIdx]
            del on_which_paths_is_this_contig[segments[delIdx]]
        else :
            delIdx += 1
            
    # print("NOT DETACHING TIPS")
    detach_and_destroy_tips(segments)

    return segments

#function that removes the links that are not supported by any path
#input : a list of segments and of paths
#output : the same list of segments, but with the links that are not supported by any path removed
def remove_unsupported_links(segments, names, lines):

    #inventory of the links in the lines
    links = set()
    for line in lines :
        cont = re.split('[><]' , line[1].rstrip())
        orientations = "".join(re.findall("[<>]", line[1]))
        del cont[0] #because the first element is always ''
        try : #TO BE REMOVE IN FUTURE VERSIONS
            contigs = [segments[names[i]] for i in cont] 
        except KeyError :
            # print([names[i] for i in cont])
            continue

        for i in range(len(contigs)-1) :
            links.add((contigs[i], "<>".index(orientations[i]), contigs[i+1], "><".index(orientations[i+1])))
            links.add((contigs[i+1], "><".index(orientations[i+1]), contigs[i], "<>".index(orientations[i])))

    #remove the links that are not supported by any path
    toRemove = set()
    for segment in segments :
        for end in range(2) :
            for n, neighbor in enumerate(segment.links[end]) :
                if (segment, end, neighbor, segment.otherEndOfLinks[end][n]) not in links :
                    toRemove.add((segment, end, neighbor, segment.otherEndOfLinks[end][n]))

    for segment, end, neighbor, otherEnd in toRemove :
        sg.delete_link(segment, end, neighbor, otherEnd, warning=False)

    return segments

#function that detach short tips
#input : a list of segments
#output : the same list of segments, but with the short tips detached
def detach_and_destroy_tips(segments):
                                              
    #detach short dead ends
    changes = True
    max_tip_length = 1000
    contig_to_delete = set()

    for s, seg in enumerate(segments):
        
        for end in range(2) :

            if len(seg.links[end]) > 1 : #one of the branch may be a short dead end 
                already_visited = set()   
                extended_lengths = [extended_length(seg.links[end][i], seg.otherEndOfLinks[end][i], max_tip_length*5, 50, already_visited) for i in range(len(seg.links[end]))]
                max_length = max(extended_lengths)
                toDelete = set()
                for n in range(len(seg.links[end])) :
                    if 5*extended_lengths[n] < max_length and max_length > 1000 : #the branch is a short dead end
                        toDelete.add((seg, end, seg.links[end][n], seg.otherEndOfLinks[end][n]))
                        contig_to_delete.add(seg.links[end][n].ID)
                for seg, end, neighbor, otherEnd in toDelete :
                    sg.delete_link(seg, end, neighbor, otherEnd, warning=False)

    #destroy the contigs that are in contig_to_delete
    delIdx = 0
    while delIdx < len(segments) :
        if segments[delIdx].ID in contig_to_delete :
            segments[delIdx].cut_all_links()
            del segments[delIdx]
        else :
            delIdx += 1

    return segments
                 

#a function returning the longest path you can get with neighbors of neighbors of neighbors of... (up to threshold) 
def extended_length(segment, end, thresholdLength, thresholdContigs, already_visited) :
        
    if thresholdContigs == 0 or thresholdLength <= 0 :
        return segment.length
    
    already_visited.add((segment, end))
    
    #start by looking down the longest contig, it will be fastest
    longestContig = [i for i in range(len(segment.links[1-end]))]
    longestContig.sort(key= lambda x : segment.links[1-end][x].length, reverse = True)
    
    #print("Longest contigs : ", longestContig, [segment.links[1-end][i].length for i in longestContig])
    
    maxLength = 0
    for n in longestContig :
        neighbor = segment.links[1-end][n]
        if (neighbor, segment.otherEndOfLinks[1-end][n]) in already_visited :
            return thresholdLength
        l = extended_length(neighbor, segment.otherEndOfLinks[1-end][n], thresholdLength-segment.length, thresholdContigs-1, already_visited)
        if l > maxLength :
            maxLength = l
        
    return maxLength + segment.length



