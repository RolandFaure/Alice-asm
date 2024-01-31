import numpy as np
import random
import sys
import re

#a segment is a supercontig
class Segment:

    def __init__(self, segNamesOfContig, segOrientationOfContigs, segLengths, segInsideCIGARs =None, segLinks = [[],[]], segOtherEndOfLinks = [[],[]], segCIGARs = [[],[]], lock = False, HiCcoverage = 0, readCoverage = []):
         
        if len(segLinks[0]) != len(segOtherEndOfLinks[0]) or len(segLinks[1]) != len(segOtherEndOfLinks[1]) :
            print('ERROR in the links while initializing a segment')
            return 0
        
        if any(i!=0 and i!=1 for i in segOtherEndOfLinks[0]) or any(i!=0 and i!=1 for i in segOtherEndOfLinks[1]):
            print('ERROR in the links while initializing a segment')
            return 0
        
        if len(segNamesOfContig) != len(segOrientationOfContigs) :
            print('ERROR in initializing the orientations of contigs within a segment')
            return 0
        
        if segInsideCIGARs == None :
            segInsideCIGARs = ['*' for i in range(len(segNamesOfContig)-1)]
            
        if segCIGARs == [[],[]] :
            segCIGARs = [['*' for i in segLinks[0]],['*' for i in segLinks[1]]]
        
        self._id = random.random() #a random number identifying the segment
        self._HiCcoverage = HiCcoverage #to keep in mind how many HiC contacts this segment has in total
        
        #this group of attributes are linked arrays : element n in one corresponds with element n in the other. Therefore they shouldn't be modified independantly
        self._namesOfContigs = segNamesOfContig.copy() #names are strings with which sequences are described in the GFA
        self._orientationOfContigs = segOrientationOfContigs.copy() #1 being '+' orientation, 0 the '-' orientation
        self._lengths = segLengths.copy()
        self._reads = [[] for i in range(len(segNamesOfContig))] #to keep in mind the reads that are associated with each subcontig
        self._sequences = [None for i in range(len(segNamesOfContig))] #to keep in mind the sequences of each subcontig, and possibly repolish them
        self._insideCIGARs = segInsideCIGARs.copy()

        if readCoverage != [] :
            self._depths = readCoverage.copy() #to keep in mind the read coverage of the contigs
        else :
            self._depths = [1 for i in range(len(self._lengths))]
            
        self._copiesOfContigs = [-1]*len(segNamesOfContig) #this is used exclusively while exporting, to indicate which copy of which contig is in the segment (copy 0 / copy 1 / copy 2 ...)
        
        #this group of attribute are linked arrays : one should never be modified without the others 
        #They are sorted by ID of the neighbor : important for handling quickly big nodes
        lists_keyed = [[(segLinks[0][i], segOtherEndOfLinks[0][i], segCIGARs[0][i]) for i in range(len(segLinks[0]))], [(segLinks[1][i], segOtherEndOfLinks[1][i], segCIGARs[1][i]) for i in range(len(segLinks[1]))]]
        lists_keyed[0].sort(key = lambda x: x[0].ID)
        lists_keyed[1].sort(key = lambda x: x[0].ID)
        self._links = [[i[0] for i in lists_keyed[0]], [i[0] for i in lists_keyed[1]]] #two lists of segments with which the segment is linked, at the left end and at the right end
        self._otherEndOfLinks = [[i[1] for i in lists_keyed[0]], [i[1] for i in lists_keyed[1]]] #for each link, indicates the side of the other segment on which the link arrives
        self._CIGARs = [[i[2] for i in lists_keyed[0]], [i[2] for i in lists_keyed[1]]] #for each link, indicates the CIGAR string found in the GFA
        
        self._trim = [0,0]

        self._freezed = [False, False] #do not duplicate from one end if frozen at this end
        self._locked = lock #That is to duplicate a contig only once in each merge_contigs
          
    # getters
    def get_id(self):
        return self._id
    
    # def __hash__(self):
    #     return self._id
    
    def get_orientations(self):
        return self._orientationOfContigs
    
    def get_insideCIGARs(self):
        return self._insideCIGARs
    
    def get_links(self):
        return self._links
    
    def get_otherEndOfLinks(self):
        return self._otherEndOfLinks
    
    def get_CIGARs(self):
        return self._CIGARs
    
    def get_lengths(self):
        return self._lengths
    
    def get_length(self):
        return np.sum(self._lengths)
    
    def get_namesOfContigs(self):
        return self._namesOfContigs
    
    def get_copiesOfContigs(self):
        return self._copiesOfContigs
    
    def get_sequences(self):
        return self._sequences
    
    def get_freezed(self):
        return self._freezed
    
    def get_locked(self):
        return self._locked
    
    def get_coverage(self):
        return self._HiCcoverage
        
    def get_depths(self):
        return self._depths
    
    def get_reads(self):
        return self._reads
    
    def get_depth(self):
        sumdepth = 0
        sumlength = 1 #1 and not 0 to be sure not to divide by 0
        for i in range(len(self._depths)) :
            sumdepth += self._depths[i]*self._lengths[i]
            sumlength += self._lengths[i]
            
        return sumdepth / sumlength
    
    def get_trim(self):
        return self._trim
    
    def full_name(self) :
        return '_'.join([self._namesOfContigs[i]+'-'+str(self._copiesOfContigs[i]) for i in range(len(self._namesOfContigs))])
    
    def print_complete(self):
        print(self._namesOfContigs, [s.names for s in self._links[0]], \
              [s.names for s in self._links[1]])
    
    # setters 
    def set_copiesNumber(self, copiesNumberForNow):
        for c, contig in enumerate(self._namesOfContigs) :
            if contig in copiesNumberForNow :
                self._copiesOfContigs[c] = copiesNumberForNow[contig]
                copiesNumberForNow[contig] += 1
            else :
                self._copiesOfContigs[c] = 0
                copiesNumberForNow[contig] = 1
    
    def set_coverage(self, newCoverage) :
        self._HiCcoverage = newCoverage
    
    def set_sequences(self, newSequences) :
        self._sequences = newSequences

    def set_orientation(self, index, newOrientation) :
        self._orientationOfContigs[index] = newOrientation
    
    def add_read(self, nameOfContig, nameOfRead):
        # print("adding ", nameOfRead, " to ", nameOfContig, " of ", self._namesOfContigs)
        self._reads[self._namesOfContigs.index(nameOfContig)].append(nameOfRead)

    def set_CIGAR(self, end, otherContig, otherEnd, newCIGAR) :
        for l in range(len(self.links[end])) :
            if self.links[end][l].ID == otherContig.ID and self.otherEndOfLinks[end][l] == otherEnd :
                self.CIGARs[end][l] = newCIGAR
                return 0
        return -1
    
    def freeze(self, endOfSegment): 
        self._freezed[endOfSegment] = True
  
    def freezeNode(self, endOfSegment):
        self._freezed[endOfSegment] = True
        for n, neighbor in enumerate(self._links[endOfSegment]) :
            neighbor.freeze(self._otherEndOfLinks[endOfSegment][n])
        
    def unfreeze(self):
        self._freezed = [False, False]
        
    def set_locked(self, b):
        self._locked = b
        
    def lockNode(self, endOfSegment):
        self._locked = True
        for i in self._links[endOfSegment]:
            i.locked = True

    def set_trim(self, trim) :
        self._trim = trim
            
    def set_id(self, newID) : #few cases where that is useful
        self._id = newID 
        
    def divide_depths(self, n) : #when duplicating a segment, you need to lower the coverage of all replicas
        for i in range (len(self._depths)) :
            self._depths[i] /= n
            
    def multiply_end_depths(self, n, end, numberOfContigs) : #this fucntion is useful when merging a wrongly duplicated dead end
        for co in range(numberOfContigs) :
           self._depths[end*(len(self._depths)-1) + (-2*end+1)*co ]
            
    def length1(self): #use this function to set length of a segment to 1 (instead of 0, mostly)
        self._lengths = [max(1, i) for i in self._lengths]
                    
    # properties
    
    ID = property(get_id, set_id)
    HiCcoverage = property(get_coverage, set_coverage)
    depths = property(get_depths)
    depth = property(get_depth)
    length = property(get_length)
    
    names = property(get_namesOfContigs)
    orientations = property(get_orientations)
    lengths = property(get_lengths)
    insideCIGARs = property(get_insideCIGARs)
    
    copiesnumber = property(get_copiesOfContigs)
    
    links = property(get_links)
    otherEndOfLinks = property(get_otherEndOfLinks)
    CIGARs = property(get_CIGARs)
    
    freezed = property(get_freezed) #segment is freezed if comparison of links was inconclusive
    locked = property(get_locked, set_locked) #segment is locked if it was already duplicated in merge_contigs and that it should not be for a second time

    #other functions that handle segments
    
    #function which goals is to return the intensity of Hi-C contacts between another segment and this one
    def interaction_with_contigs(self, segment, interactionMatrix, names, copiesnumber = None, commonContigs = set(), bestSignature = 1000):
        
        if copiesnumber == None :
            copiesnumber = [1 for i in interactionMatrix]
            
        absoluteScore = 0
        relativeScore = 0
        
        # orientation = -1 # if supercontig is directly linked to the candidates, then this variable tells us by which end
        
        # if segment in self.links[0]:
        #     orientation = 0
        # elif segment in self.links[1]:
        #     orientation = 1
        # else: 
        #     print('ERROR : trying to compute an interaction but the contigs do not touch each other')
        #     print('Looking for ', self._namesOfContigs, ' from ', segment.names)
        #     return 0, 0, 1
            
        depth = 1
        #first compute interactions with self
        for co, contig in enumerate(self.names) :
                 
            for c, contigInSegment in enumerate(segment.names):
                
                if contig not in commonContigs and copiesnumber[contigInSegment] <= bestSignature:
                    
                    absoluteScore += interactionMatrix[names[contig], names[contigInSegment]]
                    relativeScore += interactionMatrix[names[contig], names[contigInSegment]]
                else:
                    absoluteScore += interactionMatrix[names[contig], names[contigInSegment]]
                                        
        # #now compute the interaction with neighbors of self 
        # endOfSegment = 1-orientation
        # for neighbor in self.links[endOfSegment] :
        #     for co, contig in enumerate(neighbor.names) :
        #         for c, contigInSegment in enumerate(segment.names):
                
        #             if contig not in commonContigs and copiesnumber[contigInSegment] <= bestSignature:
                        
        #                 depth = 2
                        
        #                 absoluteScore += interactionMatrix[names[contigInSegment],names[contig]]
        #                 relativeScore += interactionMatrix[names[contigInSegment],names[contig]]
        #             else:
        #                 absoluteScore += interactionMatrix[names[contigInSegment],names[contig]]
                
            
        return absoluteScore, relativeScore, depth
    
    def add_link_from_GFA(self, GFAline, names, segments, leftOrRight) : #leftOrRight = 0 when the segment is at the beginning of a link (left of a GFA line), 1 otherwise
        
        l = GFAline.strip('\n').split('\t')
 
        if len(l) < 5 :
            print('ERROR : expected at least 5 fields in line ', GFAline)
        
        if l[0] != 'L':
            print('ERROR : trying to add a link from a GFA line that does not start with "L"')
        
        else :
            o1,o2 = -1, -1
            
            if l[2] == '-':
                o1 = 0
            elif l[2] == '+' :
                o1 = 1
                
            if l[4] == '-':
                o2 = 0
            elif l[4] == '+' :
                o2 = 1
                
            if o1 == -1 or o2 == -1 :
                print('ERROR while creating a link : orientations not properly given.')
                print('Problematic line : ', GFAline)   

            
            if leftOrRight == 0 and o1 == 0:

                #then comes a little test to see if the link has already been added (for example if their are several lines in the GFA describing the same link). An 'or' condition to allow it when an end of link is linked with itself
                if find_this_link(segments[names[l[3]]], 1-o2, self._links[0], self._otherEndOfLinks[0]) == -1 or (segments[names[l[3]]].ID == self._id and 1-o2 == 0):

                    index = index_at_which_new_link_should_be_inserted(segments[names[l[3]]], self._links[0], 1-o2 ,self._otherEndOfLinks[0])
                    self._links[0].insert(index, segments[names[l[3]]])
                    self._otherEndOfLinks[0].insert(index, 1-o2)
                    if len(l) > 5 :
                        self._CIGARs[0].insert(index, l[5])
                    else :
                        self._CIGARs[0].insert(index, '*')
                    
            elif leftOrRight == 0 and o1 == 1 :
                if find_this_link(segments[names[l[3]]], 1-o2, self._links[1], self._otherEndOfLinks[1]) == -1 or (segments[names[l[3]]].ID == self._id and 1-o2 == 1):
                    
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[3]]], self._links[1],  1-o2 ,self._otherEndOfLinks[1])
                    
                    self._links[1].insert(index, segments[names[l[3]]])
                    self._otherEndOfLinks[1].insert(index, 1-o2)
                    if len(l) > 5 :
                        self._CIGARs[1].insert(index, l[5])
                    else :
                        self._CIGARs[1].insert(index, '*')
                
            elif leftOrRight == 1 and o2 == 1 :
                if find_this_link(segments[names[l[1]]], o1, self._links[0], self._otherEndOfLinks[0]) == -1 or (segments[names[l[1]]].ID == self._id and o1 == 0) :
                    
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[1]]], self._links[0],  o1 ,self._otherEndOfLinks[0])
                    self._links[0].insert(index, segments[names[l[1]]])
                    self._otherEndOfLinks[0].insert(index, o1)
                    if len(l) > 5 :
                        self._CIGARs[0].insert(index, l[5])
                    else :
                        self._CIGARs[0].insert(index, '*')
                    
            elif leftOrRight == 1 and o2 == 0 :
                if find_this_link(segments[names[l[1]]], o1, self._links[1], self._otherEndOfLinks[1]) == -1 or (segments[names[l[1]]].ID == self._id and o1 == 1):
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[1]]], self._links[1],  o1 ,self._otherEndOfLinks[1])
                    self._links[1].insert(index, segments[names[l[1]]])
                    self._otherEndOfLinks[1].insert(index, o1)
                    if len(l) > 5 :
                        self._CIGARs[1].insert(index, l[5])
                    else :
                        self._CIGARs[1].insert(index, '*')
            
            else :
                print('ERROR while trying to add a new link from the gfa : could not locate a correct name')
        
    #this adds the end of a links, but only on this segment, not on the other end
    def add_end_of_link(self, endOfSegment, segment2, endOfSegment2, CIGAR = '*'):
        
        #print('A', len(segment2.otherEndOfLinks[1]), len(segment2.links[1]), len(segment2.CIGARs[1]))
        #print(self._namesOfContigs, segment2.names)
        index = index_at_which_new_link_should_be_inserted(segment2, self._links[endOfSegment], endOfSegment2, self._otherEndOfLinks[endOfSegment])

        self._links[endOfSegment].insert(index, segment2)
        #print('B', len(segment2.otherEndOfLinks[1]), len(segment2.links[1]), len(segment2.CIGARs[1]))

        self._otherEndOfLinks[endOfSegment].insert(index, endOfSegment2)
        self._CIGARs[endOfSegment].insert(index, CIGAR)        

    #this function is useful for rerouting around big nodes. It adds end of links more efficiently than if done individually, and checks for doubles
    def add_a_bunch_of_end_of_links(self, endOfSegment, listOfSegmentsToAdd, listOfEndOfSegmentsToAdd, CIGARsToAdd) :
 
        #All the list of segments being sorted by ID, this is done on the basis of the merging of merge sort.
        newLinks = []
        newEndOfLinks = []
        newCIGARs = []

        indexOfSegmentToAdd = 0
        indexOfLinks = 0
        
        while indexOfSegmentToAdd < len(listOfSegmentsToAdd) and indexOfLinks < len(self._links[endOfSegment]) :
            
            if listOfSegmentsToAdd[indexOfSegmentToAdd].ID < self._links[endOfSegment][indexOfLinks].ID :
                newLinks.append(listOfSegmentsToAdd[indexOfSegmentToAdd])
                newEndOfLinks.append(listOfEndOfSegmentsToAdd[indexOfSegmentToAdd])
                newCIGARs.append( CIGARsToAdd[indexOfSegmentToAdd] )
                indexOfSegmentToAdd += 1
                
            elif listOfSegmentsToAdd[indexOfSegmentToAdd].ID > self._links[endOfSegment][indexOfLinks].ID :
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks] )
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks] )
                indexOfLinks += 1
                
            #elsewhise the link goes toward the same segment, but maybe not the same end of this segment
            elif listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] < self._otherEndOfLinks[endOfSegment][indexOfLinks] :
                newLinks.append( listOfSegmentsToAdd[indexOfSegmentToAdd] )
                newEndOfLinks.append( listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] )
                newCIGARs.append( CIGARsToAdd[indexOfSegmentToAdd] )
                indexOfSegmentToAdd += 1
                
            elif listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] > self._otherEndOfLinks[endOfSegment][indexOfLinks] :
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks] )
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks] )
                indexOfLinks += 1
                
            else : #it means that there is twice the same link, so add it only once
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks])
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks])
                indexOfLinks += 1
                indexOfSegmentToAdd += 1
                
        if indexOfSegmentToAdd >= len(listOfEndOfSegmentsToAdd) :
            newLinks += self._links[endOfSegment][indexOfLinks:]
            newEndOfLinks += self._otherEndOfLinks[endOfSegment][indexOfLinks:]
            newCIGARs += self._CIGARs[endOfSegment][indexOfLinks:]
            
        elif indexOfLinks >= len(self._links[endOfSegment]) :
            newLinks += listOfSegmentsToAdd[indexOfSegmentToAdd:]
            newEndOfLinks += listOfEndOfSegmentsToAdd[indexOfSegmentToAdd:]
            newCIGARs += CIGARsToAdd[indexOfSegmentToAdd:]
        
        self._links[endOfSegment] = newLinks
        self._otherEndOfLinks[endOfSegment] = newEndOfLinks
        self._CIGARs[endOfSegment] = newCIGARs

    def remove_end_of_link(self, endOfSegment, segmentToRemove, endOfSegmentToRemove = None, warning = True): #endOfSegmentToRemove is there in case there exists two links between self[endOfSegment] and segment to remove. Needed for extra security
        
        #first determine the index of the segment to remove
        #print('Removing ', segmentToRemove.names, endOfSegmentToRemove, ' from ', self._namesOfContigs)
        #print('Among these links :', [i.names for i in self._links[endOfSegment]], self._otherEndOfLinks[endOfSegment])
        index = find_this_link(segmentToRemove, endOfSegmentToRemove, self._links[endOfSegment], self._otherEndOfLinks[endOfSegment], warning = warning)
        #index = self._links[endOfSegment].index(segmentToRemove)
   
        #then remove the end of unwanted link in all attributes
        if index != -1 :
            del self._links[endOfSegment][index]
            del self._otherEndOfLinks[endOfSegment][index]
            del self._CIGARs[endOfSegment][index]
            return True
        elif index == -1 and warning:
             print('Trying unsuccesfully to remove ', segmentToRemove.names, ' from ', self._namesOfContigs)
             return False
        
    #returns two contigs, equal to this contig but split at axis, corresponding to the number of contigs left of the junction
    def break_contig(self, axis) :
        
        newSegment1 = Segment(self._namesOfContigs[:axis], self._orientationOfContigs[:axis], self._lengths[:axis], self._insideCIGARs[:axis-1], [self._links[0], []], [self._otherEndOfLinks[0], []], [self._CIGARs[0], []], readCoverage = self._depths[:axis])
        
        newSegment2 = Segment(self._namesOfContigs[axis:], self._orientationOfContigs[axis:], self._lengths[axis:], self._insideCIGARs[axis:], [[], self._links[1]], [[], self._otherEndOfLinks[1]], [[], self._CIGARs[1]], readCoverage = self._depths[axis:])
        
        return newSegment1, newSegment2
     
    #function to be used on small loops only
    def flatten(self, replicas) :
        if self not in self._links[0] :
            
            print('ERROR : in segment.flatten, trying to flatten something that is not a loop')
            
        else :
            
            for i in range(len(self._depths)) :
                self._depths[i] /= replicas+1
            
            newName = self._namesOfContigs.copy()
            newOrientations = self._orientationOfContigs.copy()
            newLengths = self._lengths.copy()
            newinsideCIGARs = self._insideCIGARs.copy()
            newCopies = self._copiesOfContigs.copy()
            newDepths = self._depths.copy()
            for i in range(replicas) :
                newName += self._namesOfContigs
                newOrientations += self._orientationOfContigs
                newLengths += self._lengths
                newCopies += self._copiesOfContigs
                newinsideCIGARs += [self._CIGARs[0][self._links[0].index(self)]] + self._insideCIGARs
                newDepths += self._depths
            
            self._namesOfContigs = newName
            self._orientationOfContigs = newOrientations
            self._lengths = newLengths
            self._copiesOfContigs = newCopies
            self._insideCIGARs = newinsideCIGARs
            self._depths = newDepths
            
           # print('In segment.flatten : ', self._namesOfContigs, self._insideCIGARs, [self._CIGARs[0][self._links[0].index(self)]])
            #print('Links before any removal, ', [i.names for i in self._links[0]], '\n')
            self.remove_end_of_link(0, self)
            self.remove_end_of_link(1, self)
  
    #a function to delete all the links of a segment (typically before deleting it)
    def cut_all_links(self) :
        
        for end in range(2) :
            
            for n, neighbor in enumerate(self._links[end]) :
                
                otherEnd = self._otherEndOfLinks[end][n]
                neighbor.remove_end_of_link(otherEnd, self, end)
        
        self._links = [[],[]]
        self._otherEndOfLinks = [[],[]]
        self._CIGARs = [[],[]]
        
#This function is OUTSIDE the class. It takes two segments and the end of the first segment which is linked to the second. It appends a merged contig to the listOfSegments, without modifying the two inputed segments
def merge_two_segments(segment1, endOfSegment1, segment2, listOfSegments):
    
    if segment1.links[endOfSegment1].count(segment2) > 1 : #this means a loop
        return 0
      
    # creating a new segment
    orientation1 = endOfSegment1*2-1
    endOfSegment2 = segment1.otherEndOfLinks[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    orientation2 = -2*endOfSegment2+1
    
    orientationOfContigs1 = segment1.orientations
    orientationOfContigs2 = segment2.orientations
    
    if orientation1 == -1 : #then change the orientation of all the contigs within the segment, since the segment will be mirrored in the new supersegment
        orientationOfContigs1 = [1-i for i in orientationOfContigs1]
    if orientation2 == -1 :
        orientationOfContigs2 = [1-i for i in orientationOfContigs2]
        
    CIGAR = segment1.CIGARs[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    
    newSegment = Segment(segment1.names[::orientation1] + segment2.names[::orientation2],\
                                orientationOfContigs1[::orientation1]+orientationOfContigs2[::orientation2],\
                                segment1.lengths[::orientation1]+segment2.lengths[::orientation2], \
                                segment1.insideCIGARs[::orientation1] + [CIGAR] + segment2.insideCIGARs[::orientation2],\
                                segLinks = [segment1.links[1-endOfSegment1], \
                                segment2.links[1-endOfSegment2]], \
                                segOtherEndOfLinks = [segment1.otherEndOfLinks[1-endOfSegment1], \
                                segment2.otherEndOfLinks[1-endOfSegment2]],\
                                segCIGARs = [segment1.CIGARs[1-endOfSegment1], \
                                segment2.CIGARs[1-endOfSegment2]],
                                lock = True,
                                HiCcoverage = segment1.HiCcoverage + segment2.HiCcoverage,
                                readCoverage = segment1.depths[::orientation1] + segment2.depths[::orientation2])
            
    listOfSegments.append(newSegment)
    
    #building the other end of links with the new segment
    self_loop_CIGAR = ''
    for n, neighbor in enumerate(newSegment.links[0]) :
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[0][n], newSegment, 0, CIGAR = newSegment.CIGARs[0][n])
        if neighbor.ID == segment2.ID and newSegment.otherEndOfLinks[0][n] == 1 - endOfSegment2 : #check if the new contig should loop back on itself
            self_loop_CIGAR = newSegment.CIGARs[0][n]

    for n, neighbor in enumerate(newSegment.links[1]) :
        #print(len(newSegment.otherEndOfLinks[1]), len(newSegment.links[1]), len(newSegment.CIGARs[1]))
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[1][n], newSegment, 1, CIGAR = newSegment.CIGARs[1][n])
        
    if self_loop_CIGAR != '' :
        add_link(newSegment, 0, newSegment, 1, self_loop_CIGAR)

#function creating a link between two ends of contigs, OUTSIDE of the class
def add_link(segment1, end1, segment2, end2, CIGAR = '*'):
    segment1.add_end_of_link(end1, segment2, end2, CIGAR)
    segment2.add_end_of_link(end2, segment1, end1, CIGAR)
    
def delete_link(segment1, end1, segment2, end2, warning = True) :
    success1 = segment1.remove_end_of_link(end1, segment2, end2, warning = warning)
    success2 = segment2.remove_end_of_link(end2, segment1, end1, warning = warning)
    return success1 and success2
           
def compute_copiesNumber(listOfSegments):
    cn = {}
    for s in listOfSegments :
        s.set_copiesNumber(cn)
        
    return cn

#returns the position of the link pointing towards segment and its endOfSegment
def find_this_link(segment, endOfSegment, listOfLinks, listOfEndsOfLinks, warning = False) :

    lo = 0
    hi = len(listOfLinks)

    while lo < hi:
        mid = (lo+hi)//2
        if segment.ID < listOfLinks[mid].ID:
            hi = mid
        elif segment.ID > listOfLinks[mid].ID:
            lo = mid+1
        else :
            #print('Found : ', endOfSegment , listOfEndsOfLinks[mid])
            if endOfSegment == None :
                return mid
                
            elif endOfSegment == listOfEndsOfLinks[mid] :
                return mid
                
            elif endOfSegment > listOfEndsOfLinks[mid] :
                mid += 1
                while mid < len(listOfLinks) and listOfLinks[mid].ID == segment.ID :
                    if endOfSegment == listOfEndsOfLinks[mid] :
                        return mid
                    mid += 1
                    
                break
                
            elif endOfSegment < listOfEndsOfLinks[mid] :
                mid -= 1
                while mid >= 0 and listOfLinks[mid].ID == segment.ID :
                    if endOfSegment == listOfEndsOfLinks[mid] :
                        return mid
                    mid -= 1
                break
    
    if not warning :
        return -1
        
    print('In find_this_link : did not find the link')
    #print([[listOfLinks[se].names, listOfEndsOfLinks[se]] for se in range(len(listOfLinks))])
    print('Did not find ', segment.names , endOfSegment, ' among ', [i.names for i in listOfLinks], listOfEndsOfLinks)
    
    return -1

#returns the index at which a segment should be inserted in a list sorted by ID : useful because links[0] and links[1] are kept sorted at all times
def index_at_which_new_link_should_be_inserted(segment, listOfSegments, endOfLink, listOfEndOfLinks) :
    lo = 0
    hi = len(listOfSegments)

    while lo < hi:
        mid = (lo+hi)//2
        if segment.ID < listOfSegments[mid].ID or (segment.ID == listOfSegments[mid].ID and endOfLink < listOfEndOfLinks[mid]):
            hi = mid
        else:
            lo = mid+1
            
    while lo < len(listOfSegments) and listOfSegments[lo].ID == segment.ID and endOfLink == listOfEndOfLinks[lo] :
        lo += 1
    return lo

def check_if_all_links_are_sorted(listOfSegments) :
    
    for segment in listOfSegments :
        for endOfSegment in range(2) :
            for n in range(len(segment.links[endOfSegment])-1) :
                if segment.links[endOfSegment][n].ID > segment.links[endOfSegment][n+1].ID :
                    print('Problem in the links of ', segment.names, ' : ', [s.ID for s in segment.links[endOfSegment]])
    
    for segment in listOfSegments :
        for endOfSegment in range(2) :
            for n, neighbor in enumerate(segment.links[endOfSegment]) :
                if segment not in neighbor.links[segment.otherEndOfLinks[endOfSegment][n]] :
                    print('Non-reciprocal links : ', segment.names, segment.ID, neighbor.names, neighbor.ID)
                
#funtion to delete links that are present twice in the graph (often because they are present twice in the gfa)
def delete_links_present_twice(segments):
    

    for segment in segments :
        toBeRemoved = []
        for endOfSegment in range(2) :
            
            for n1 in range(len(segment.links[endOfSegment])-1) :
                
                for n2 in range(n1+1, len(segment.links[endOfSegment])) :
                    
                    if segment.links[endOfSegment][n1].ID == segment.links[endOfSegment][n2].ID and segment.otherEndOfLinks[endOfSegment][n1] == segment.otherEndOfLinks[endOfSegment][n2] and segment.links[endOfSegment][n1].ID != segment.ID:
                        
                        segment.links[endOfSegment][n2].remove_end_of_link(segment.otherEndOfLinks[endOfSegment][n2], segment, endOfSegment)
                        toBeRemoved += [[endOfSegment, segment.links[endOfSegment][n2], segment.otherEndOfLinks[endOfSegment][n2]]]

        for r in toBeRemoved :
            segment.remove_end_of_link(r[0], r[1], r[2])

# input : one supercontig to be joined with a neighbor at one end
# output : actualized listOfSegments with the two contigs merged
def merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments):

    if len(segment.links[endOfSegment]) != 1:
        print("ERROR : trying to merge simply two contigs that cannot be merged simply")
        return -1, -1

    neighbor = segment.links[endOfSegment][0]
    endOfSegmentNeighbor = segment.otherEndOfLinks[endOfSegment][0]

    if len(neighbor.links[endOfSegmentNeighbor]) != 1 :
        print('ERROR : trying to merge simply two contigs that cannot be merged simply')
        return -1,-1
        
    if neighbor == segment :  # then do not merge a contig with itself
        return -1, -1

    # add the new segment
    merge_two_segments(segment, endOfSegment, neighbor, listOfSegments)
    

    # delete links going towards the two ex-segments
    otherEnd = 1 - endOfSegment
    otherEndNeighbor = 1 - endOfSegmentNeighbor

    for i, n in enumerate(segment.links[otherEnd]) :
        n.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)
        
    for i, n in enumerate(neighbor.links[otherEndNeighbor]) :
        #print('Removing ', neighbor.names, ' from ', n.names, ' and adding the new contig',listOfSegments[-1].names, ' at end ', neighbor.otherEndOfLinks[otherEndNeighbor][i])
        n.remove_end_of_link(neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor, otherEndNeighbor)

    # delete the ex-segments    
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)
    
    return listOfSegments

# a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):

    goOn = True
    while goOn:
        goOn = False
        for segment in listOfSegments:

            alreadyDidThisOne = False # if the segment is deleted when looking at its first end, you don't want it to look at its other end, since it does not exist anymore
            for endOfSegment in range(2):
                
                if not alreadyDidThisOne:
                    if len(segment.links[endOfSegment]) == 1 and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]])== 1 :  # then merge
                        alreadyDidThisOne = True
                        if segment.ID != segment.links[endOfSegment][0].ID:
                            goOn = True
                            listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)


    return listOfSegments

# input : a list of segments
# output : a list of segments, with the segments that can be duplicated duplicated
def duplicate_contigs(segments):

    print("Duplicating contigs")

    continueDuplication = True
    while continueDuplication:
        continueDuplication = False
        toDelete = []
        alreadyDuplicated = set()
        originalLength = len(segments)
        for se in range (originalLength) :
            contig = segments[se]
            #check if the segment can be duplicated by one of its end
            for end in range(2):
                if contig.ID not in alreadyDuplicated:

                    if len(contig.links[end]) > 1 \
                        and all([len(contig.links[end][i].links[contig.otherEndOfLinks[end][i]]) == 1 for i in range(len(contig.links[end]))]) \
                        and (contig.depth > 0.7*np.sum([contig.links[end][i].depth for i in range(len(contig.links[end]))]) or contig.length < 1000) \
                        and contig.ID not in [i.ID for i in contig.links[end]] :

                        #then duplicate the segment !
                        alreadyDuplicated.add(contig.ID)
                        numberofcopies = len(contig.links[end])
                        toDelete += [se]
                        totalNeighborCoverage = np.sum([contig.links[end][i].depth for i in range(len(contig.links[end]))])
                        if totalNeighborCoverage == 0:
                            totalNeighborCoverage = 1
                        for n in range(numberofcopies):
                            percentageOfDepth = contig.links[end][n].depth/totalNeighborCoverage
                            newSegment = Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i*percentageOfDepth for i in contig.depths])
                            segments.append(newSegment)
                            add_link(segments[-1], end, contig.links[end][n] , contig.otherEndOfLinks[end][n], "0M")
                            for on, otherneighbor in enumerate(contig.links[1-end]) :
                                add_link(newSegment, 1-end, otherneighbor, contig.otherEndOfLinks[1-end][on], "0M")

        for s in sorted(toDelete, reverse=True):
            continueDuplication = True
            segments[s].cut_all_links()
            del segments[s]

        segments = merge_adjacent_contigs(segments)

    return segments

#input: segments
#output : segments trimmed to reduce the number of overlaps
def trim_overlaps(segments):

    something_changes = True
    while something_changes :
        something_changes = False

        for s in segments :
            min_overlap_left = 0
            if len(s.CIGARs[0]) != 0 :
                #check that there are nothing else than Ms in the CIGAR before trimming
                only_M = True
                for cig in s.CIGARs[0] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    min_overlap_left = min([int(i.strip('M')) for i in s.CIGARs[0]])
            min_overlap_right = 0
            if len(s.CIGARs[1]) != 0 :
                only_M = True
                for cig in s.CIGARs[1] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    min_overlap_right = min([int(i.strip('M')) for i in s.CIGARs[1]])

            max_overlap_left = 0
            if len(s.CIGARs[0]) != 0 :
                only_M = True
                for cig in s.CIGARs[0] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    max_overlap_left = max([int(i.strip('M')) for i in s.CIGARs[0]]+[0])
            max_overlap_right = 0
            if len(s.CIGARs[1]) != 0 :
                only_M = True
                for cig in s.CIGARs[1] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    max_overlap_right = max([int(i.strip('M')) for i in s.CIGARs[1]]+[0])

            already_trimmed = s.get_trim()

            trim_left = min(min_overlap_left, s.length-max_overlap_right)
            trim_right = min(min_overlap_right, s.length-max_overlap_left)


            s.set_trim([already_trimmed[0]+trim_left, already_trimmed[1]+trim_right])

            if (trim_left > 0 or trim_right > 0) :
                something_changes = True

            leftCIGARs = [str(int(i.strip('M'))-trim_left)+'M' for i in s.CIGARs[0]]
            for l in range(len(s.links[0])) :
                newcigar = leftCIGARs[l]
                s.links[0][l].set_CIGAR(s.otherEndOfLinks[0][l], s, 0, newcigar)
                s.set_CIGAR(0, s.links[0][l], s.otherEndOfLinks[0][l], newcigar)

            rightCIGARs = [str(int(i.strip('M'))-trim_right)+'M'for i in s.CIGARs[1]]
            for l in range(len(s.links[1])) :
                newcigar = rightCIGARs[l]
                s.links[1][l].set_CIGAR(s.otherEndOfLinks[1][l], s, 1, newcigar)
                s.set_CIGAR(1, s.links[1][l], s.otherEndOfLinks[1][l], newcigar)

            # if (trim_left > 0 or trim_right > 0) and s.names == ['edge_165@0@0', 'edge_167@0@0', 'edge_169@0@0'] :
            #     print("qfdjqsdjlm after changeds : ", s.CIGARs[0], " ", s.CIGARs[1])

            

    return segments

#input : contig ID and fasta file
#output : sequence
def get_contig_FASTA(fastaFile, contig, firstline=0):

    with open(fastaFile) as f:

        lookAtNextLine = False
        linenumber = 0
        for line in f:
            if linenumber >= firstline:
                if lookAtNextLine:
                    return line
                target = ">" + str(contig)
                if target in line:
                    lookAtNextLine = True

            linenumber += 1
    return "In get_contig : the contig you are seeking is not in the fasta file"

#input : contig ID, gfa file and contigOffset, the position of the contig in the GFA file
#output : sequence, and if it is present, the sequencing depth of the contig and the rest of the optional tags that could be present in the input gfa
def get_contig_GFA(gfaFile, contig, contigOffset):
       
    with open(gfaFile) as f:

        f.seek(contigOffset)
        line = f.readline()         
        sline = line.strip('\n').split('\t')
            
        if len(sline) >= 3 and sline[0] == 'S' and (contig in sline[1]):
            extra_tags = ''
            depth = ''
            tags = []
            for i in range(3, len(sline)) :
                tags += sline[i].split()
            for f in tags :
                if 'dp' in f or 'DP' in f or 'km' in f or 'RC' in f:
                    depth = f
                else :
                    extra_tags += f + ' '
            extra_tags = extra_tags.strip(" ")
            
            sequence = sline[2]
            if sequence == "*" :
                sequence = ""
            return sequence, depth, extra_tags

        else :
            print('ERROR : Problem in the offset file, not pointing to the right lines')

    return "In get_contig : the contig you are seeking is not in the gfa file"


# Input :
#   offset file is for speeding up exportation
#   merge_adjacent_contig is to produce a GFA with contigs merged
def export_to_GFA(listOfSegments, gfaFile="", exportFile="results/newAssembly.gfa", offsetsFile = "", merge_adjacent_contigs = False, rename_contigs = False): 
    
    #compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig
    noOffsets = False
    #print('Offsets : ', offsetsFile)
    if offsetsFile == "" :
        noOffsets = True
        offsetsFile = gfaFile.strip('.gfa') + '_offsets.pickle'
        
    if gfaFile != "" and noOffsets:
        #print("coucou")
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile :
            for line in gfafile:
                sline = line.strip('\n').split('\t')
                if sline[0] == 'S' :
                    #print('In export_to_GFA : exporting ', sline[1])
                    line_offset[sline[1]] = offset #adds pair sline[1]:offset to the dict
                    
                offset += len(line)
            
        # with open(offsetsFile, 'wb') as o:
        #     pickle.dump(line_offset, o)
    
    # if gfaFile != '' :
    #     with open(offsetsFile, 'rb') as o:
    #         line_offset = pickle.load(o)
        
        #print(line_offset)
 
    #print('Line_offsets computed, launching proper writing of the new GFA')
    #Now that the preliminary work is done, start writing the new gfa file    

    # now sort the segments by length, to output at the beginning of the files the longests fragments
    listOfSegments.sort(key = lambda x : x.length, reverse = True)


    f = open(exportFile, "w")

    #write the sequences and the links within the supercontigs
    
    if merge_adjacent_contigs == False :
        for s, segment in enumerate(listOfSegments):
            sequences = segment.get_sequences()
            if len(sequences) != len(segment.names) :
                sequences = [None for i in range(len(segment.names))]
    
            for c, contig in enumerate(segment.names):
                
                f.write("S\t" + contig + "-" + str(segment.copiesnumber[c]) + "\t")
                if gfaFile != "":
                    sequence, depth, extra_tags = get_contig_GFA(gfaFile, contig, line_offset[contig])
                    if sequences[c] != None :
                        sequence = sequences[c]
                    #print("Here is the depth I got : ", depth)
                    if extra_tags != "" :
                        extra_tags = "\t"+extra_tags
                    if depth == '':
                        f.write(sequence + extra_tags +"\n")
                    else :
                        newdepth = str(float(depth.split(':')[-1])/copies[contig])
                        f.write(sequence + '\t' + ":".join(depth.split(':')[:-1]) + ":" + newdepth +  extra_tags + '\n')
                else:
                    f.write("*\n")
    
                if c > 0:
                    
                    f.write("L\t"+ segment.names[c-1]+ "-"+ str(segment.copiesnumber[c-1]))
                    
                    if segment.orientations[c-1] == 1 :                    
                        f.write("\t+\t")
    
                    elif segment.orientations[c-1] == 0:
                        f.write("\t-\t")
                 
                    f.write(contig + "-"+ str(segment.copiesnumber[c]))
                    
                    if segment.orientations[c] == 1 :                    
                        f.write("\t+\t")
    
                    elif segment.orientations[c] == 0:
                        f.write("\t-\t")
                        
                    f.write(segment.insideCIGARs[c-1]+'\n')
    
        print('Done exporting sequences, just a little more time...')
        #then write in the gfa file the links between the ends of supercontigs
    
        for s, segment in enumerate(listOfSegments):
                
            for endOfSegment in range(2):
                for l, neighbor in enumerate(segment.links[endOfSegment]):
                    
                    if segment.ID <= neighbor.ID : #that is to ensure each link is written only once
                    
                            
                        endOfNeighbor = segment.otherEndOfLinks[endOfSegment][l]
                        orientation1, orientation2 = '-', '-'
                        
                        if segment.orientations[-endOfSegment] == endOfSegment :
                            orientation1 = '+'
                            
                        if neighbor.orientations[-endOfNeighbor] != endOfNeighbor :
                            orientation2 = '+'
                            
            
                        f.write("L\t"+segment.names[-endOfSegment] +"-"+ str(segment.copiesnumber[-endOfSegment]) + '\t' \
                                + orientation1 + '\t' +\
                                    neighbor.names[-endOfNeighbor] +"-"+ str(neighbor.copiesnumber[-endOfNeighbor])+'\t'\
                                +orientation2+'\t'+segment.CIGARs[endOfSegment][l]+'\n')

    # in the case the user prefers having merged contigs as an output
    else : #if merge_adjacent_contigs == True
        #open a file recording which contigs correspond to which supercontigs (with lines such as supercontig_1 contig_A_contig_B_contig_C). Also store that information in a dictionary
        if rename_contigs :
            splitName = exportFile.split('/')[:-1]
            if len(splitName) > 0 :
                fcontigs = open('/'.join(splitName)+'/supercontigs.txt', 'w') 
            else :
                fcontigs = open('supercontigs.txt', 'w') 

            supercontigs = {}
            for s, segment in enumerate(listOfSegments):
                supercontigs[segment.full_name()] = "supercontig_"+ str(s)
            
        for s, segment in enumerate(listOfSegments):

            sequences = segment.get_sequences()
            if len(sequences) != len(segment.names) :
                sequences = [None for i in range(len(segment.names))]
            
            if rename_contigs :
                f.write("S\t" + "supercontig_"+ str(s) + "\t") #the name of the contigs are supercontig_i
                fcontigs.write("supercontig_"+ str(s) + "\t"+segment.full_name()+"\n")
            else :
                f.write("S\t" + segment.full_name() + "\t") #the name of the contigs are supercontig_i            
            
            fullDepth = 0
            
            if gfaFile != "":
                
                sequence = ""
                all_sequences = []
                for c, contig in enumerate(segment.names) :
                    seq, depth, extra_tags = get_contig_GFA(gfaFile, contig, line_offset[contig])
                    if sequences[c] != None :
                        seq = sequences[c]
                    if segment.orientations[c] == 0 :
                        seq = seq[::-1]
                        complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
                        seq = ''.join([complement_dict[base] for base in seq])
                    if c > 0 :
                        CIGARlength = np.sum([int(i) for i in re.findall(r'\d+', segment.insideCIGARs[c-1])])
                        # print("in input_output.py dqkdmjc ", CIGARlength, " ", segment.insideCIGARs, " ", segment.names)
                        seq = seq[CIGARlength:]
                    if depth != '' :
                        fullDepth += float(depth.split(':')[-1]) * len(seq)
                    
                    all_sequences += [seq]
                
                sequence = "".join(all_sequences)

                trimmed_ends = segment.get_trim()
                sequence = sequence[trimmed_ends[0]:len(sequence)-trimmed_ends[0]-trimmed_ends[1]]
                    
                if fullDepth == 0 or len(sequence) == 0:
                    f.write(sequence + "\n")
                else :
                    newdepth = str(fullDepth/len(sequence))
                    f.write(sequence + '\tDP:f:'+newdepth + '\n')
                
            else:
                f.write("*\n")
                
            for endOfSegment in range(2) :
                for n, neighbor in enumerate(segment.links[endOfSegment]):
                    if segment.ID < neighbor.ID : #to write each link just one
                        orientation1, orientation2 = '+', '+'
                        if endOfSegment == 0 :
                            orientation1 = '-'
                        if segment.otherEndOfLinks[endOfSegment][n] == 1 :
                            orientation2 = '-'
                        
                        if not rename_contigs :
                            f.write("L\t"+segment.full_name()+'\t'+orientation1+'\t'+neighbor.full_name()+\
                                    '\t'+orientation2+'\t'+ segment.CIGARs[endOfSegment][n]+'\n')
                        else :
                            f.write("L\t"+supercontigs[segment.full_name()]+'\t'+orientation1+'\t'+supercontigs[neighbor.full_name()]+\
                                    '\t'+orientation2+'\t'+ segment.CIGARs[endOfSegment][n]+'\n')
                                
def export_to_fasta(listOfSegments, gfaFile, exportFile="results/newAssembly.fasta", rename_contigs = False): 
    
    #compute the offsetfile : it will be useful for speeding up exportation. It will enable get_contig not to have to look through the whoooooole file each time to find one contig

    t = 0
    noOffsets = True
    offsetsFile = gfaFile.strip('.gfa') + '_offsets.pickle'
        
    if gfaFile != "" and noOffsets:
        #print("coucou")
        line_offset = {}
        offset = 0
        with open(gfaFile) as gfafile :
            for line in gfafile:
                sline = line.strip('\n').split('\t')
                if sline[0] == 'S' :
                    #print('In export_to_GFA : exporting ', sline[1])
                    line_offset[sline[1]] = offset #adds pair sline[1]:offset to the dict
                    
                offset += len(line)
            
 
    #print('Line_offsets computed, launching writing of the fasta')
    #Now that the preliminary work is done, start writing the new fasta file    

    f = open(exportFile, "w")
    
    #compute the copiesnumber
    copies = compute_copiesNumber(listOfSegments)


    # now sort the segments by length, to output at the beginning of the files the longests fragments
    listOfSegments.sort(key = lambda x : x.length, reverse = True)
    
    
    #Finally, write the sequences
    for s, segment in enumerate(listOfSegments):

        sequences = segment.get_sequences()
        if len(sequences) != len(segment.names) :
            sequences = [None for i in range(len(segment.names))]
        
        if rename_contigs :
            f.write(">supercontig_" + str(s+1) + "\n")
        else :
            f.write(">"+segment.full_name() + "\n")
        
        fullDepth = 0
        

        
        sequence = ''
        for c, contig in enumerate(segment.names) :
            seq, depth, extra_contigs = get_contig_GFA(gfaFile, contig, line_offset[contig])
            if sequences[c] != None :
                seq = sequences[c]
            if segment.orientations[c] == 0 :
                seq = seq[::-1]
                complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
                seq = ''.join([complement_dict[base] for base in seq])
            if c > 0 :
                CIGARlength = np.sum([int(i) for i in re.findall(r'\d+', segment.insideCIGARs[c-1])])
                
                seq = seq[CIGARlength:]
            if depth != '' :
                fullDepth += ( float(depth.split(':')[-1])/copies[contig] ) * len(seq)
            
            sequence += seq

        trimmed_ends = segment.get_trim()
        sequence = sequence[trimmed_ends[0]:len(sequence)-trimmed_ends[0]-trimmed_ends[1]]
            
        f.write(sequence + "\n")


# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
# Also returns the list of the contig's names
def load_gfa(file):

    gfa_read = open(file, "r")

    segments = []
    
    index = 0
    names = {} # names is a dictionary that associates the name of each contig in the gfa with an index (which will correspond later to the one in interactionMatrix and copiesnumber)
    
    for line in gfa_read:
        if line[0] == "S":

            l = line.strip('\n').split("\t")
            if line[4] == '\t' :
                l = l[:4] + [''] + l[4:]
            cov = 0
            
            for element in l :
                if 'dp' in element[:2] or 'DP' in element[:2] or 'km' in element[:2] or 'rc' in element[:2] :
                    try :
                       cov = float(element.split(":")[-1])
                    except:
                        pass
                        
                elif 'RC' in element[:2] or 'KC' in element[:2] :
                    try :
                       cov = float(element.split(":")[-1])/len(l[2])
                    except:
                        pass
            
            s = Segment([l[1]], [1], [len(l[2])], readCoverage = [cov])
            segments.append(s)
            names[s.names[0]] = index #now this contig (identified by its name) is attached to index
            index += 1            

    gfa_read = open(file, "r")
        
    cov = 1
    for line in gfa_read:
        if line[0] == "L":

            l = line.strip('\n').split("\t")
            
            segments[names[l[1]]].add_link_from_GFA(line, names, segments, 0)
            segments[names[l[3]]].add_link_from_GFA(line, names, segments, 1)

    gfa_read.close()
    
    delete_links_present_twice(segments)

    return segments, names


def main():

    gfa_in = sys.argv[1]
    gfa_out = sys.argv[2]

    segments, names = load_gfa(gfa_in)
    merge_adjacent_contigs(segments)
    # trim_overlaps(segments)
    export_to_GFA(segments, gfaFile = gfa_in, exportFile = gfa_out, merge_adjacent_contigs=True, rename_contigs=False)

if __name__ == "__main__" :
    main()