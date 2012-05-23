#!/usr/bin/python
"""**Functions for testing position membership and processing exons.**"""

from numpy import array

def owningStructure(position,strand,starts,ends,strands,overlapFlag):
    """
    **Find the structure(s) wich contains the tested position.**
    ::
    
        Arguments:
            position ----- coordinate of the tested position
            strand ------- strandedness (-1 or -1) of the tested position
            starts ------- array of starting coordinates of the scanned structures
            ends --------- array of ending coordinates of the scanned structures
            strands ------ array of strandedness of the scanned structures
            overlapFlag -- allows to retrieve several overlapping structures
            
        Return: 
            indexes -- array of indexe(s) of structure(s) containing the tested position   
         
    """         
  
    greater = starts<=position
    lesser = ends>=position
    strandened = strands==strand
    # Get an array of boolean which intersect the three previous boolean arrays
    isOwned = greater & lesser & strandened            
    indexes=[]            
    for i in range(len(isOwned)):
        # test if the tested position is in structure of index i                        
        if isOwned[i]==True:                                            
            indexes.append(i)
            if overlapFlag==0:
                # the position is unique and/or it is not necessary to test overlapping structures
                break                    
    return indexes
    
def isInStructure(start,end,strand,starts,ends,strands,overlapFlag):
    """
    **Find the structure(s) wich contains the tested structure.**
    ::
    
        Arguments:
            start    ----- starting position of the tested structure
            end      ----- ending position of the tested structure
            strand ------- strandedness (-1 or -1) of the tested structure
            starts ------- array of starting coordinates of the scanned structures
            ends --------- array of ending coordinates of the scanned structures
            strands ------ array of strandedness of the scanned structures
            overlapFlag -- allows to retrieve several overlapping structures
        
        Return: 
            indexes -- array of indexe(s) of structure(s) containing the tested structure
            
    """
    greater = starts<=start
    lesser = ends>=end
    strandened = strands==strand
    # Get an array of boolean which intersect the three previous boolean arrays
    isOwned = greater & lesser & strandened            
    indexes=[]                    
    for i in range(len(isOwned)):
        # test if the tested structure is in structure of index i                            
        if isOwned[i]==True:                                                            
            indexes.append(i)
            if overlapFlag==0:
                # the position is unique and/or it is not necessary to test overlapping structures
                break                    
    return indexes  
  
def isOverlapStructure(start,end,strand,starts,ends,strands,overlapFlag):    
    """
    **Find the structure(s) wich overlap the tested structure.**
    ::
    
        Arguments:
            start -------- starting position of the tested structure
            end ---------- ending position of the tested structure
            strand ------- strandedness (-1 or -1) of the tested structure
            starts ------- array of starting coordinates of the scanned structures
            ends --------- array of ending coordinates of the scanned structures
            strands ------ array of strandedness of the scanned structures
            overlapFlag -- allows to retrieve several overlapping structures
        
        Return: 
            indexes -- array of indexe(s) of structure(s) overlaping the tested structure
            
    """            
                        
    smallerStarts = starts<=start
    greaterStarts = ends>=start        
    greaterEnds = ends>=end
    smallerEnds = starts<=end        
    strandened = strands==strand
    # Get an array of boolean which intersect the five previous boolean arrays
    overlaps = ((smallerStarts & greaterStarts)^(smallerEnds & greaterEnds)) & strandened            
    indexes=[]                    
    for i in range(len(overlaps)):
        # test if the tested structure is overlaped by structure of index i                            
        if overlaps[i]==True:                            
            # the position is unique and it is allowed to test overlapping structures                
            indexes.append(i)
            if overlapFlag==0:
                break                    
    return indexes  
      

def process_exons(exonIDs=[], exonStarts=[],exonEnds=[],exonStrands=[],transcriptIDs=[]):
    """
    **Order exons and eliminate doublons.**
    ::
     
        Keywords arguments:
            exonIDs -------- list of exon IDs
            exonStarts ----- list of exon starting positions
            exonEnds ------- list of exon ending positions
            exonStrands ---- list of exon strandedness
            transcriptIDs -- list of transcript IDs
            
        Return ordered list without doublons:
            exonIDList --------- list of exon IDs
            exonStarts --------- array of exon starting positions
            exonEnds ----------- array of exon ending positions
            exonStrands -------- array of exon strandedness
            transcriptsByExon -- list of sets of transcripts
            transcriptIDs ------ list of transcript IDs
            transcriptStarts --- array of transcript starts (smallest targeting exon start) 
            transcriptEnds ----- array of transcript ends (largest targeting exon end)
            exonIndexes -------- array of indexes indicating the rank of exon start positions relative
                                 to the 5' gene start position
            exonGroups --------- indicate for each exon, the group it belongs to
                                 (a chain of overlaping exons is a group)
            intronStarts ------- start positions of introns defined as inter-group intevals
            intronEnds --------- end positions of introns defined as inter-group intevals
             
    """            

    if transcriptIDs==[]:
        transcriptIDs=array([1]*len(exonStarts))        
    if exonIDs==[]:
        exonIDs=array([1]*len(exonStarts))        
    # order the exon limit            
    sortIndex=exonEnds.argsort(kind='mergesort')
    exonIDs=exonIDs[sortIndex]
    exonStarts=exonStarts[sortIndex]
    exonEnds=exonEnds[sortIndex]
    exonStrands=exonStrands[sortIndex]     
    transcriptIDs=transcriptIDs[sortIndex]
      
    sortIndex=exonStarts.argsort(kind='mergesort')
    exonIDs=exonIDs[sortIndex]
    exonStarts=exonStarts[sortIndex]
    exonEnds=exonEnds[sortIndex]
    exonStrands=exonStrands[sortIndex]     
    transcriptIDs=transcriptIDs[sortIndex]  
         
    # eliminate exon doublons    
    Start=exonStarts[0]
    End=exonEnds[0]
    #list of exon IDs 
    exonIDList=[exonIDs[0]]
    #list of exon starts    
    exonStartList=[Start]
    #list of exon ends
    exonEndList=[End]
    #list of exon strands
    exonStrandList=[exonStrands[0]]  
    #list of sets of IDs of the transcript(s) targeted by each exon  
    transcriptsByExon=[set([transcriptIDs[0]])]
    #set of IDs of all targeted transcripts    
    transcriptIDSet=set([transcriptIDs[0]])
    #dictionnary of transcripts with the smallest starting value of targeting exons 
    transcriptStartsDict={}
    #dictionnary of transcripts with the largest  of targeting exonsnding value
    transcriptEndsDict={}
    for transcript in transcriptsByExon[0]:        
        transcriptStartsDict[transcript]=exonStarts[0]
        transcriptEndsDict[transcript]=exonEnds[0]
    
    exonPos=0    
    for i in range(1,len(exonStarts)):
        if exonStarts[i]!=Start or exonEnds[i]!=End:            
            exonPos+=1
            exonIDList.append(exonIDs[i])
            exonStartList.append(exonStarts[i])
            exonEndList.append(exonEnds[i])
            exonStrandList.append(exonStrands[0])
            transcriptsByExon.append(set([transcriptIDs[i]]))
            transcriptIDSet=transcriptIDSet.union(set([transcriptIDs[i]]))
            Start=exonStarts[i]
            End=exonEnds[i]        
            
        else:
            transcriptsByExon[exonPos]=transcriptsByExon[exonPos].union(set([transcriptIDs[i]]))
            transcriptIDSet=transcriptIDSet.union(set([transcriptIDs[i]]))
        for transcript in [transcriptIDs[i]]:                
            try:                    
                transcriptStartsDict[transcript]=min(transcriptStartsDict[transcript],exonStarts[i])
                transcriptEndsDict[transcript]=max(transcriptEndsDict[transcript],exonEnds[i])
            except:                                                
                transcriptStartsDict[transcript]=exonStarts[i]
                transcriptEndsDict[transcript]=exonEnds[i]   
                
    exonStarts=array(exonStartList)   
    exonEnds=array(exonEndList)
    exonStrands=array(exonStrandList)    
    transcriptIDs=list(transcriptIDSet) 
        
    # define start and end of transcripts by taking the smallest and greatest
    # start and end of targeting exons                                    
    transcriptStarts=[]
    transcriptEnds=[]
    for transcript in transcriptIDs:
        transcriptStarts.append(transcriptStartsDict[transcript])
        transcriptEnds.append(transcriptEndsDict[transcript])
    transcriptStarts=array(transcriptStarts)    
    transcriptEnds=array(transcriptEnds)
                    
    # groups overlaping exons                               
    exonGroups=[-1]*len(exonStarts)
    group=-1                        
    for i in range(len(exonStarts)):  
        if  exonGroups[i]==-1:
            #construct a new group of exons
            group+=1
            exonGroups[i]=group
            endGroup=exonEnds[i]
            for j in range(i+1,len(exonStarts)):
                if exonStarts[j]<=endGroup:
                    exonGroups[j]=group
                    endGroup=max(endGroup,exonEnds[j])
                                            
    exonIndexes=range(len(exonStarts))                  
    if exonStrands[0]==-1:
        #reverse index if strand == -1
        exonIndexes=exonIndexes[::-1]
        for i in range(len(exonGroups)):
            exonGroups[i]=group-exonGroups[i]                                              
                
    #find limits of exon groups        
    exonGroupStarts=[exonStarts[0]]
    exonGroupEnds=[exonEnds[0]]                                                     
    group=exonGroups[0]                      
    for i in range(1,len(exonStarts)):
        if  exonGroups[i]!=group:                                                                                                                         
            exonGroupEnds.append(exonEnds[i])
            exonGroupStarts.append(exonStarts[i])
            group= exonGroups[i]                               
    exonGroupEnds.append(exonEnds[i])
    exonGroupStarts=array(exonGroupStarts)
    exonGroupEnds=array(exonGroupEnds)
    
    #find limits of introns defined as intervals between exon groups
    intronStarts=[]
    intronEnds=[]    
    if len(exonGroupStarts)>1:
        for i in range(len(exonGroupStarts)-1):    
            intronStarts.append(exonGroupEnds[i]+1)    
            intronEnds.append(exonGroupStarts[i+1]-1)
        intronStarts=array(intronStarts)
        intronEnds=array(intronEnds)     
    return [exonIDList, exonStarts, exonEnds, exonStrands, transcriptsByExon, transcriptIDs, transcriptStarts, transcriptEnds, exonIndexes, exonGroups, intronStarts, intronEnds]               

