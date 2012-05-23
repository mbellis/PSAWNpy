#!/usr/bin/python
'''
Created on 26 mai 2009

@author: bellis
'''
# Copyright (c) 2009 the authors listed at the following URL, and/or
# the authors of referenced articles or incorporated external code:
# http://en.literateprograms.org/Boyer-Moore_string_search_algorithm_(Python)?action=history&offset=20081102100510
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Retrieved from: http://en.literateprograms.org/Boyer-Moore_string_search_algorithm_(Python)?oldid=15398


def match(pattern, text,mismatchNb,singleHit=0):
    matches = []
    mismatches = []    
    m = len(text)
    n = len(pattern)

    rightMostIndexes = preprocessForBadCharacterShift(pattern)
    
    alignedAt = 0
    found=0
    while alignedAt + (n - 1) < m:
#        print '\naligneAt = %u'%alignedAt
        currMismatchNb=0
        for indexInPattern in xrange(n-1, -1, -1):
            indexInText = alignedAt + indexInPattern
            x = text[indexInText]
            y = pattern[indexInPattern]

            if indexInText >= m:
                break
            print x,' ',y
                        
            if x != y:
                currMismatchNb+=1
                jump=1
                if currMismatchNb>mismatchNb:
                    jump=0
                if jump==0:                        
                                            
    #                print 'indexInPattern = %u'%indexInPattern
    #                print pattern
    #                print text [indexInText-indexInPattern:indexInText-indexInPattern+n]
    
                    r = rightMostIndexes.get(x)                
                    print 'p=%s, x=%s'%(y,x)
                    print 'rightMostIndex ',r
    
                    if not r:
                        alignedAt = indexInText + 1
    
                    else:
                        shift = indexInText - (alignedAt + r)
                        print 'shift = %u'%shift
                        print 'plus =  %u'%(shift > 0 and shift or 1)
                        alignedAt += (shift > 0 and shift or 1)
    
                    break
                elif indexInPattern == 0:
                    matches.append(alignedAt)
                    mismatches.append(currMismatchNb)
                    alignedAt += 1
                    print mismatches
                    print matches
            elif indexInPattern == 0:
                matches.append(alignedAt)
                mismatches.append(currMismatchNb)
                found=1                                
                alignedAt += 1
                print mismatches
                print matches
        if found==1 and singleHit==1:
            break            

    return [mismatches,matches]

def preprocessForBadCharacterShift(pattern):
    map = { }
    charNb=0
    for i in xrange(len(pattern)-2, -1, -1):
        c = pattern[i]
        if c not in map:
            map[c] = i
            charNb+=1
#        if charNb==4:
#            break    

    return map


def match_with_miss(pattern='',seq='',mismatchNb=0):
    matches = []
    mismatches = []
    seqLength = len(seq)
    patternLength = len(pattern)

    rightmostIndex = preprocess_for_shift(pattern)    
    seqStartIndex = 0    
    while seqStartIndex + patternLength <= seqLength:
#        print '\nseqStart = %u'%seqStartIndex
        currMismatchNb=0
            
        if patternLength%2==0:
            patternIndexes=range(patternLength-2, -1, -2)
        else:
            patternIndexes=range(patternLength-2, 0, -2)            
            patternIndexes.append(0)                      
        for patternIndex in patternIndexes:
            seqIndex = seqStartIndex + patternIndex
            if patternIndex>0:
                x = seq[seqIndex:seqIndex+2]
                y = pattern[patternIndex:patternIndex+2]               
            else:
                if patternLength%2==1:
                    x = seq[seqIndex:seqIndex+1]
                    y = pattern[patternIndex:patternIndex+1]
                else:                    
                    x = seq[seqIndex:seqIndex+2]
                    y = pattern[patternIndex:patternIndex+2]
            if x != y:                            
                #count the mismatches only if it is in first position
#                print 'x = %s y = %s patternIndex = %u'%(x,y,patternIndex)
                try:
                    if x[0]!=y[0] and x[1]!=y[1]:   
                        currMismatchNb+=2
                    else:            
                        currMismatchNb+=1
                except:                       
                    currMismatchNb+=1
                                                                                                                                                                    
                dropOut=0                
                if currMismatchNb>mismatchNb:
                    dropOut=1
#                else:
#                    print 'p=%s, x=%s'%(y,x)
                if dropOut==1:                                                                    
#                    print 'patternIndex = %u'%patternIndex
#                    print pattern
#                    print seq [seqIndex-patternIndex:seqIndex-patternIndex+patternLength]
    
                    r = rightmostIndex.get(x)                
#                    print 'p=%s, x=%s'%(y,x)
#                    print 'rightMostIndex ',r
    
                    if not r:  
                        if seqIndex>seqStartIndex:                       
                            seqStartIndex=seqIndex
                        else:
                            seqStartIndex+=1
                            
#                        print 'seqIndex = %u'%seqIndex
#                        print 'seqStartIndex = %u'%seqStartIndex                           
    
                    else:
                        if len(x)==2: 
                            shift = seqIndex - (seqStartIndex + r)-1
                        else:
                            shift = seqIndex - (seqStartIndex + r)
#                        print 'seqStartIndex = %u'%seqStartIndex
#                        print 'seqIndex = %u'%seqIndex
#                        print 'shift = %u'%shift    
                        shift=shift > 0 and shift or 1
#                        print 'shift = %u'%shift                        
                        seqStartIndex += shift
    
                    break
                
                elif patternIndex == 0:
                    matches.append(seqStartIndex)
                    mismatches.append(currMismatchNb)
                    seqStartIndex += 1
            elif patternIndex == 0:
                matches.append(seqStartIndex)
                mismatches.append(currMismatchNb)                
                seqStartIndex += 1


    return [mismatches,matches]

def preprocess_for_shift(pattern):
    map = { }                    
    for i in range(len(pattern)-3, -1, -1):
        c = pattern[i:i+2]                
        if c not in map:
            map[c] = i
    if len(pattern)%2==1:
        charNb=0
        for i in range(len(pattern)-2, -1, -1):
            c = pattern[i]                   
            if c not in map:
                charNb=0
                map[c] = i
            if charNb==4:
                break 
#    print pattern                              
#    for key in map.keys():
#        print '%s => %u '%(key,map[key])    
    return map


if __name__ == "__main__":
    
#    [mismatchNb,position]=match('ATTG','GATTGCCACTGCCCATTGCT',mismatchNb=1,singleHit=1)#
    [mismatchNb,position]=match_with_miss(pattern='CAAACAGCCCAGAAGACCAGCAGCA',seq='AGAAGACCTAGAGGGGAAGGAGCAGGACCAAACAGCCCATCAGACCAGCAGCATTGGCAACATTC',mismatchNb=1)#
    print mismatchNb
    print position
#    CAAACAGCCCAGAAGACCAGCAGCA
#    TAAACAGCCCAGAAGACCAGAAGCA
