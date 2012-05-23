#!/usr/bin/python
"""**Functions for displaying results.**"""
import bsddb, cPickle, os, sys
import setenviron
import getopt
import bintools
from numpy import where, array
setenviron.setenviron()


def test_chip(chipIDs='',chipRange='',testFile='',tutoFlag=0):
    path=os.environ['COMMONDB']+'/chip.bkdb'
    chipDB=bsddb.btopen(path)
    dbKeys=chipDB.keys()
    for i in range(2):
        write=0
        if i==0:
            if chipIDs!='':
                write=1
        else:
            if chipRange!='':
                write=1
                if chipRange==0:
                    chipIDs=dbKeys
                else:
                    chipIDs=list()                      
                    for chipRank in chipRange:
                        chipIDs.append(dbKeys[chipRank])                    
        if write:
            for chipID in chipIDs:    
                chip=cPickle.loads(chipDB[chipID])
                testFile.write('\t%s chip %s (rank %u)\n'%('*'*5,chipID,dbKeys.index(chipID)))
                testFile.write('\t\tmyName: %s\n'%chip.myName)
                testFile.write('\t\tname: %s\n'%chip.name)
                testFile.write('\t\tshortName: %s\n'%chip.shortName)
                testFile.write('\t\tmySpeciesName: %s\n'%chip.mySpeciesName)
                testFile.write('\t\tprobesetNb: %u\n'%chip.probesetNb)
                testFile.write('\t\tprobeNb: %u\n'%chip.probeNb)
                testFile.write('\t\tcompName: %s\n'%chip.compName)
                testFile.write('\t\tens47Name: %s\n'%chip.ens47Name)
                testFile.write('\t\tens48Name: %s\n'%chip.ens48Name)
                testFile.write('\t\tgeoName: %s\n'%chip.geoName)
                testFile.write('\t\n')            
    if tutoFlag:
        testFile.write('\t%s\n'%('.'*5))                                         
    chipDB.close()        
        
def test_species(speciesIDs='',speciesRange='',testFile='',tutoFlag=0):
    path=os.environ['COMMONDB']+"/species.bkdb"
    speciesDB=bsddb.btopen(path)    
    dbKeys=speciesDB.keys()
    for i in range(2):
        write=0
        if i==0:
            if speciesIDs!='':
                write=1
        else:
            if speciesRange!='':
                write=1
                if speciesRange==0:
                    speciesIDs=dbKeys
                else:
                    speciesIDs=list()                      
                    for speciesRank in speciesRange:
                        speciesIDs.append(dbKeys[speciesRank])                    
        if write:
            for speciesID in speciesIDs:    
                species=cPickle.loads(speciesDB[speciesID])
                testFile.write("\t%s species %s (rank %u)\n"%("*"*5,speciesID,dbKeys.index(speciesID)))
                testFile.write("\t\tmyName: %s\n"%species.myName)
                testFile.write("\t\tofficialName: %s\n"%species.officialName)
                testFile.write('\t\n')    
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))                        
    speciesDB.close()
    
def test_probeset(action,species,chipName,testFile='',probeSetIDs='',probeSetRange='',tutoFlag=0):
                    
    path=os.environ['PYDATA']+"/"+species+"/"+chipName+"_probeset.bkdb"   
    probeSetDB=bsddb.btopen(path)    
    dbKeys=probeSetDB.keys()
    for i in range(2):
        write=0
        if i==0:
            if probeSetIDs!='':
                write=1
        else:
            if probeSetRange!='':
                write=1
                if probeSetRange==0:
                    probeSetIDs=dbKeys
                else:
                    probeSetIDs=list()                      
                    for probeSetRank in probeSetRange:
                        probeSetIDs.append(dbKeys[probeSetRank])                    
        if write:
            for probeSetID in probeSetIDs:                                                                            
                probeset=cPickle.loads(probeSetDB[probeSetID])
                testFile.write('\t%s probeset %s (rank %u)\n\n'%('+'*5,probeSetID,dbKeys.index(probeSetID)))                                                
                testFile.write('\t\tprobesetID: %s\n'%probeset.probesetID)
                testFile.write('\t\tprobesetIndex: %s\n'%probeset.probesetIndex)
                if action==5:
                    testFile.write('\t\tensGeneIDs: %s\n'%probeset.ensGeneIDs)                                                                 
                if action >=6:            
                    if probeset.affyGeneIDs!=None:
                        testFile.write('\t\taffyGeneIDs: %s\n'%probeset.affyGeneIDs)
                    else:
                        testFile.write('\t\taffyGeneIDs: none\n')                            
                    testFile.write('\t\tprobeNb: %u\n'%probeset.probeNb)
                    testFile.write('\t\tprobeIDs: %s\n'%str(probeset.probeIDs))
                    testFile.write('\t\tensemblProbeIDs: %s\n'%str(probeset.ensProbeIDs))
                    testFile.write('\t\tprobeIndexes: %s\n'%str(probeset.probeIndexes))        
                    testFile.write('\t\tprobesetTargetLength: %u\n'%probeset.probesetTargetLength)
                    testFile.write('\t\tprobesetTargetStart: %u\n'%probeset.probesetTargetStart)        
                    testFile.write('\t\tprobesetTargetEnd: %u\n'%probeset.probesetTargetEnd)
                if action <7:
                    testFile.write('\t\n')                    
                if action>=7:   
                    testFile.write('\t\tensemblExonGeneNbs: %s\n'%probeset.ensExonGeneNbs)
                    testFile.write('\t\tensemblIntronGeneNbs: %s\n'%probeset.ensIntronGeneNbs)
                    testFile.write('\t\tensemblUpGeneNbs: %s\n'%probeset.ensUpGeneNbs)
                    testFile.write('\t\tensemblDownGeneNb: %s\n'%probeset.ensDownGeneNbs)
                    testFile.write('\t\tensemblOutProbeNbs: %s\n'%probeset.ensOutProbeNbs)
                    testFile.write('\t\tensemblNisProbeNbs: %s\n'%probeset.ensNisProbeNbs)
                    testFile.write('\t\tourEnsemblGeneIDs: %s\n'% probeset.ourEnsGeneIDs)
                    testFile.write('\t\tourEnsemblGeneNb: %s\n'%probeset.ourEnsGeneNb)
                    testFile.write('\t\tourEnsemblProbeNb: %s\n'%probeset.ourEnsProbeNb)
                    testFile.write('\t\tEns: specific to ensemblGeneIds: %s\n'%probeset.ensGeneIDs4Ens)
                    testFile.write('\t\tEns: specific to ourEnsemblGeneIds: %s\n'%probeset.ourEnsGeneIDs4Ens)
                    testFile.write('\t\tEns: common to both: %s\n'%probeset.commonGeneIDs4Ens)
                    testFile.write('\t\taceExonGeneNbs: %s\n'%probeset.aceExonGeneNbs)
                    testFile.write('\t\taceIntronGeneNbs: %s\n'%probeset.aceIntronGeneNbs)
                    testFile.write('\t\taceUpGeneNbs: %s\n'%probeset.aceUpGeneNbs)
                    testFile.write('\t\taceDownGeneNb: %s\n'%probeset.aceDownGeneNbs)
                    testFile.write('\t\taceOutProbeNbs: %s\n'%probeset.aceOutProbeNbs)
                    testFile.write('\t\taceNisProbeNbs: %s\n'%probeset.aceNisProbeNbs)
                    testFile.write('\t\tourAceGeneIDs: %s\n'% probeset.ourAceGeneIDs)
                    testFile.write('\t\tourAceGeneNb: %s\n'%probeset.ourAceGeneNb)
                    testFile.write('\t\tourAceProbeNb: %s\n'%probeset.ourAceProbeNb)     
                    testFile.write('\t\tourAceToEnsGeneIDs: %s\n'%probeset.ourAceToEnsGeneIDs)
                    testFile.write('\t\tourAceToEnsGeneNbs: %s\n'%probeset.ourAceToEnsGeneNbs)
                    testFile.write('\t\tAce: specific to ensemblGeneIds: %s\n'%probeset.ensGeneIDs4Ace)
                    testFile.write('\t\tAce: specific to ourEnsemblGeneIds: %s\n'%probeset.ourEnsGeneIDs4Ace)
                    testFile.write('\t\tAce: common to both: %s\n'%probeset.commonGeneIDs4Ace)
                                    
                    testFile.write('\n\t\tINFORMATION ON TARGETED GENES:\n')    
                    if probeset.ensGenesByProbeNb is not None:
                        testFile.write('\n\t\tENSEMBL:\n')                
                        for key in probeset.ensGenesByProbeNb.keys():            
                            testFile.write('\t\t\n')
                            testFile.write('\t\tkey: %s\n'%key)            
                            for geneID in probeset.ensGenesByProbeNb[key].keys():
                                testFile.write('\n\t\t%s gene: %s\n'%('*'*5,geneID))                
                                testFile.write('\t\tprobeIndexes: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probeIndexes)        
                                testFile.write('\t\tprobeLocalisations: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probeLocalisations)        
                                testFile.write('\t\tfirstStructureIndexes: %s\n'%probeset.ensGenesByProbeNb[key][geneID].firstStructureIndexes)        
                                testFile.write('\t\tfirstStructureGroups: %s\n'%probeset.ensGenesByProbeNb[key][geneID].firstStructureGroups)
                                testFile.write('\t\tfirstStructureIDs: %s\n'%probeset.ensGenesByProbeNb[key][geneID].firstStructureIDs)
                                testFile.write('\t\tfirstStructureGroups: %s\n'%probeset.ensGenesByProbeNb[key][geneID].firstStructureGroups)
                                testFile.write('\t\tsndStructureIndexes: %s\n'%probeset.ensGenesByProbeNb[key][geneID].sndStructureIndexes)        
                                testFile.write('\t\tprobePositions: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probePositions)        
                                testFile.write('\t\tprobeRepetitionNbs: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probeRepetitionNbs)        
                                testFile.write('\t\tprobeStrands: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probeStrands)        
                                testFile.write('\t\tprobeMismatchNbs: %s\n'%probeset.ensGenesByProbeNb[key][geneID].probeMismatchNbs)        
                                testFile.write('\t\tin ExonProbeNb: %s\n'%probeset.ensGenesByProbeNb[key][geneID].inExonProbeNb)        
                                testFile.write('\t\tinSpliceProbeNb: %s\n'%probeset.ensGenesByProbeNb[key][geneID].inSpliceProbeNb)        
                                testFile.write('\t\tinIntronProbeNb: %s\n'%probeset.ensGenesByProbeNb[key][geneID].inIntronProbeNb)        
                                testFile.write('\t\tupProbeNb: %s\n'%probeset.ensGenesByProbeNb[key][geneID].upProbeNb)        
                                testFile.write('\t\tdownProbeNb: %s\n'%probeset.ensGenesByProbeNb[key][geneID].downProbeNb)          
                                testFile.write('\t\texonSet: %s\n'%probeset.ensGenesByProbeNb[key][geneID].exonSet)                               
                                testFile.write('\t\tgroupSet: %s\n'%probeset.ensGenesByProbeNb[key][geneID].groupSet)
                                testFile.write('\t\ttargetedTranscripts: %s\n'%probeset.ensGenesByProbeNb[key][geneID].targetedTranscripts)                
                                testFile.write('\t\tnotTargetedtranscripts: %s\n'%probeset.ensGenesByProbeNb[key][geneID].notTargetedTranscripts)
                                testFile.write('\t\tensemblGenes: %s\n'%probeset.ensGenesByProbeNb[key][geneID].ensemblGenes)
                                testFile.write('\t\t\n')
                    if probeset.aceGenesByProbeNb is not None:
                        testFile.write('\n\t\tACEVIEW:\n')
                        for key in probeset.aceGenesByProbeNb.keys():                                        
                            testFile.write('\t\t\n')
                            testFile.write('\t\tkey: %s\n'%key)            
                            for geneID in probeset.aceGenesByProbeNb[key].keys():                                
                                testFile.write('\n\t\t%s gene: %s\n'%('*'*5,geneID))               
                                testFile.write('\t\tprobeIndexes: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probeIndexes)        
                                testFile.write('\t\tprobeLocalisations: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probeLocalisations)        
                                testFile.write('\t\tfirstStructureIndexes: %s\n'%probeset.aceGenesByProbeNb[key][geneID].firstStructureIndexes)        
                                testFile.write('\t\tfirstStructureGroups: %s\n'%probeset.aceGenesByProbeNb[key][geneID].firstStructureGroups)
                                testFile.write('\t\tfirstStructureIDs: %s\n'%probeset.aceGenesByProbeNb[key][geneID].firstStructureIDs)
                                testFile.write('\t\tfirstStructureGroups: %s\n'%probeset.aceGenesByProbeNb[key][geneID].firstStructureGroups)
                                testFile.write('\t\tsndStructureIndexes: %s\n'%probeset.aceGenesByProbeNb[key][geneID].sndStructureIndexes)        
                                testFile.write('\t\tprobePositions: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probePositions)        
                                testFile.write('\t\tprobeRepetitionNbs: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probeRepetitionNbs)        
                                testFile.write('\t\tprobeStrands: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probeStrands)        
                                testFile.write('\t\tprobeMismatchNbs: %s\n'%probeset.aceGenesByProbeNb[key][geneID].probeMismatchNbs)        
                                testFile.write('\t\tin ExonProbeNb: %s\n'%probeset.aceGenesByProbeNb[key][geneID].inExonProbeNb)        
                                testFile.write('\t\tinSpliceProbeNb: %s\n'%probeset.aceGenesByProbeNb[key][geneID].inSpliceProbeNb)        
                                testFile.write('\t\tinIntronProbeNb: %s\n'%probeset.aceGenesByProbeNb[key][geneID].inIntronProbeNb)        
                                testFile.write('\t\tupProbeNb: %s\n'%probeset.aceGenesByProbeNb[key][geneID].upProbeNb)        
                                testFile.write('\t\tdownProbeNb: %s\n'%probeset.aceGenesByProbeNb[key][geneID].downProbeNb)          
                                testFile.write('\t\texonSet: %s\n'%probeset.aceGenesByProbeNb[key][geneID].exonSet)                               
                                testFile.write('\t\tgroupSet: %s\n'%probeset.aceGenesByProbeNb[key][geneID].groupSet)
                                testFile.write('\t\ttargetedTranscripts: %s\n'%probeset.aceGenesByProbeNb[key][geneID].targetedTranscripts)                
                                testFile.write('\t\tnotTargetedtranscripts: %s\n'%probeset.aceGenesByProbeNb[key][geneID].notTargetedTranscripts)
                                testFile.write('\t\tensemblGenes: %s\n'%probeset.aceGenesByProbeNb[key][geneID].ensemblGenes)
                                testFile.write('\t\t\n')
    if tutoFlag and action<6:
        testFile.write('\n\t%s\n'%('.'*5))                        
    probeSetDB.close()
    
def test_probe(action,species,chipName,dbType,testFile='',probeIDs='',tutoFlag=0):
        
    path=os.environ['PYDATA']+"/"+species+"/"+dbType+"/"+chipName+"_probe.bkdb"    
    probeDB=bsddb.btopen(path)
#    dbKeys=probeDB.keys()
#    for i in range(2):
#        write=0
#        if i==0:
#            if probeIDs!='':
#                write=1
#        else:
#            if probeRange!='':
#                write=1
#                if probeRange==0:
#                    probeIDs=dbKeys
#                else:
#                    probeIDs=list()                      
#                    for probeRank in probeRange:
#                        probeIDs.append(dbKeys[probeRank])   
#                        probeDB.                 
#        if write:            
    for probeID in probeIDs:
        probe=cPickle.loads(probeDB[probeID])                        
        testFile.write('\t%s probe %s\n'%('*'*5,probeID))                                                                
        testFile.write('\t\tprobeID: %s\n'%probe.probeID)        
        testFile.write('\t\tensemblID: %s\n'%probe.ensemblID)        
        testFile.write('\t\tensProbesetID: %s\n'%probe.ensProbesetID)                            
        if action>=6:            
            testFile.write('\t\tprobesetID: %s\n'%probe.probesetID)                                    
            testFile.write('\t\tindex: %u\n'%probe.index)
            testFile.write('\t\tsequence: %s\n'%probe.sequence)
            testFile.write('\t\ttargetPosition: %u\n'%probe.targetPosition)
        if action<=6:
            testFile.write('\t\n')    
        if action>=7:                      
            keys=probe.genes4exon.keys()
            if len(keys)>0:
                testFile.write('\t\n')
                testFile.write('\t\tgenes targeted in exons: %s\n'%keys)
                testFile.write('\t\tfirst gene targeted in exons: %s\n'%keys[0])
                testFile.write('\t\trepetitions: %s\n'%probe.genes4exon[keys[0]].repetition)
                testFile.write('\t\tIDs: %s\n'%probe.genes4exon[keys[0]].IDs)
                testFile.write('\t\tindexes: %s\n'%probe.genes4exon[keys[0]].indexes)
                testFile.write('\t\tgroups: %s\n'%probe.genes4exon[keys[0]].groups)
                testFile.write('\t\tpositions: %s\n'%probe.genes4exon[keys[0]].positions)
                testFile.write('\t\tstrands: %s\n'%probe.genes4exon[keys[0]].strands)
                testFile.write('\t\tmismatchNbs: %s\n'%probe.genes4exon[keys[0]].mismatchNbs)
             
            keys=probe.genes4intron.keys()
            if len(keys)>0:
                testFile.write('\t\n')
                testFile.write('\t\tgenes targeted in introns: %s\n'%keys)
                testFile.write('\t\tfirst gene targeted in introns: %s\n'%keys[0])
                testFile.write('\t\trepetitions: %s\n'%probe.genes4intron[keys[0]].repetition)
                testFile.write('\t\tIDs: %s\n'%probe.genes4intron[keys[0]].IDs)
                testFile.write('\t\tindexes: %s\n'%probe.genes4intron[keys[0]].indexes)
                testFile.write('\t\tgroups: %s\n'%probe.genes4intron[keys[0]].groups)
                testFile.write('\t\tpositions: %s\n'%probe.genes4intron[keys[0]].positions)
                testFile.write('\t\tstrands: %s\n'%probe.genes4intron[keys[0]].strands)
                testFile.write('\t\tmismatchNbs: %s\n'%probe.genes4intron[keys[0]].mismatchNbs)
             
            keys=probe.genes4splice.keys()
            if len(keys)>0:
                testFile.write('\t\n')
                testFile.write('\t\tgenes targeted in splices: %s\n'%keys)
                testFile.write('\t\tfirst gene targeted in splices: %s\n'%keys[0])
                testFile.write('\t\trepetitions: %s\n'%probe.genes4splice[keys[0]].repetition)
                testFile.write('\t\tIDs: %s\n'%probe.genes4splice[keys[0]].IDs)
                testFile.write('\t\tindexes: %s\n'%probe.genes4splice[keys[0]].indexes)
                testFile.write('\t\tgroups: %s\n'%probe.genes4splice[keys[0]].groups)
                testFile.write('\t\tpositions: %s\n'%probe.genes4splice[keys[0]].positions)
                testFile.write('\t\tstrands: %s\n'%probe.genes4splice[keys[0]].strands)
                testFile.write('\t\tmismatchNbs: %s\n'%probe.genes4splice[keys[0]].mismatchNbs)
             
            keys=probe.genes4up.keys()
            if len(keys)>0:
                testFile.write('\t\n')
                testFile.write('\t\tgenes targeted in up: %s\n'%keys)
                testFile.write('\t\tfirst gene targeted in up: %s\n'%keys[0])
                testFile.write('\t\trepetitions: %s\n'%probe.genes4up[keys[0]].repetition)
                testFile.write('\t\tIDs: %s\n'%probe.genes4up[keys[0]].IDs)
                testFile.write('\t\tindexes: %s\n'%probe.genes4up[keys[0]].indexes)
                testFile.write('\t\tgroups: %s\n'%probe.genes4up[keys[0]].groups)
                testFile.write('\t\tpositions: %s\n'%probe.genes4up[keys[0]].positions)
                testFile.write('\t\tstrands: %s\n'%probe.genes4up[keys[0]].strands)
                testFile.write('\t\tmismatchNbs: %s\n'%probe.genes4up[keys[0]].mismatchNbs)
             
            keys=probe.genes4down.keys()
            if len(keys)>0:
                testFile.write('\t\n')
                testFile.write('\t\tgenes targeted in down: %s\n'%keys)
                testFile.write('\t\tfirst gene targeted in down: %s\n'%keys[0])
                testFile.write('\t\trepetitions: %s\n'%probe.genes4down[keys[0]].repetition)
                testFile.write('\t\tIDs: %s\n'%probe.genes4down[keys[0]].IDs)
                testFile.write('\t\tindexes: %s\n'%probe.genes4down[keys[0]].indexes)
                testFile.write('\t\tgroups: %s\n'%probe.genes4down[keys[0]].groups)
                testFile.write('\t\tpositions: %s\n'%probe.genes4down[keys[0]].positions)
                testFile.write('\t\tstrands: %s\n'%probe.genes4down[keys[0]].strands)
                testFile.write('\t\tmismatchNbs: %s\n'%probe.genes4down[keys[0]].mismatchNbs)
            testFile.write('\t\n')
            testFile.write('\t\toutOfGeneNb: %s\n'%probe.outOfGeneNb)
            testFile.write('\t\tnotInSequence:  %s\n'%probe.notInSequence)
            testFile.write('\t\n')
            
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))        
    probeDB.close()
    

def test_transcripts_by_exon(action,species,dbType,testFile='',exonIDs=[''],exonRange='',tutoFlag=0):
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_transcripts_by_exon.bkdb"
    exonDB=bsddb.btopen(path,'r')
    dbKeys=exonDB.keys()
    for i in range(2):
        write=0
        if i==0:
            if exonIDs!='':
                write=1
        else:
            if exonRange!='':
                exonIDs=list()                      
                for exonRank in exonRange:
                    exonIDs.append(dbKeys[exonRank])
                write=1
        if write:                
            for exonID in exonIDs:                            
                transcript=cPickle.loads(exonDB[exonID])                        
                testFile.write('\n\t%s exon %s\n'%(('*'*5),exonID))
                if dbType=='ensembl':                
                    testFile.write('\t\tTranscript IDS: %s\n'%transcript.IDs)                    
                    if action>7:
                        testFile.write('\t\tTranscripts starts: %s\n'%transcript.starts)
                        testFile.write('\t\tTranscripts ends: %s\n'%transcript.ends)
                        testFile.write('\t\tTranscripts strands: %s\n'%transcript.strands)
                elif dbType=='aceview':    
                    testFile.write('\t\tTranscript IDS: %s\n'%transcript)
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))                                                                                       
    exonDB.close()
    
def test_exons_by_gene(action,species,testFile='',geneIDs='',geneRange='',tutoFlag=0,aceFlag=0):
    if aceFlag==0:
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_exons_by_gene.bkdb"
    else:
        path=os.environ['PYDATA']+"/"+species+"/aceview/"+species+"_exons_by_gene.bkdb"
    print path
    print geneIDs
    geneDB=bsddb.btopen(path,'r')  
    for i in range(2):
        write=0
        if i==0:
            if geneIDs!='':
                write=1
        else:
            if geneRange!='':
                geneIDs=list()
                dbKeys=geneDB.keys()      
                for geneRank in geneRange:
                    geneIDs.append(dbKeys[geneRank])
                write=1
        if write: 
            print geneDB.keys()[0]               
            for geneID in geneIDs:  
                testFile.write('\t%s gene %s:\n'%('*'*5,geneID))                                                    
                exons=cPickle.loads(geneDB[geneID])
                testFile.write('\t\texon.IDs = %s\n'%exons.IDs)
                testFile.write('\t\texon.indexes = %s\n'%exons.indexes)
                testFile.write('\t\texon.groups = %s\n'%exons.groups)
                testFile.write('\t\texon.exonStarts = %s\n'%exons.exonStarts)
                testFile.write('\t\texon.exonEnds = %s\n'%exons.exonEnds)
                testFile.write('\t\texon.strands = %s\n'%exons.strands)
                testFile.write('\t\texon.intronStarts = %s\n'%exons.intronStarts)
                testFile.write('\t\texon.intronEnds = %s\n'%exons.intronEnds)
                testFile.write('\t\texon.transcriptIDs = %s\n'%exons.transcriptIDs)
                testFile.write('\t\texon.transcriptStarts = %s\n'%exons.transcriptStarts)
                testFile.write('\t\texon.transcriptEnds = %s\n'%exons.transcriptEnds)            
                if tutoFlag:
                    for i in range(min(3,len(exons.transcriptsByExon))):
                        testFile.write('\t\texon.transcriptList of exon %u = %s\n'%(i,exons.transcriptsByExon[i]))
                    testFile.write('\t\t%s\n'%('.'*5))
                else:
                    for i in range(len(exons.transcriptsByExon)):
                        testFile.write('\t\texon.transcriptList of exon %u = %s\n'%(i,exons.transcriptsByExon[i]))
                testFile.write('\t\n')            
    geneDB.close()
        
def test_genes_by_region(action,species,testFile='',geneIDs='',region='',geneRange='',tutoFlag=0):    
    
    if geneIDs!='':        
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_region_by_gene.bkdb"
        geneDB=bsddb.btopen(path,'r')
        dbKeys=geneDB.keys()
        print dbKeys[0]
        for geneID in geneIDs:
            print geneID            
            testFile.write('\t%s gene %s (rank %u)\n'%('*'*5,geneID,dbKeys.index(geneID)))                           
            testFile.write('\t\tregion: %s\n'%geneDB[geneID])
        if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))                    
    if region!='' and geneRange!='':                 
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_genes_by_region.bkdb"
        regionDB=bsddb.btopen(path,'r')
        dbKeys=regionDB.keys()    
        testFile.write('\t\n\t%s region %s (rank %u)\n'%('*'*5,region,dbKeys.index(region)))                
        genes=cPickle.loads(regionDB[str(region)])
        testFile.write('\t\tgene rank\tgene name\tstart\tend\tstrand\n')        
        for gene in geneRange:
            testFile.write('\t\t%u\t%s\t%s\t%s\t%s\n'%(gene,genes.IDs[gene],genes.starts[gene],genes.ends[gene],genes.strands[gene]))
        if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))            
        regionDB.close()        
        
def test_exons_by_predicted_transcript(action,species,testFile='',tscriptIDs='',tscriptRange='',tutoFlag=0):
    if tscriptRange!=''or tscriptIDs!='':
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_exons_by_predicted_transcript.bkdb"            
        tscriptDB=bsddb.btopen(path,'r')
        dbKeys=tscriptDB.keys()        
        for i in range(2):
            write=0
            if i==0:
                if tscriptIDs!='':
                    write=1
            else:
                if tscriptRange!='':
                    tscriptIDs=list()                          
                    for tscriptRank in tscriptRange:
                        tscriptIDs.append(dbKeys[tscriptRank])
                    write=1
            if write:                
                for tscriptID in tscriptIDs:                                                                                                          
                    testFile.write('\t%s predicted transcript %s (rank %u)\n'%('*'*5,tscriptID,dbKeys.index(tscriptID)))            
                    exons=cPickle.loads(tscriptDB[tscriptID])
                    testFile.write('\t\texon.IDs = %s\n'%exons.IDs)
                    testFile.write('\t\texon.indexes = %s\n'%exons.indexes)
                    testFile.write('\t\texon.exonStarts = %s\n'%exons.exonStarts)
                    testFile.write('\t\texon.exonEnds = %s\n'%exons.exonEnds)
                    testFile.write('\t\texon.strands = %s\n'%exons.strands)      
                    testFile.write('\t\n')    
        if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))                                          
        tscriptDB.close()
            
def test_predicted_transcripts_by_region(action,species,testFile='',region='',tscriptIDs='',tscriptRange='',tutoFlag=0):                    
    if region!='' and (tscriptRange!='' or tscriptIDs!=''):
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_predicted_transcripts_by_region.bkdb"
        regionDB=bsddb.btopen(path,'r')
        tscripts=cPickle.loads(regionDB[str(region)])                    
        for i in range(2):
            write=0
            if i==0:
                if tscriptIDs!='':
                    tscriptRange=list()
                    for tscriptID in tscriptIDs:
                        tscriptRange.append(tscripts.IDs.index(tscriptID))
                    write=1
            else:
                if tscriptRange!='':                    
                    write=1
            if write:                                                
                for tscriptRank in tscriptRange:                                
                    testFile.write('\t%s predicted transcript %s\n'%('*'*5,tscripts.IDs[tscriptRank]))
                    testFile.write('\t\tregion = %s\n'%str(region))                              
                    testFile.write('\t\ttranscript.start = %s\n'%tscripts.starts[tscriptRank])
                    testFile.write('\t\ttranscript.end = %s\n'%tscripts.ends[tscriptRank])
                    testFile.write('\t\ttranscript.strand = %s\n'%tscripts.strands[tscriptRank])
                    testFile.write('\t\n')
                if i==0:
                    tscriptRange=''
        if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))               
        regionDB.close()

def test_transcripts_by_gene(action,species,dbType,testFile='',geneIDs='',geneRange='',tutoFlag=0):
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_transcripts_by_gene.bkdb"    
    geneDB=bsddb.btopen(path,'r')
    dbKeys=geneDB.keys()    
    for i in range(2):
        write=0        
        if i==0:
            if geneIDs!='':
                write=1
        else:
            if geneRange!='':
                geneIDs=list()                      
                for geneRank in geneRange:                    
                    geneIDs.append(dbKeys[geneRank])
                write=1              
        if write:                
            for geneID in geneIDs:                       
                #try:                     
                    transcripts=cPickle.loads(geneDB[geneID])                               
                    testFile.write('\t%s gene %s (rank %u)\n'%('*'*5,geneID,dbKeys.index(geneID)))
                    if dbType=='aceview':
                        testFile.write('\t\ttranscript.IDs = %s\n'%transcripts)
                        testFile.write('\t\n')                        
                    elif dbType=='ensembl':                    
                        testFile.write('\t\ttranscript.IDs = %s\n'%transcripts.IDs)
                        if action<=12:
                            testFile.write('\t\n')
                        if action>12:
                            testFile.write('\t\ttranscript.starts = %s\n'%transcripts.starts)
                            testFile.write('\t\ttranscript.ends = %s\n'%transcripts.ends)
                            testFile.write('\t\ttranscript.strands = %s\n'%transcripts.strands)
                            testFile.write('\t\n')                                    
                #except:
                #    testFile.write('\n\t\t%s gene %s does not exist\n'%('*'*5,geneID))
    if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))                   
    geneDB.close()

def test_chromosome_by_region(action,species,testFile='',tutoFlag=0):    
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_chromosome_by_region.bkdb"
    regionDB=bsddb.btopen(path,'r')
    dbKeys=regionDB.keys()
    if tutoFlag:
        while len(dbKeys)>3:
            dbKeys.pop()  
    testFile.write('\t%s region => chromosome\n'%('*'*5))
    for key in dbKeys:            
        testFile.write('\t\t%s => %s\n'%(key,regionDB[key]))
    if tutoFlag:
            testFile.write('\t\t%s\n'%('.'*5))                        
    testFile.write('\t\n')
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_region_by_chromosome.bkdb"
    chromosomeDB=bsddb.btopen(path,'r')    
    dbKeys=chromosomeDB.keys()
    if tutoFlag:
        while len(dbKeys)>3:
            dbKeys.pop()    
    testFile.write('\t%s chromosome => region\n'%('*'*5))
    for key in range(len(dbKeys)):            
        testFile.write('\t\t%s => %s\n'%(dbKeys[key],chromosomeDB[dbKeys[key]]))
    if tutoFlag:
            testFile.write('\t\t%s\n'%('.'*5))
    regionDB.close()        
    chromosomeDB.close()                        
    
def test_transcript_sequence(action,species,dbType,testFile='',tscriptIDs='',tscriptRange='',tutoFlag=0):
    path=os.environ['PYDATA']+"/"+species+"/"+dbType+"/"+species+"_transcript_sequence.bkdb"    
    transcriptDB=bsddb.btopen(path,'r')    
    dbKeys=transcriptDB.keys()
    for i in range(2):
        write=0        
        if i==0:
            if tscriptIDs!='':
                write=1
        else:
            if tscriptRange!='':
                tscriptIDs=list()                      
                for tscriptRank in tscriptRange:                    
                    tscriptIDs.append(dbKeys[tscriptRank])
                write=1              
        if write:                
            for tscriptID in tscriptIDs:                                                
                transcript=cPickle.loads(transcriptDB[tscriptID])
                testFile.write('\t%s transcript %s (rank %u)\n'%('*'*5,tscriptID,dbKeys.index(tscriptID)))
                testFile.write('\t\ttranscript.ID=%s\n'%transcript.transcriptID)
                testFile.write('\t\ttranscript.geneID=%s\n'%transcript.geneID)
                testFile.write('\t\ttranscript.start=%s\n'%transcript.start)
                testFile.write('\t\ttranscript.end=%s\n'%transcript.end)
                testFile.write('\t\ttranscript.region=%s\n'%transcript.region)
                testFile.write('\t\ttranscript.chromosome=%s\n'%transcript.chromosome)
                testFile.write('\t\ttranscript.strand=%s\n'%transcript.strand)
                testFile.write('\t\ttranscript.sequence=%s ...\n'%transcript.sequence[:30])    
                testFile.write('\t\n')
    if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))                            
    transcriptDB.close()
       
def test_exons_by_acegene(action,species,testFile='',geneIDs='',geneRange='',tutoFlag=0):
    path=os.environ['PYDATA']+'/'+species+"/aceview/"+species+'_exons_by_gene.bkdb'
    geneDB=bsddb.btopen(path)
    dbKeys=geneDB.keys()
    for i in range(2):
        write=0        
        if i==0:
            if geneIDs!='':
                write=1
        else:
            if geneRange!='':
                geneIDs=list()                      
                for geneRank in geneRange:                    
                    geneIDs.append(dbKeys[geneRank])
                write=1              
        if write:                
            for geneID in geneIDs:
                try:    
                    gene=cPickle.loads(geneDB[geneID])
                    testFile.write('\t%s gene %s (rank %u)\n'%('*'*5,geneID,dbKeys.index(geneID)))                                                                    
                    testFile.write('\t\texonStarts= %s\n'%gene.exonStarts)
                    testFile.write('\t\texonEnds= %s\n'%gene.exonEnds)
                    testFile.write('\t\tstrands= %s\n'%gene.strands)
                    testFile.write('\t\tindexes = %s\n'%gene.indexes)
                    testFile.write('\t\tgroups= %s\n'%gene.groups)
                    testFile.write('\t\tintronStarts= %s\n'%gene.intronStarts)
                    testFile.write('\t\tintronEnds= %s\n'%gene.intronEnds)
                    testFile.write('\t\ttranscriptIDs= %s\n'%gene.transcriptIDs)      
                    testFile.write('\t\ttranscriptStarts= %s\n'%gene.transcriptStarts)
                    testFile.write('\t\ttranscriptEnds= %s\n'%gene.transcriptEnds)  
                    for i in range(len(gene.transcriptsByExon)):
                        testFile.write('\t\texon.transcriptList of exon %u = %s\n'%(i,gene.transcriptsByExon[i]))                        
                except:
                    testFile.write('\n\t\t%s gene %s does not exist\n'%('*'*5,geneID))
                testFile.write('\t\n')    
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))        
    geneDB.close()
  
def test_aceview_genes_by_ensembl_region(action,species,testFile='',regionRank='',geneIDs='',geneRange='',tutoFlag=0):
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_genes_by_ensembl_region.bkdb'
    regionDB=bsddb.btopen(path,'r')
    dbKeys=regionDB.keys()        
    region=cPickle.loads(regionDB[dbKeys[int(regionRank)]])
    for i in range(2):
        write=0
        if i==1:            
            if geneIDs!='':
                geneRange=list()
                for geneID in geneIDs:
                    geneRange.append(region.IDs.index(geneID))
                write=1                
        else:
            if geneRange!='':                              
                write=1
        if write:
            if i==0:    
                testFile.write('\t%s Ensembl region %s\n'%('*'*5,regionRank))
            testFile.write('\t\tgene ranks: %s\n'%repr(geneRange))
            testFile.write('\t\tgene IDs:')            
            for geneRank in geneRange:  
                testFile.write(' %s'%region.IDs[geneRank])
            testFile.write('\t\t\n')    
            testFile.write('\t\tstarts:')
            for geneRank in geneRange:        
                testFile.write(' %s'%repr(region.starts[geneRank]))
            testFile.write('\t\t\n')
            testFile.write('\t\tends:')
            for geneRank in geneRange:    
                testFile.write(' %s'%repr(region.ends[geneRank]))
            testFile.write('\t\t\n')
            testFile.write('\t\tstrands:')        
            for geneRank in geneRange:    
                testFile.write(' %s'%repr(region.strands[geneRank]))
            testFile.write('\t\n\t\n')
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))                            
    regionDB.close()                
    
def test_ensembl_genes_by_gene(action,species,testFile='',geneIDs='',geneRange='',tutoFlag=0):    
    #open ace genes by- ensembl region    
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_ensembl_genes_by_gene.bkdb'          
    geneDB=bsddb.btopen(path,'w')
    dbKeys=geneDB.keys()           
    for i in range(2):
        write=0        
        if i==0:
            if geneIDs!='':
                write=1
        else:
            if geneRange!='':
                geneIDs=list()                      
                for geneRank in geneRange:                    
                    geneIDs.append(dbKeys[geneRank])
                write=1              
        if write:                
            for geneID in geneIDs:               
                try:
                    testFile.write('\t%s gene %s (rank %u)\n'%('*'*5,geneID,dbKeys.index(geneID)))    
                    testFile.write('\t\t%s\n'%repr(cPickle.loads(geneDB[geneID])))
                except:
                    print geneID+'does not exist\n'
                testFile.write('\t\n')
    if tutoFlag:
        testFile.write('\t%s\n'%('.'*5))                
    geneDB.close()                                              
        

def test_chip_positions_by_region(action,species,chipName,region,probeIDs='',probeRange='',testFile='',tutoFlag=0):
    #open position
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_positions_by_region.bkdb"
    regionDB=bsddb.btopen(path)    
    positionList=cPickle.loads(regionDB[str(region)])        
    for i in range(2):
        write=0
        if i==1:            
            if probeIDs!='':
                probeRange=list()
                for probeID in probeIDs:                                 
                    probeRange.append(positionList.probeIDs.index(probeID))
                write=1                                               
        else:                       
            if probeRange!='':                                          
                write=1                
        if write:    
            if i==0:        
                testFile.write('\t%s region %s\n\t\n'%('+'*5,region))    
            for rank in probeRange:            
                testFile.write('\t\t%s position.probeID = %s (rank %s)\n'%('*'*5,positionList.probeIDs[rank],rank))                                    
                testFile.write('\t\tposition.median = %s\n'%str(positionList.medians[rank]))
                testFile.write('\t\tposition.strand: %s\n'%str(positionList.strands[rank]))
                testFile.write('\t\tposition.mismatch: %s\n'%str(positionList.mismatches[rank]))
                testFile.write('\t\n')    
    if tutoFlag:
        testFile.write('\t\t%s\n'%('.'*5))                                        
    regionDB.close()            
            
def test_species_positions_by_region(action,species,region,probeRange='',testFile='',tutoFlag=0):
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_positions_by_region.bkdb"
    regionDB=bsddb.btopen(path,'r')
    positionList=cPickle.loads(regionDB[str(region)])
    
    write=0
    if probeRange!='':                                          
        write=1 
    if write:                               
        testFile.write('\t%s region %s\n\t\n'%('+'*5,region))                    
        for  rank in probeRange:           
            testFile.write('\t\t%s position.probeID = %s (rank %s)\n'%('*'*5,positionList.probeIDs[rank],rank))
            testFile.write('\t\tposition.chips = %s\n'%positionList.chips[rank])                                                
            testFile.write('\t\tposition.median = %s\n'%positionList.medians[rank])
            testFile.write('\t\tposition.strand = %s\n'%positionList.strands[rank])
            testFile.write('\t\tposition.mismatch = %s\n'%positionList.mismatches[rank])
            testFile.write('\t\n')        
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))                
    regionDB.close()    
    

def test_gop_by_region(action,species,chipName,region,gopIDs='',gopRange='',testFile='',tutoFlag=0):    
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_gops_by_region.bkdb" 
    print path       
    regionDB=bsddb.btopen(path,'r')               
    gops=cPickle.loads(regionDB[str(region)])            
    for i in range(2):
        write=0
        if i==1:            
            if gopIDs!='':
                gopRange=list()
                for gopID in gopIDs:
                    try:                                                             
                        gopRange.append(int(where(gops.IDs==gopID)[0][0]))
                    except:
                        pass                                        
                write=1                                                               
        else:                       
            if gopRange!='':                                          
                write=1                 
        if write:                       
            if i==0:        
                testFile.write('\t%s region %s\n\t\n'%('+'*5,region)) 
            print gopRange               
            for rank in gopRange:                                    
                testFile.write('\t\t%s gop.ID %s (rank %u)\n'%('+'*5,gops.IDs[rank],rank))                    
                testFile.write('\t\tgop.start = %s\n'%gops.starts[rank])
                testFile.write('\t\tgop.end = %s\n'%gops.ends[rank])
                testFile.write('\t\tgop.strand = %s\n'%gops.strands[rank])
                testFile.write('\t\tgop.probeNb = %s\n'%gops.positionNbs[rank])                                
                try:
                    testFile.write('\t\tgop.upGeneIDs = %s\n'%gops.upGeneIDs[rank])
                    testFile.write('\t\tgop.upGeneDistances = %s\n'%gops.upGeneDistances[rank])
                except:
                    testFile.write('\t\tgop.upGeneIDs =-\n')
                    testFile.write('\t\tgop.upGeneDistances =-\n')
                try:     
                    testFile.write('\t\tgop.downGeneIDs = %s\n'%gops.downGeneIDs[rank])
                    testFile.write('\t\tgop.downGeneDistances = %s\n'%gops.downGeneDistances[rank])
                except:
                    testFile.write('\t\tgop.downGeneIDs = -\n')
                    testFile.write('\t\tgop.downGeneDistances =-\n')
    if tutoFlag:
        testFile.write('\n\t%s\n'%('.'*5))                                        
    regionDB.close()    
 
    
if __name__=='__main__':    
    #control argument list
    if len(sys.argv)==1:
        sys.exit("pscontrol_chip needs at least actions parameters")
     
    dbType=''    
    chipIDs=''
    chipRange=''
    gopIDs=''
    gopRange=''
    myChipName=''
    mySpeciesName=''
    speciesIDs=''
    speciesRange=''             
    probeSetIDs=''
    probeSetRange=''
    probeSetRank=''    
    probeIDs=''
    probeRange=''    
    exonRange=''
    geneIDs=''
    geneRange=''
    exonIDs=''    
    region=''    
    tscriptIDs=''
    tscriptRange=''
    tutoFlag=0
    writeFlag=0    
    print gopIDs
    usage='test uses -a -b -c -i -e -q -r -s -t -w -h or --action --sub --chip --ids -ensembl --range1 --range2 --species --tuto --write --help'
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:i:e:q:r:s:t:w:h", ["action=","sub=","chip=","ids=","ensembl=","species=","range1=","range2=","tuto=","write=","help"])        
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print(usage)
        sys.exit(2)            
    for opt, arg in opts:        
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit()
        elif opt in ("-a", "--action"):
            action= int(arg)
        elif opt in ("-b", "--sub"):
            subAction=arg    
        elif opt in ("-c", "--chip"):
            myChipName=arg
            """
            Recover Ensembl chip name which is in general equal to  
            chip.shortName, the bracketed name of the corresponding
            GPL platform in GEO.        
            """                
            path=os.environ['COMMONDB']+"/chip.bkdb"
            chip=bsddb.btopen(path,'r')                
            currChip=cPickle.loads(chip[myChipName])
            mySpeciesName=currChip.mySpeciesName        
        elif opt in ("-i", "--ids"):
            if action==1:
                chipIDs=eval(arg)
            elif action==2:
                speciesIDs=eval(arg)            
            elif action==3:                            
                probeSetIDs=eval(arg)
            elif action==4:                      
                probeIDs=eval(arg)
            elif action==5:                            
                probeSetIDs=eval(arg)
            elif action==6:                                                                           
                probeSetID=[arg]
            elif action==7:                                                                           
                exonIDs=eval(arg)            
            elif action==8:                        
                geneIDs=eval(arg)
            elif action==9:                        
                geneIDs=eval(arg)
            elif action==10:                        
                tscriptIDs=eval(arg)
            elif action==11:                        
                tscriptIDs=eval(arg)
            elif action==12:                        
                geneIDs=eval(arg)        
            elif action==14: 
                tscriptIDs=eval(arg)
            elif action == 15:                        
                tscriptIDs=eval(arg)                                       
            elif action==16:
                if subAction=='a': 
                    geneIDs=eval(arg)
                elif subAction=='b': 
                    geneIDs=eval(arg)    
                elif subAction=='c': 
                    tscriptIDs=eval(arg)
                elif subAction=='d': 
                    geneIDs=eval(arg)    
                elif subAction=='e': 
                    exonIDs=eval(arg)
            elif action==17:                        
                geneIDs=eval(arg)    
            elif action==18:                        
                probeIDs=eval(arg)                                                   
            elif action==19:                        
                probeIDs=eval(arg)                
            elif action==20:                        
                gopIDs=eval(arg)
            elif action==22:                        
                gopIDs=eval(arg)                
            elif action in [21,23,24,25]:
                if subAction=='a':
                    geneIDs=eval(arg)                    
                elif subAction=='b':                
                    probeSetID=[arg]                                                        
        elif opt in ("-q","--range1"):
            if action==4:           
                probeRange=eval(arg)
            elif action==9:
                region=arg                
            elif action==11:
                region=arg                    
            elif action==16:
                region=arg                
            elif action==18:
                region=arg                
            elif action==19:
                region=arg                
            elif action==20:
                region=arg                                                                    
        elif opt in ("-r","--range2"):            
            if arg!='':
                probeSetRange=eval(arg)
                if action==1:
                    chipRange=eval(arg)
                elif action==2:
                    speciesRange=eval(arg)
                elif action==3:                            
                    probeSetRange=eval(arg)            
                elif action==5:                            
                    probeSetRange=eval(arg)
                elif action==6:                                                                           
                    probeSetRank=int(arg)
                elif action==7:
                    exonRange=eval(arg)                                
                elif action==8:
                    geneRange=eval(arg)                
                elif action==9:
                    geneRange=eval(arg)                        
                elif action==10:
                    tscriptRange=eval(arg)
                elif action==11:
                    tscriptRange=eval(arg)
                elif action==12:
                    geneRange=eval(arg)   
                elif action==14:
                    tscriptRange=eval(arg)     
                elif action==15:
                    tscriptRange=eval(arg)
                elif action==16:
                    if subAction=='a':
                        geneRange=eval(arg)
                    elif subAction=='b':
                        geneRange=eval(arg)
                    elif subAction=='c':
                        tscriptRange=eval(arg)
                    elif subAction=='d':                
                        geneRange=eval(arg)
                    elif subAction=='e':                
                        exonRange=eval(arg)
                elif action==17:           
                    geneRange=eval(arg)
                elif action==18:
                    probeRange=eval(arg)     
                elif action==19:
                    probeRange=eval(arg)
                elif action==20:
                    gopRange=eval(arg)
                elif action in [21,23,24,25]:
                    if subAction=='a': 
                        geneRange=eval(arg)
                    elif subAction=='b':
                        probeSetRank=eval(arg)    
        elif opt in ("-s", "--species"):            
            mySpeciesName=arg                        
        elif opt in ("-t","--tuto"):
            tutoFlag=int(arg)    
        elif opt in ("-w","--write"):
            writeFlag=int(arg)    
        elif opt in ("-e", "--ensembl"):
            if myChipName=='':
                print("chip must be specified if ensembl is used")
                sys.exit()            
            ensVersion=arg            
            if int(ensVersion.split("_")[0])<48:
                ensChipName=currChip.ens47Name
                ensPort=3306
                martPort=3316
            else:
                ensChipName=currChip.ens48Name
                ensPort=5306
                martPort=5316
            
        elif opt in ("-v", "--aceview"):
            aceVersion=arg                                     
        else:            
            assert False, "unhandled option %s"%opt
    
    if mySpeciesName!='':        
        """
        Recover the species name used in Ensembl, which is
        the lower case official species name with underscores
        """                            
        path=os.environ['COMMONDB']+"/species.bkdb"        
        speciesDB=bsddb.btopen(path,'r')                
        currSpeciesName=cPickle.loads(speciesDB[mySpeciesName])            
        ensSpeciesName=currSpeciesName.officialName.lower()
        ensSpeciesName=ensSpeciesName.split(" ")[0][0]+ensSpeciesName.split(" ")[1]
        ensSpeciesName=ensSpeciesName.lower()
        
    os.chdir(os.environ['PYDATA'])    
    if writeFlag:
        testFile=open('test.txt','a')
    elif tutoFlag:
        testFile=open('tuto.txt','a')                
    else:        
        testFile=sys.stdout 
                        
    if tutoFlag:
        if action ==1:
            if len(chipIDs)+len(chipRange)<=3:
                tutoFlag=0
            while len(chipIDs)>3:
                chipIDs.pop()
            if len(chipIDs)<3:
                while len(chipRange)>3-len(chipIDs):
                    chipRange.pop()
            else:
                chipRange=''
        
        if action ==2:    
            if len(speciesIDs)+len(speciesRange)<=3:
                tutoFlag=0
            while len(speciesIDs)>3:
                speciesIDs.pop()
            if len(speciesIDs)<3:
                while len(speciesRange)>3-len(speciesIDs):
                    speciesRange.pop()
            else:
                speciesRange=''
        
        if action in [3,5,21,23,24,25]:
            limit=3    
            if action in [21,23,24,25]:
                limit=1
            else:
                if len(probeSetIDs)+len(probeSetRange)<=3:
                    tutoFlag=0                                            
            while len(probeSetIDs)>limit:
                probeSetIDs.pop()
            if len(probeSetIDs)<limit:
                while len(probeSetRange)>limit-len(probeSetIDs):
                    probeSetRange.pop()
            else:
                probeSetRange=''
        
        if action in [4,18,19,23,25]: 
            if action in [4,18,19]:               
                if len(probeIDs)+len(probeRange)<=3:
                    tutoFlag=0
            while len(probeIDs)>3:
                probeIDs.pop()
            if len(probeIDs)<3:
                while len(probeRange)>3-len(probeIDs):
                    probeRange.pop()
            else:
                probeRange=''
        
        if action in [7,16]:    
            if len(exonIDs)+len(exonRange)<=3:
                tutoFlag=0
            while len(exonIDs)>3:
                exonIDs.pop()
            if len(exonIDs)<3:
                while len(exonRange)>3-len(exonIDs):
                    exonRange.pop()
            else:
                exonRange=''
                     
        if action in [8,9,12,16,17,21,24]:
            limit=3    
            if action in [8,21,24]:
                limit=1  
            else:            
                if len(geneIDs)+len(geneRange)<=limit:
                    tutoFlag=0
            while len(geneIDs)>limit:
                geneIDs.pop()
            if len(geneIDs)<limit:
                while len(geneRange)>limit-len(geneIDs):
                    geneRange.pop()
            else:
                geneRange=''
                
        if action in [10,11,14,15,16]:    
            if len(tscriptIDs)+len(tscriptRange)<=3:
                tutoFlag=0
            while len(tscriptIDs)>3:
                tscriptIDs.pop()
            if len(tscriptIDs)<3:
                while len(tscriptRange)>3-len(tscriptIDs):
                    tscriptRange.pop()
            else:
                tscriptRange=''
                
        if action==20:    
            if len(gopIDs)+len(gopRange)<=3:
                tutoFlag=0
            while len(gopIDs)>3:
                gopIDs.pop()
            if len(gopIDs)<3:
                while len(gopRange)>3-len(gopIDs):
                    gopRange.pop()
            else:
                gopRange=''    
                                                        
       
    if action==1:                
        test_chip(chipIDs=chipIDs,chipRange=chipRange,testFile=testFile,tutoFlag=tutoFlag)
    if action==2:                        
        test_species(speciesIDs=speciesIDs,speciesRange=speciesRange,testFile=testFile,tutoFlag=tutoFlag)
    if action==3:                     
        test_probeset(3,mySpeciesName,myChipName,testFile=testFile,probeSetIDs=probeSetIDs,probeSetRange=probeSetRange,tutoFlag=tutoFlag)        
    if action==4:                                                                 
        test_probe(4,mySpeciesName,myChipName,'ensembl',testFile=testFile,probeIDs=probeIDs,tutoFlag=tutoFlag)                
    if action==5:                                                             
        test_probeset(5,mySpeciesName,myChipName,testFile=testFile,probeSetIDs=probeSetIDs,probeSetRange=probeSetRange)                        
    if action==6:    
        path=os.environ['COMMONDB']+"/chip_probeset.bkdb"           
        chipDB=bsddb.btopen(path)                                                        
        if probeSetRank!='':        
            probeSetIDs=cPickle.loads(chipDB[myChipName])                                            
            probeSetID=[probeSetIDs[probeSetRank]]                
        test_probeset(6,mySpeciesName,myChipName,testFile=testFile,probeSetIDs=probeSetID)
        chipDB.close()        
        path=os.environ['PYDATA']+"/"+mySpeciesName+"/"+myChipName+"_probeset.bkdb"    
        probeSetDB=bsddb.btopen(path)
        probeSet=cPickle.loads(probeSetDB[probeSetID[0]])
        probeSetDB.close()
        testFile.write('\t\n')
        if tutoFlag:
            probeNb=3
        else:
            probeNb=len(probeSet.probeIDs)        
        for probe in probeSet.probeIDs[:probeNb]:
            test_probe(6,mySpeciesName,myChipName,'ensembl',testFile=testFile,probeIDs=[probe])
        if tutoFlag:
            testFile.write('\n\t%s\n'%('.'*5))                                    
    if action==7:
        test_transcripts_by_exon(7,mySpeciesName,'ensembl',testFile=testFile,exonIDs=exonIDs,exonRange=exonRange,tutoFlag=tutoFlag)                            
    if action==8:            
        test_exons_by_gene(8,mySpeciesName,testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=1)
    if action==9:              
        test_genes_by_region(9,mySpeciesName,testFile=testFile,geneIDs=geneIDs,region=region,geneRange=geneRange,tutoFlag=tutoFlag)    
    if action==10:            
        test_exons_by_predicted_transcript(10,mySpeciesName,testFile=testFile,tscriptIDs=tscriptIDs,tscriptRange=tscriptRange,tutoFlag=tutoFlag)    
    if action==11:               
        test_predicted_transcripts_by_region(11,mySpeciesName,testFile=testFile,region=region,tscriptIDs=tscriptIDs,tscriptRange=tscriptRange,tutoFlag=tutoFlag)    
    if action==12:                    
        test_transcripts_by_gene(12,mySpeciesName,'ensembl',testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)        
    if action==13:                    
        test_chromosome_by_region(13,mySpeciesName,testFile=testFile,tutoFlag=tutoFlag)    
    if action==14:                        
        test_transcript_sequence(14,mySpeciesName,'ensembl',testFile=testFile,tscriptIDs=tscriptIDs,tscriptRange=tscriptRange,tutoFlag=tutoFlag)    
    if action==15:                            
        test_transcript_sequence(15,mySpeciesName,'aceview',testFile=testFile,tscriptIDs=tscriptIDs,tscriptRange=tscriptRange)        
    if action==16:
        if subAction=='a':                            
            test_exons_by_acegene(16,mySpeciesName,testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)        
        elif subAction=='b':        
            test_aceview_genes_by_ensembl_region(16,mySpeciesName,testFile=testFile,regionRank=region,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)    
        elif subAction=='c':                
            test_transcript_sequence(16,mySpeciesName,'aceview',testFile=testFile,tscriptIDs=tscriptIDs,tscriptRange=tscriptRange,tutoFlag=tutoFlag)
        elif subAction=='d':                               
            test_transcripts_by_gene(16,mySpeciesName,'aceview',testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)
        elif subAction=='e':                        
            test_transcripts_by_exon(16,mySpeciesName,'aceview',testFile=testFile,exonIDs=exonIDs,exonRange=exonRange,tutoFlag=tutoFlag)                                
    if action==17:        
        test_ensembl_genes_by_gene(17,mySpeciesName,testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)
    if action==18:            
        test_chip_positions_by_region(18,mySpeciesName,myChipName,region,probeIDs=probeIDs,probeRange=probeRange,testFile=testFile,tutoFlag=tutoFlag)
    if action==19:        
        test_species_positions_by_region(19,mySpeciesName,region,probeRange=probeRange,testFile=testFile,tutoFlag=tutoFlag)        
    if action==20:        
        test_gop_by_region(20,mySpeciesName,myChipName,region,gopIDs=gopIDs,gopRange=gopRange,testFile=testFile,tutoFlag=tutoFlag)    
    if action in [21,23,24,25]:
        if subAction=='a':
            if action==21:
                test_exons_by_gene(action,mySpeciesName,testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)
            else:#action = 24                
                test_exons_by_acegene(action,mySpeciesName,testFile=testFile,geneIDs=geneIDs,geneRange=geneRange,tutoFlag=tutoFlag)                
        if subAction=='b':
            #action = 21,23,24 or 25
            path=os.environ['COMMONDB']+"/chip_probeset.bkdb"           
            chipDB=bsddb.btopen(path)                                                        
            if probeSetRank!='':        
                probeSetIDs=cPickle.loads(chipDB[myChipName])
                probeSetID=[probeSetIDs[probeSetRank]]                                                                                                  
            test_probeset(action,mySpeciesName,myChipName,testFile=testFile,probeSetIDs=probeSetID)
            chipDB.close()       
            if action in [23,25]: 
                path=os.environ['PYDATA']+"/"+mySpeciesName+"/"+myChipName+"_probeset.bkdb"    
                probeSetDB=bsddb.btopen(path)
                probeSet=cPickle.loads(probeSetDB[probeSetID[0]])
                testFile.write('\t\n')
                if tutoFlag:
                    probeNb=3
                else:
                    probeNb=len(probeSet.probeIDs)                        
                for probe in probeSet.probeIDs[:probeNb]:
                    if action==21:
                        test_probe(action,mySpeciesName,myChipName,'ensembl',testFile=testFile,probeIDs=[probe])
                    else:                                        
                        test_probe(action,mySpeciesName,myChipName,'aceview',testFile=testFile,probeIDs=[probe])