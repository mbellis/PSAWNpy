#!/usr/bin/python
"""**Classes used by psawnPy.**"""
from numpy import array,uint32,int8,uint8

class GeneList(object):
    """
    **Define a list of genes.**
    ::
    
        Data attributes:
            IDs        
            starts
            ends
            strands
    
    """
    
    def __init__(self,idList=[],startArray=array((),dtype=uint32),endArray=array((),dtype=uint32),strandArray=array((),dtype=int8)):
        self.IDs=idList        
        self.starts=startArray
        self.ends=endArray
        self.strands=strandArray

class GopList(object):
    """
    **Define a list of group of probes (GOP).**
    ::
    
        Data attributes:
            IDs
            region        
            starts
            ends
            strands
            upGeneDistances
            downGeneDistances
            upGeneIDs
            downGeneIDs        
            positionNbs
    
    """
        
    def __init__(self,idArray=array(()),regionArray=array(()),startArray=array((),dtype=uint32),endArray=array((),dtype=uint32),strandArray=array((),dtype=int8),upGeneDistanceArray=array((),dtype=uint8),downGeneDistanceArray=array((),dtype=uint8),upGeneIDArray=array(()),downGeneIDArray=array(()),positionNbArray=array((),dtype=uint8)):
        self.IDs=idArray
        self.region=regionArray        
        self.starts=startArray
        self.ends=endArray
        self.strands=strandArray
        self.upGeneDistances=upGeneDistanceArray
        self.downGeneDistances=downGeneDistanceArray
        self.upGeneIDs=upGeneIDArray
        self.downGeneIDs=downGeneIDArray        
        self.positionNbs=positionNbArray
                     
class StructureList(object):
    """
    **Define a list of genomic structures.**
    ::
    
        Data attributes:    
            IDs=idList        
            starts=startArray
            ends=endArray
            strands=strandArray
        
    """
    
    def __init__(self,idList=[],startArray=array((),dtype=uint32),endArray=array((),dtype=uint32),strandArray=array((),dtype=int8)):
        self.IDs=idList        
        self.starts=startArray
        self.ends=endArray
        self.strands=strandArray 
         
class ExonList(object):
    """
    **Define a list of exons.**
    ::
    
        data attributes:
            exonID ------ list of ensembl exon stable id
            exonStart --- list of exon start positions
            exonEnd ----- list of exon end positions
            exonStrand -- list of strand orientation        
            groups
            intronStarts
            intronEnds
            transcriptsByExon
            transcriptIDs      
            transcriptStarts
            transcriptEnds
            indexes
            
    """
    
    def __init__(self,exonIndexes=[],idList=[],groupList=[],exonStartArray=array((),dtype=uint32),\
        exonEndArray=array((),dtype=uint32),strandArray=array((),dtype=int8),\
        intronStartArray=array((),dtype=uint32),intronEndArray=array((),dtype=uint32),\
        transcriptsByExon=[],transcriptIDs=[],transcriptStarts=[],transcriptEnds=[]):        
        self.IDs=idList        
        self.exonStarts=exonStartArray
        self.exonEnds=exonEndArray
        self.strands=strandArray
        self.groups=groupList
        self.intronStarts=intronStartArray
        self.intronEnds=intronEndArray
        self.transcriptsByExon=transcriptsByExon
        self.transcriptIDs=transcriptIDs      
        self.transcriptStarts=transcriptStarts
        self.transcriptEnds=transcriptEnds
        if exonIndexes==[]:
            if self.strands[0]==1:
                self.indexes=range(len(self.strands))                    
            else:
                self.indexes=range(len(self.strands)-1,-1,-1)
        else:
            self.indexes=exonIndexes
        
class TranscriptList(object):           
    """
    **Define a list of transcripts.**
    ::
    
        Data attributes:
            IDs=idList
            starts=startArray
            ends=endArray
            strands=strandArray
        
    """
    
    def __init__(self,idList=[],startArray=array((),dtype=uint32),endArray=array((),dtype=uint32),strandArray=array((),dtype=int8)):
        self.IDs=idList
        self.starts=startArray
        self.ends=endArray
        self.strands=strandArray
        
class PositionList(object):
    """
    **Define a list of genomic positions.**
    ::
    
        Data attributes:    
            probeIDs=probeIDList
            chips=chipsList                
            medians=medianArray        
            strands=strandArray
            mismatches=mismatchNbArray
        
    """
    
    def __init__(self,probeIDList=[],chipsList=[],medianArray=array((),dtype=uint32),strandArray=array((),dtype=int8),mismatchNbArray=array((),dtype=uint8)):
        self.probeIDs=probeIDList
        self.chips=chipsList                
        self.medians=medianArray        
        self.strands=strandArray
        self.mismatches=mismatchNbArray

class TargetingProbeset(object):
    """
    **Define a dictionnary of targets** indexed on the number of targeting probes**
    ::
    
        Data attributes:
            probesetNames
            probesetIndexes
            targetedExons
            targetedGroups
            targetedTranscripts
            notTargetedTranscripts
            inGeneProbeNbs
            notInExonProbeNbs
    
    """
    
    def __init__(self,probesetNames={},probesetIndexes={},targetedExons={},targetedGroups={},targetedTranscripts={},notTargetedTranscripts={},inGeneProbeNbs={},notInExonProbeNbs={}):
        self.probesetNames=probesetNames
        self.probesetIndexes=probesetIndexes
        self.targetedExons=targetedExons
        self.targetedGroups=targetedGroups
        self.targetedTranscripts=targetedTranscripts
        self.notTargetedTranscripts=notTargetedTranscripts
        self.inGeneProbeNbs=inGeneProbeNbs
        self.notInExonProbeNbs=notInExonProbeNbs
                                
class TargetedGene(object):
    """
    **Define a list of probes targeting a gene.**
    ::
    
        Data attributes:    
            probeIndexes
            probeLocalisations
            firstStructureIndexes
            firstStructureGroups
            sndStructureIndexes
            sndStructureGroups
            firstStructureIDs
            probePositions
            probeRepetitionNbs
            probeStrands
            probeMismatchNbs
            inExonProbeNb
            inSpliceProbeNb
            inIntronProbeNb
            upProbeNb
            downProbeNb  
            exonSet                       
            groupSet
            targetedTranscripts        
            notTargetedTranscripts
            ensemblGenes
         
        Function:
            fill_info     
    
    """
    
    def __init__(self, probeIndexes=[],probeLocalisations=[],firstStructureIndexes=[],firstStructureGroups=[],sndStructureIndexes=[],sndStructureGroups=[],firstStructureIDs=[],probePositions=[],probeMismatchNbs=[],probeRepetitionNbs=[],probeStrands=0,inExonProbeNb=0,inSpliceProbeNb=0,inIntronProbeNb=0,upProbeNb=0,downProbeNb=0,targetedTranscripts=None,exonSet=None,groupSet=None,notTargetedTranscripts=None, ensemblGenes=None):
        self.probeIndexes=probeIndexes
        self.probeLocalisations=probeLocalisations
        self.firstStructureIndexes=firstStructureIndexes
        self.firstStructureGroups=firstStructureGroups
        self.sndStructureIndexes=sndStructureIndexes
        self.sndStructureGroups=sndStructureGroups
        self.firstStructureIDs=firstStructureIDs
        self.probePositions=probePositions
        self.probeRepetitionNbs=probeRepetitionNbs
        self.probeStrands=probeStrands
        self.probeMismatchNbs=probeMismatchNbs
        self.inExonProbeNb=inExonProbeNb
        self.inSpliceProbeNb=inSpliceProbeNb
        self.inIntronProbeNb=inIntronProbeNb
        self.upProbeNb=upProbeNb
        self.downProbeNb=downProbeNb  
        self.exonSet=exonSet                       
        self.groupSet=groupSet
        self.targetedTranscripts=targetedTranscripts        
        self.notTargetedTranscripts=notTargetedTranscripts
        self.ensemblGenes=ensemblGenes        
        
    def fill_info(self):
        probeNb=len(self.probeIndexes)
        self.inExonProbeNb=[0]*probeNb            
        self.inSpliceProbeNb=[0]*probeNb
        self.inIntronProbeNb=[0]*probeNb
        self.upProbeNb=[0]*probeNb
        self.downProbeNb=[0]*probeNb
        self.exonSet=set()
        self.groupSet=set()
        for probeIndex in range(len(self.probeLocalisations)):
            if self.probeLocalisations[probeIndex]==[]:
                    self.probeLocalisations[probeIndex]=['o']
            else:        
                for j in range(len(self.probeLocalisations[probeIndex])):
                    currStructure=self.probeLocalisations[probeIndex][j]                        
                    if ['e','i'].count(currStructure):                                            
                        self.exonSet.add(self.firstStructureIndexes[probeIndex][j])                        
                        self.groupSet.add(self.firstStructureGroups[probeIndex][j])                                            
                    if self.probeLocalisations[probeIndex][j]==[]:
                        self.probeLocalisations[probeIndex][j]=['o']                             
                    if currStructure=='e':
                        self.inExonProbeNb[probeIndex]+=self.probeRepetitionNbs[probeIndex][j]
                    elif currStructure=='s':
                        self.inSpliceProbeNb[probeIndex]+=self.probeRepetitionNbs[probeIndex][j]
                    elif currStructure=='i':
                        self.inIntronProbeNb[probeIndex]+=self.probeRepetitionNbs[probeIndex][j]
                    elif currStructure=='u':
                        self.upProbeNb[probeIndex]+=self.probeRepetitionNbs[probeIndex][j]
                    elif currStructure=='d':
                        self.downProbeNb[probeIndex]+=self.probeRepetitionNbs[probeIndex][j]             
       
class Chip(object):
    """            
    **define a Chip model.**
    ::
    
        Data attributes:    
            myName --------- local chip name ('m2')        
            name ----------- chip name ('Human Genome U95A Array')
            shortName ------ chip short name ('HG_U95A') 
            mySpeciesName -- species name used inside program 
            probesetNb------ number of probe sets           
            probeNb--------- median probe number per (pair of) probe set
            compName ------- compagny name ('Affymetrix')
            ens47Name ------ chip name in Ensembl version 1 to 47
            ens48Name ------ chip name from Ensembl version 48
            geoName -------- GEO platform name (GPL)
            
    """
    
    def __init__(self,myName=None,name=None,shortName=None,mySpeciesName=None,probesetNb=None,probeNb=None,compName=None,ens47Name=None,ens48Name=None,geoName=None):
        self.myName=myName        
        self.name=name
        self.shortName=shortName        
        self.mySpeciesName=mySpeciesName
        self.probesetNb=int(probesetNb)        
        self.probeNb=int(probeNb)
        self.compName=compName
        self.ens47Name=ens47Name
        self.ens48Name=ens48Name
        self.geoName=geoName

class Species(object):
    """
    **Define a Species.** 
    ::
       
        Data attributes:    
            myName ----- -- species name used inside program
            officialName -- scientific latin name
            ensName ------- name used in Ensembl tables
             
    """
    
    def __init__(self,myName=None,officialName=None,ensName=None):
        self.myName=myName
        self.officialName=officialName
        self.ensName=ensName
                                        
class Transcript(object):
    """
    **Define a transcript.**
    ::
    
        Data attributes:
            transcriptID -- transcript ID
            geneID -------- gene ID
            start --------- genomic start position of transcript
            end ----------- genomic end position of transcript
            region -------- Ensembl region
            chromosome ---- chromosome        
            sequence ------ transcript sequence
            strand -------- gene strandedness
    
    """
    def __init__(self,ID=None,geneID=None,start=None,end=None,region=None,chromosome=None,sequence=None,strand=None):
        self.transcriptID=ID
        self.geneID=geneID
        self.start=start
        self.end=end
        self.region=region
        self.chromosome=chromosome        
        self.sequence=sequence
        self.strand=strand
        
class Probe(object):
    """
    **Define a probe.**
    ::
    
        Data attributes:
            probeID
            xPosition
            yPosition
            ensemblID        
            probesetID        
            ensProbesetID
            index
            sequence
            targetPosition               
            genes4exon
            genes4intron
            genes4splice        
            genes4up
            genes4down
            outOfGeneNb
            notInSequence
    
    """
    
    def __init__(self,probeID=None,xPosition=None,yPosition=None,ensemblID=None,probesetID=None,ensProbesetID=None,index=None,sequence=None,targetPosition=None,genes4exon={},genes4splice={},genes4intron={},genes4up={},genes4down={},outOfGeneNb=0,notInSequence=None):
        self.probeID=probeID
        self.xPosition=xPosition
        self.yPosition=yPosition
        self.ensemblID=ensemblID        
        self.probesetID=probesetID        
        self.ensProbesetID=ensProbesetID
        self.index=index
        self.sequence=sequence
        self.targetPosition=targetPosition               
        self.genes4exon=genes4exon
        self.genes4intron=genes4intron
        self.genes4splice=genes4splice        
        self.genes4up=genes4up
        self.genes4down=genes4down
        self.outOfGeneNb=outOfGeneNb
        self.notInSequence=notInSequence
        
            
class ProbesetList(object):
    """
    **Define a list of probe sets.**
    ::    
    
        Data attributes:
            probesetNb=len(probesetID)
            probesetID=probesetID
            probesetIndex=probesetIndex
            probeNb=probeNb
            ensGeneIDnb=ensGeneIDnb
            affyGeneIDnb=affyGeneIDnb                
            ensExonGeneNbs=ensExonGeneNbs
            ensIntronGeneNbs=ensIntronGeneNbs
            ensUpGeneNbs=ensUpGeneNbs
            ensDownGeneNbs=ensDownGeneNbs
            ensOutProbeNbs=ensOutProbeNbs
            ensNisProbeNbs=ensNisProbeNbs        
            ourEnsGeneNb=ourEnsGeneNb
            ourEnsProbeNb=ourEnsProbeNb
            ensGeneIDs4EnsNb=ensGeneIDs4EnsNb
            ourEnsGeneIDs4EnsNb=ourEnsGeneIDs4EnsNb
            commonGeneIDs4EnsNb=commonGeneIDs4EnsNb
            ourEnsTargetedExons=ourEnsTargetedExons
            ourEnsTargetedGroups=ourEnsTargetedGroups
            ourEnsTargetedTranscripts=ourEnsTargetedTranscripts
            ourEnsNotTargetedTranscripts=ourEnsNotTargetedTranscripts          
            aceExonGeneNbs=aceExonGeneNbs
            aceIntronGeneNbs=aceIntronGeneNbs
            aceUpGeneNbs=aceUpGeneNbs
            aceDownGeneNbs=aceDownGeneNbs
            aceOutProbeNbs=aceOutProbeNbs
            aceNisProbeNbs=aceNisProbeNbs              
            ourAceGeneNb=ourAceGeneNb
            ourAceProbeNb=ourAceProbeNb                            
            ensGeneIDs4AceNb=ensGeneIDs4AceNb
            ourEnsGeneIDs4AceNb=ourEnsGeneIDs4AceNb
            commonGeneIDs4AceNb=commonGeneIDs4AceNb        
            ourAceTargetedExons=ourAceTargetedExons
            ourAceTargetedGroups=ourAceTargetedGroups
            ourAceTargetedTranscripts=ourAceTargetedTranscripts
            ourAceNotTargetedTranscripts=ourAceNotTargetedTranscripts
        
    """
    
    def __init__(self,probesetID=None,probesetIndex=None,probeNb=None,ensGeneIDnb=None,affyGeneIDnb=None,\
                 ensExonGeneNbs=None,ensIntronGeneNbs=None,ensUpGeneNbs=None,ensDownGeneNbs=None,ensOutProbeNbs=None,ensNisProbeNbs=None,\
                 aceExonGeneNbs=None,aceIntronGeneNbs=None,aceUpGeneNbs=None,aceDownGeneNbs=None,aceOutProbeNbs=None,aceNisProbeNbs=None,\
                 ourEnsGeneNb=None,ourEnsProbeNb=None,\
                 ourAceGeneNb=None,ourAceProbeNb=None,\
                 ensGeneIDs4EnsNb=None, ourEnsGeneIDs4EnsNb=None, commonGeneIDs4EnsNb=None,\
                 ourEnsGeneIDs4AceNb=None, ensGeneIDs4AceNb=None, ourEnsGeneIDs4AceNB=None, commonGeneIDs4AceNb=None,\
                 ourEnsTargetedExons=None,\
                 ourEnsTargetedGroups=None,\
                 ourEnsTargetedTranscripts=None,ourEnsNotTargetedTranscripts=None,\
                 ourAceTargetedExons=None,\
                 ourAceTargetedGroups=None,\
                 ourAceTargetedTranscripts=None,ourAceNotTargetedTranscripts=None):                 
                 
        self.probesetNb=len(probesetID)
        self.probesetID=probesetID
        self.probesetIndex=probesetIndex
        self.probeNb=probeNb
        self.ensGeneIDnb=ensGeneIDnb
        self.affyGeneIDnb=affyGeneIDnb                
        self.ensExonGeneNbs=ensExonGeneNbs
        self.ensIntronGeneNbs=ensIntronGeneNbs
        self.ensUpGeneNbs=ensUpGeneNbs
        self.ensDownGeneNbs=ensDownGeneNbs
        self.ensOutProbeNbs=ensOutProbeNbs
        self.ensNisProbeNbs=ensNisProbeNbs        
        self.ourEnsGeneNb=ourEnsGeneNb
        self.ourEnsProbeNb=ourEnsProbeNb
        self.ensGeneIDs4EnsNb=ensGeneIDs4EnsNb
        self.ourEnsGeneIDs4EnsNb=ourEnsGeneIDs4EnsNb
        self.commonGeneIDs4EnsNb=commonGeneIDs4EnsNb
        self.ourEnsTargetedExons=ourEnsTargetedExons
        self.ourEnsTargetedGroups=ourEnsTargetedGroups
        self.ourEnsTargetedTranscripts=ourEnsTargetedTranscripts
        self.ourEnsNotTargetedTranscripts=ourEnsNotTargetedTranscripts          
        self.aceExonGeneNbs=aceExonGeneNbs
        self.aceIntronGeneNbs=aceIntronGeneNbs
        self.aceUpGeneNbs=aceUpGeneNbs
        self.aceDownGeneNbs=aceDownGeneNbs
        self.aceOutProbeNbs=aceOutProbeNbs
        self.aceNisProbeNbs=aceNisProbeNbs              
        self.ourAceGeneNb=ourAceGeneNb
        self.ourAceProbeNb=ourAceProbeNb                            
        self.ensGeneIDs4AceNb=ensGeneIDs4AceNb
        self.ourEnsGeneIDs4AceNb=ourEnsGeneIDs4AceNb
        self.commonGeneIDs4AceNb=commonGeneIDs4AceNb        
        self.ourAceTargetedExons=ourAceTargetedExons
        self.ourAceTargetedGroups=ourAceTargetedGroups
        self.ourAceTargetedTranscripts=ourAceTargetedTranscripts
        self.ourAceNotTargetedTranscripts=ourAceNotTargetedTranscripts
             
class Probeset(object):
    """
    **Define a probe set.**
    ::
    
        Data attributes:
            probesetID
            probesetIndex
            ensGeneIDs
            affyGeneIDs
            probeNb        
            probeIDs
            ensProbeIDs
            probeIndexes        
            probesetTargetLength
            probesetTargetStart        
            probesetTargetEnd
            ensExonGeneNbs
            ensIntronGeneNbs
            ensUpGeneNbs
            ensDownGeneNbs
            ensOutProbeNbs
            ensNisProbeNbs
            ourEnsGeneIDs
            ourEnsGeneNb
            ourEnsProbeNb
            ensGeneIDs4Ens
            ourEnsGeneIDs4Ens
            commonGeneIDs4Ens
            ensGenesByProbeNb  
            aceExonGeneNbs
            aceIntronGeneNbs
            aceUpGeneNbs
            aceDownGeneNbs
            aceOutProbeNbs
            aceNisProbeNbs      
            ourAceGeneIDs
            ourAceGeneNb
            ourAceProbeNb
            aceGenesByProbeNb
            ourAceToEnsGeneNbs
            ourAceToEnsGeneIDs
            ensGeneIDs4Ace
            ourEnsGeneIDs4Ace
            commonGeneIDs4Ace
            ourGeneIDs
            ourGeneNb
            ourProbeNb        
        
    """
    
    def __init__(self,probesetID=None,ensGeneIDs=None,affyGeneIDs=None,probesetIndex=None,probeIDs=None,\
                 ensProbeIDs=None,probeIndexes=None,probeNb=None,\
                 probesetTargetLength=None, probesetTargetStart=None,probesetTargetEnd=None,\
                 ensExonGeneNbs=None,ensIntronGeneNbs=None,ensUpGeneNbs=None,ensDownGeneNbs=None,ensOutProbeNbs=None,ensNisProbeNbs=None,\
                 aceExonGeneNbs=None,aceIntronGeneNbs=None,aceUpGeneNbs=None,aceDownGeneNbs=None,aceOutProbeNbs=None,aceNisProbeNbs=None,\
                 ourEnsGeneIDs=None, ourEnsGeneNb=None,ourEnsProbeNb=None, ensGenesByProbeNb=None,\
                 ensGeneIDs4Ens=None, ourEnsGeneIDs4Ens=None, commonGeneIDs4Ens=None,\
                 ourAceGeneIDs=None, ourAceGeneNb=None,ourAceProbeNb=None,aceGenesByProbeNb=None,ourAceToEnsGeneIDs=None,\
                 ensGeneIDs4Ace=None, ourEnsGeneIDs4Ace=None, commonGeneIDs4Ace=None,\
                 ourAceToEnsGeneNbs=None,ourGeneIDs=None,ourGeneNb=None,ourProbeNb=None):
        self.probesetID=probesetID
        self.probesetIndex=probesetIndex
        self.ensGeneIDs=ensGeneIDs
        self.affyGeneIDs=affyGeneIDs
        self.probeNb=probeNb        
        self.probeIDs=probeIDs
        self.ensProbeIDs=ensProbeIDs
        self.probeIndexes=probeIndexes        
        self.probesetTargetLength=probesetTargetLength
        self.probesetTargetStart=probesetTargetStart        
        self.probesetTargetEnd=probesetTargetEnd
        self.ensExonGeneNbs=ensExonGeneNbs
        self.ensIntronGeneNbs=ensIntronGeneNbs
        self.ensUpGeneNbs=ensUpGeneNbs
        self.ensDownGeneNbs=ensDownGeneNbs
        self.ensOutProbeNbs=ensOutProbeNbs
        self.ensNisProbeNbs=ensNisProbeNbs
        self.ourEnsGeneIDs=ourEnsGeneIDs
        self.ourEnsGeneNb=ourEnsGeneNb
        self.ourEnsProbeNb=ourEnsProbeNb
        self.ensGeneIDs4Ens=ensGeneIDs4Ens
        self.ourEnsGeneIDs4Ens=ourEnsGeneIDs4Ens
        self.commonGeneIDs4Ens=commonGeneIDs4Ens
        self.ensGenesByProbeNb=ensGenesByProbeNb  
        self.aceExonGeneNbs=aceExonGeneNbs
        self.aceIntronGeneNbs=aceIntronGeneNbs
        self.aceUpGeneNbs=aceUpGeneNbs
        self.aceDownGeneNbs=aceDownGeneNbs
        self.aceOutProbeNbs=aceOutProbeNbs
        self.aceNisProbeNbs=aceNisProbeNbs      
        self.ourAceGeneIDs=ourAceGeneIDs
        self.ourAceGeneNb=ourAceGeneNb
        self.ourAceProbeNb=ourAceProbeNb
        self.aceGenesByProbeNb=aceGenesByProbeNb
        self.ourAceToEnsGeneNbs=ourAceToEnsGeneNbs
        self.ourAceToEnsGeneIDs=ourAceToEnsGeneIDs
        self.ensGeneIDs4Ace=ensGeneIDs4Ace
        self.ourEnsGeneIDs4Ace=ourEnsGeneIDs4Ace
        self.commonGeneIDs4Ace=commonGeneIDs4Ace
        self.ourGeneIDs=ourGeneIDs
        self.ourGeneNb=ourGeneNb
        self.ourProbeNb=ourProbeNb        
        
       

class Position(object):
    """
    **Define a position.**
    ::
    
        Data attributes:
            probeID        
            median        
            strand
            mismatchNb
        
        Functions:        
            owningStructure -- get the surrounding element (gene or exon)
            exon ------------- test if the probe is in an exon
            intron ----------- test if the probe is in an intron
            up_or_down ------- test if the probe is up or down from a gene
    """
    def __init__(self,probeID=None,regionID=None,chromosome=None,median=None,strand=None,mismatchNb=None,geneIDs=None,updownGeneIDs=None,exonIDs=None,predictedTranscriptIDs=None,updownPredictedTranscriptIDs=None,predictedExonIDs=None,localisation=None):
        self.probeID=probeID        
        self.median=median        
        self.strand=strand
        self.mismatchNb=mismatchNb
           
          
    # Test if the probe is in a specific structure (gene, exon, transcript)    
    def owningStructure(self,arrays,overlapFlag):
        # array is a list of a given structure (gene, predicted transcript, exon)
        # the current structure is described by
        # its start position, end position, strand and ID
        # get the median index of the position
        #if overlaFlag==1, overlaping of structures is allowed        
        # test if the probe is in the current structure        
        
        
                        
        greater = arrays.starts<=self.median        
        lesser = arrays.ends>=self.median        
        strandened = arrays.strands==self.strand        
        isOwned = greater & lesser & strandened                                    
        indexes=[]                       
        for i in range(len(isOwned)):
            # test if the position is in structure of index structureIDL                            
            if isOwned[i]==True:                
                # Return the structure identifier and stop the search
                # the position is unique and it is allowed to test overlapping structures        
                indexes.append(i)
                if overlapFlag==0:
                    break                        
        return indexes
           
    def exon(self,arrays):
        """Test if the probe is in an exon"""  
        newArrays=StructureList(idList=arrays.IDs,startArray=arrays.exonStarts,endArray=arrays.exonEnds,strandArray=arrays.strands)  
        indexes=self.owningStructure(newArrays,1)
        if indexes!=[]:            
            return [indexes,'e']
        else:              
            return[[],'']
                                
    def intron(self,arrays,geneStart,geneEnd):
        """Test if the probe is in an intron"""
        idList=[]
        startArray=[]
        endArray=[]
        strandArray=[]  
        if geneStart<arrays.exonStarts[0]:    
            idList=['start'];
            startArray=[geneStart]
            endArray=[arrays.exonStarts[0]]
            strandArray=[arrays.strands[0]]
        idList.extend(arrays.IDs[0:len(arrays.intronStarts)])
        startArray.extend(arrays.intronStarts)
        endArray.extend(arrays.intronEnds)
        strandArray.extend(arrays.strands[0:len(arrays.intronStarts)])
        if geneEnd>arrays.exonEnds[-1]:                
            idList.append('end')                             
            startArray.append(arrays.exonEnds[-1])            
            endArray.append(geneEnd)            
            strandArray.append(arrays.strands[0])                                        
#        if len(arrays.exonStarts)>1:        
        startArray=array(startArray)
        endArray=array(endArray)    
        strandArray=array(strandArray)
        if len(startArray)>1:
            newArrays=StructureList(idList,startArray,endArray,strandArray)            
            indexes=self.owningStructure(newArrays,1)         
            if indexes!=[]:                                
                return [indexes,'i']
            else:
                return[[],'']
        else:                        
            return[[],'']            
            
    def up_or_down(self,arrays):
        """Test if the probe is up or down from a gene"""        
        localisations=[]
        IDs=[]             
        newArrays=GeneList(idList=arrays.IDs,startArray=arrays.starts-2000,endArray=arrays.starts+1,strandArray=arrays.strands)                        
        indexes=self.owningStructure(newArrays,1)
        if indexes!=[]:            
            for index in indexes:
                IDs.append(arrays.IDs[index])                
            if self.strand==1:
                localisations=['d']*len(indexes)
            else:
                localisations=['u']*len(indexes)                
        newArrays=GeneList(idList=arrays.IDs,startArray=arrays.ends-1,endArray=arrays.ends+2000,strandArray=arrays.strands)                
        indexes=self.owningStructure(newArrays,1)        
        if indexes!=[]:            
            for index in indexes:                
                IDs.append(arrays.IDs[index])                                        
            if self.strand==1:
                localisations.extend(['u']*len(indexes))
            else:
                localisations.extend(['d']*len(indexes))                                                 
        return [IDs,localisations]
    

class Targets(object):
    """
    **Define a target.**
    ::
    
        Data attributes:
            repetition
            IDs
            indexes
            groups
            positions
            strands
            mismatchNbs
        
    """


    def __init__(self,repetition=[],IDs=[],indexes=[],groups=[],mismatchNbs=[],positions=[],strands=[]):                                     
        self.repetition=repetition
        self.IDs=IDs
        self.indexes=indexes
        self.groups=groups
        self.positions=positions
        self.strands=strands
        self.mismatchNbs=mismatchNbs
