#!/usr/bin/python
"""**Process AceView data.**"""
import bsddb, cPickle,os,time
import ps_class
from datetime import date
from numpy import array
import bintools
import setenviron
setenviron.setenviron()


def transcript_sequence(species,aceVersion,log=0):
    """
    **Read AceView files containing transcript sequences for each chromosomes.**
    ::
    
        Arguments:
            species ----- species
            aceVersion -- AceView version
            log --------- handle of a log file for recording messages
            
        Input:
            %species_region_by_chromosome.bkdb
            
        Output:
            creates %species_transcripts_sequence.bkdb
                key ==== AceView transcript id
                value == Transcript object
         
    """
    
    os.chdir(os.environ['PYDATA']+'/%s/log'%species)
    logFile=open('%s_ace_transcripts.txt'%species,'w')
    t1=time.time()
    #create ace transcript_sequence
    path=os.environ['PYDATA']+"/"+species+"/aceview/"+species+"_transcript_sequence.bkdb"
    if os.path.exists(path):
        os.remove(path)
    transcriptDB=bsddb.btopen(path,'w')
    
    #test if mRNAs sequences are in one file or in several chromosome files
    try:
        sequenceFile = open('%s/%s_%s/AceView.ncbi_37.all_mrnas_dna.fasta' %(os.environ['ACEDATA'],species,aceVersion.lower()),'r')
        chrFlag=0        
    except:      
        chrFlag=1  
    
    if chrFlag:    
        #open database for relation between chromosome and Ensembl region
        path=os.environ['PYDATA']+'/'+species+'/ensembl/'+species+'_region_by_chromosome.bkdb'
        chrDB=bsddb.btopen(path,'r')
        chromosomes=chrDB.keys()
        tscriptNb=0          
        for chromosome in chromosomes:
            print 'processing chromosome: '+chromosome
            try:     
                sequenceFile = open('%s/%s_%s/x1.all_mrnas_fasta.%s.fasta' %(os.environ['ACEDATA'],species,aceVersion.lower(),chromosome),'r')                        
                region=chrDB[chromosome]                 
                geneName='' 
                transcriptName=''
                sequence=''       
                for lines in sequenceFile:
                    tscriptNb=tscriptNb+1                
                    line = lines.split('\n')[0]
                    if not line:
                        #save last transcript
                        if geneName != '':                    
                            #save previous transcript
                            transcriptDB[transcriptName]=cPickle.dumps(ps_class.Transcript(ID=transcriptName,geneID=geneName,chromosome=chromosome,region=region,sequence=sequence),protocol=-1)
                        break
                    # get some informations 
                    if line[0]=='>':                                
                        if geneName != '':                    
                            #save previous transcript
                            transcriptDB[transcriptName]=cPickle.dumps(ps_class.Transcript(ID=transcriptName,geneID=geneName,chromosome=chromosome,region=region,sequence=sequence),protocol=-1)
                        transcriptName = line.split(':')[1]            
                        shortName=transcriptName.split(aceVersion)[0]                                                                                                                                                                                                                                         
                        transcriptLetter=shortName.split('.')[-1]
                        geneName=shortName.split('.'+transcriptLetter)[0]                
                        sequence=''                                                    
                    else:
                        # Construct sequence
                        sequence=sequence+line
            except:
                logFile.write('no AceView files %s/x1.all_mrnas_fasta.%s.fasta' %(os.environ['ACEDATA'],chromosome))                                
        transcriptDB.close()
        chrDB.close()
    else:        
        tscriptNb=0                              
        sequenceFile = open('%s/%s_%s/AceView.ncbi_37.all_mrnas_dna.fasta' %(os.environ['ACEDATA'],species,aceVersion.lower()),'r')                                                
        geneName='' 
        transcriptName=''
        sequence=''       
        for lines in sequenceFile:
            tscriptNb=tscriptNb+1                
            line = lines.split('\n')[0]
            if not line:
                #save last transcript
                if geneName != '':                    
                    #save previous transcript
                    transcriptDB[transcriptName]=cPickle.dumps(ps_class.Transcript(ID=transcriptName,geneID=geneName,sequence=sequence),protocol=-1)
                break
            # get some informations 
            if line[0]=='>':                                
                if geneName != '':                    
                    #save previous transcript
                    transcriptDB[transcriptName]=cPickle.dumps(ps_class.Transcript(ID=transcriptName,geneID=geneName,sequence=sequence),protocol=-1)
                transcriptName = line.split(':')[1]            
                shortName=transcriptName.split(aceVersion)[0]                                                                                                                                                                                                                                         
                transcriptLetter=shortName.split('.')[-1]
                geneName=shortName.split('.'+transcriptLetter)[0]                
                sequence=''                                                    
            else:
                # Construct sequence
                sequence=sequence+line                            
        transcriptDB.close()
        
    t2=time.time()
    if log!=0:
        log.write('\t%u\t%.2f\n'%(tscriptNb,t2-t1))

                    
def genes_exons(species,aceVersion,log=0):
    """
    **Process AceView gene files.**
    ::
    
        Arguments:
            species    -- species
            aceVersion -- AceView version
            log --------- handle of a log file for recording messages
        
        Input:
            %species_region_by_chromosome.bkdb
            x1.genes_gff.%chromosome.gff
        
        Output:
            Create %species_genes_by_ensembl_region.bkdb
                key ==== Ensembl region id
                value == GeneList object
            Create %species_exons_by_gene.bkdb
                key ==== AceView gene id
                value == ExonList   
            Update %species_transcript_sequence.bkdb
                key ==== AceView transcript id
                value == Transcript object       
            Create %species_transcripts_by_gene.bkdb
                key ==== AceView gene id
                value == list of AceView transcript ids       
            Create %species_transcripts_by_exon.bkdb
                key ==== AceView exon id
                value == list of AceView transcript ids       
    
    """
    
    t1=time.time()
    os.chdir(os.environ['PYDATA']+'/%s/log'%species)
    logFile=open('%s_ace_exons.txt'%species,'w')
    #open database for relation between chromosome and Ensembl region
    path=os.environ['PYDATA']+'/'+species+'/ensembl/'+species+'_region_by_chromosome.bkdb'
    chrDB=bsddb.btopen(path,'r')
    
    #open ace genes by- ensembl region
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_genes_by_ensembl_region.bkdb'
    if os.path.exists(path):
        os.remove(path)
    regionDB=bsddb.btopen(path,'w')
    regionNb=0
    
       
    #open ace exon_by_gene
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_exons_by_gene.bkdb'
    if os.path.exists(path):
        os.remove(path)
    exonsByGeneDB=bsddb.btopen(path,'w')
    exonGeneNb=0
    
    #open ace transcript_sequence for updating
    path=os.environ['PYDATA']+"/"+species+'/aceview/'+species+'_transcript_sequence.bkdb'
    transcriptDB=bsddb.btopen(path,'w')
    transcriptNb=0
    
    #create ace transcript by gene
    path=os.environ['PYDATA']+"/"+species+'/aceview/'+species+'_transcripts_by_gene.bkdb'
    if os.path.exists(path):
        os.remove(path)
    transcriptsByGeneDB=bsddb.btopen(path,'w')
    transcriptGeneNb=0
    
    #create ace transcript by exon
    path=os.environ['PYDATA']+"/"+species+'/aceview/'+species+'_transcripts_by_exon.bkdb'
    if os.path.exists(path):
        os.remove(path)
    transcriptsByExonDB=bsddb.btopen(path,'w')
    transcriptExonNb=0
    
    
    chromosomes=chrDB.keys()            
    for chromosome in chromosomes:
        try:                            
            region=chrDB[chromosome]
            print 'processing chromosome %s (Ensembl region %s)'%(chromosome,region)
            try:
                fileName = 'x1.genes_gff.%s.gff' %(chromosome)
                os.chdir('%s/%s_%s'%(os.environ['ACEDATA'],species,aceVersion.lower()))
                chrFile = open(fileName)
                chrFlag=1
            except:
                fileName = 'AceView.ncbi_37.genes_gff.gff'
                os.chdir('%s/%s_%s'%(os.environ['ACEDATA'],species,aceVersion.lower()))
                chrFile = open(fileName)
                chrFlag=0
            
            #complete list of AceView genes located in chromosomes common with Ensembl data                   
            geneIDList=[]
            geneStartList=[]
            geneEndList=[]
            geneStrandList=[]
                    
            #initiated to empty to mark the start of the process
            geneID=''
            #information recovered on each gene
            exonIDList=[]                
            exonStartList=[]
            exonEndList=[]
            exonStrandList=[]
            transcriptLetterList=[]
            read=1            
            lineRank=-1
            while read:            
                line = chrFile.readline()                    
                lineRank+=1
                if not line:
                    read=0                
                proceed=1    
                if read:                                     
                    items=line.split('\t')
                    if items[0]!=chromosome:
                        if chrFlag:                                                                                    
                            print 'abnormal line %u for chromosome %s'%(lineRank,chromosome)                        
                            print line
                            logFile.write('abnormal line %u for chromosome %s\n'%(lineRank,chromosome))
                            logFile.write('%s\n'%line)
                        proceed=0
                    else:                    
                        geneInfo = items[8].split(';')                                        
                if proceed==1: 
                    #test if the current line belong to the current gene, or if the end of file is reached       
                    if geneInfo[0].split(' ')[1]!=geneID or read==0:
                        #test if it is the first line (start of the process)        
                        if geneID=='': 
                            geneID=geneInfo[0].split(' ')[1]
                        else:                        
                            #a gene has been processed and information must be stored
                            #update list of all AceView genes
                            geneIDList.append(geneID)                           
                            geneStartList.append(min(exonStartList))
                            geneEndList.append(max(exonEndList))                        
                            geneStrandList.append(exonStrandList[0])
                            
                            # process exon information                                                                            
                            [exonIDList,exonStarts, exonEnds, exonStrands, transcriptsByExon, transcriptIDs, transcriptStarts, transcriptEnds,\
                              exonIndexes, exonGroups, intronStarts, intronEnds]=\
                              bintools.process_exons(exonStarts=array(exonStartList),exonEnds=array(exonEndList),\
                                                    exonStrands=array(exonStrandList),transcriptIDs=array(transcriptLetterList))
                            '''SAVE LIST OF TRANSCRIPTS in transcriptsByGeneDB'''                    
                            transcriptsByGeneDB[geneID]=cPickle.dumps(list(transcriptIDs),protocol=-1)
                            transcriptGeneNb=transcriptGeneNb+1
                            
                                                
                            '''SAVE LIST OF TRANSCRIPTS in transcriptsByExonDB'''
                            for i in range(len(exonStarts)):                        
                                transcriptsByExonDB[geneID+'.exon'+str(exonIndexes[i])]=cPickle.dumps(list(transcriptsByExon[i]),protocol=-1)                                
                                transcriptExonNb=transcriptExonNb+1
                                                                                                                            
                            '''SAVE LIST OF EXONS in exonByGeneDB'''                   
                            exonsByGeneDB[geneID]=cPickle.dumps(ps_class.ExonList(exonIndexes=exonIndexes,groupList=exonGroups,exonStartArray=exonStarts,exonEndArray=exonEnds,\
                                                                                 strandArray=exonStrands,intronStartArray=intronStarts,intronEndArray=intronEnds,\
                                                                                 transcriptsByExon=transcriptsByExon,transcriptIDs=transcriptIDs,\
                                                                                 transcriptStarts=transcriptStarts, transcriptEnds=transcriptEnds),protocol=-1)
                            exonGeneNb=exonGeneNb+1                            
                            '''UPDATE transcriptDB'''                                                             
                            for i in range(len(transcriptIDs)):
                                transcript=None
                                transcriptID=geneID+'.'+transcriptIDs[i]+aceVersion                          
                                try:                   
                                    transcript=cPickle.loads(transcriptDB[transcriptID])
                                except:
                                    try:
                                        transcriptID=transcriptID+'-unspliced'
                                        transcript=cPickle.loads(transcriptDB[transcriptID])
                                    except:
                                        logFile.write('transcript %s (chromosome %s) not in %s_transcripts_sequence.bkdb'%(transcriptID,chromosome,species)+'\n')        
                                if transcript is not None:                                                                
                                    transcript.start=transcriptStarts[i]
                                    transcript.end=transcriptEnds[i]
                                    transcript.strand=exonStrands[0]
                                    transcriptDB[transcriptID]=cPickle.dumps(transcript,protocol=-1)
                                    transcriptNb=transcriptNb+1                                                                                                                                                    
                                
                            #geneID of the current gene
                            geneID=geneInfo[0].split(' ')[1]
                                                    
                        exonStartList=[]
                        exonEndList=[]
                        exonStrandList=[]
                        transcriptLetterList=[]                    
                    #test if the current line contain information on exons                                                    
                    if read==1 and items[2]=='exon':                                            
                        exonStartList.append(int(items[3]))
                        exonEndList.append(int(items[4]))
                        if items[6]=='+':
                            exonStrandList.append(1)
                        else:
                            exonStrandList.append(-1)                                                    
                        transcriptName = geneInfo[2].split(' ')[2].split('.')                                 
                        transcriptLetter=transcriptName[len(transcriptName)-1].split(aceVersion)                                                                                                
                        transcriptLetterList.append(transcriptLetter[0])                                        
            '''SAVE LIST OF GENES in regionDB'''
            regionDB[region]=cPickle.dumps(ps_class.GeneList(idList=geneIDList,startArray=array(geneStartList),endArray=array(geneEndList),strandArray=array(geneStrandList)),protocol=-1)
            regionNb=regionNb+1
        except:            
            logFile.write('no AceView file %s\n'%fileName)
            pass                               
    chrDB.close()
    regionDB.close()
    exonsByGeneDB.close()
    transcriptDB.close()
    transcriptsByGeneDB.close()
    transcriptsByExonDB.close()        
    logFile.close()
    t2=time.time()
    log.write('%s\t%s\t\t16\t%s_genes_by_ensembl_region.bkdb\tace_%s\t%u\t%.2f\n'%(date.today(),species,species,aceVersion,regionNb,t2-t1))
    log.write('%s\t%s\t\t16\t%s_exons_by_gene.bkdb\tace_%s\t%u\t%.2f\n'%(date.today(),species,species,aceVersion,exonGeneNb,t2-t1))
    log.write('%s\t%s\t\t16\t%s_transcripts_by_gene.bkdb\tace_%s\t%u\t%.2f\n'%(date.today(),species,species,aceVersion,transcriptGeneNb,t2-t1))
    log.write('%s\t%s\t\t16\t%s_transcripts_by_exon.bkdb\tace_%s\t%u\t%.2f\n'%(date.today(),species,species,aceVersion,transcriptExonNb,t2-t1))
    log.write('%s\t%s\tupdate\t16\t%s_transcript_sequence.bkdb\tace_%s\t%u\t%.2f\n'%(date.today(),species,species,aceVersion,transcriptNb,t2-t1))
    
    
def ensembl_genes(species,log=0):
    """
    **Find overlaping genes between Ensembl and AceView.**
    ::
    
        Arguments:
            species -- species    
            log ------ handle of a log file for recording messages
        
        Input:
            %species_genes_by_ensembl_region.bkdb
            %species_genes_by_region.bkdb
        
        Output:
            species_ensembl_genes_by_gene.bkdb
                key ==== AceView gene id
                value == list of Ensembl genes ids
        
    """
    
    t1=time.time()
    #open ace genes by- ensembl region
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_genes_by_ensembl_region.bkdb'
    aceGeneDB=bsddb.btopen(path,'r')    
            
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+'_genes_by_region.bkdb'
    otherGenesDB=bsddb.btopen(path,'r')
    '''create ensembl genes by ace gene'''
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_ensembl_genes_by_gene.bkdb'
    otherGenesByGeneDB=bsddb.btopen(path,'w')
    geneNb=0

    for region in aceGeneDB.keys():
        if region not in otherGenesDB.keys():
            print 'region %s not in ensembl/%s_genes_by_region.bkdb'%(region,species)
        else:        
            print 'processing region ',region
            aceGenes=cPickle.loads(aceGeneDB[region])
            otherGenes=cPickle.loads(otherGenesDB[region])
            #construct dictionnary of ace genes containing ensembl genes
            ace={}
            '''process each ensembl gene in the current region'''
            for ensemblGeneIndex in range(len(otherGenes.IDs)):                    
                '''find the index of the ace genes that contains either the start or the end of the current ensembl gene'''
                aceGeneIndexes=[]
                indexes=bintools.owningStructure(otherGenes.starts[ensemblGeneIndex],otherGenes.strands[ensemblGeneIndex],aceGenes.starts,aceGenes.ends,aceGenes.strands,1)
                aceGeneIndexes.extend(indexes)            
                indexes=bintools.owningStructure(otherGenes.ends[ensemblGeneIndex],otherGenes.strands[ensemblGeneIndex],aceGenes.starts,aceGenes.ends,aceGenes.strands,1)
                aceGeneIndexes.extend(indexes)            
                '''register the current ensembl gene index in the ace dictionnary'''
                for aceGeneIndex in aceGeneIndexes:                
                    try:
                        ace[aceGenes.IDs[aceGeneIndex]].append(ensemblGeneIndex)
                    except:                    
                        ace[aceGenes.IDs[aceGeneIndex]]=[ensemblGeneIndex]                    
            
            '''process each ace gene in the current region'''    
            for aceGeneIndex in range(len(aceGenes.IDs)):                                
                '''find the index of the ensembl genes that contains either the start or the end of the current ace gene'''
                ensemblGeneIndexes=[]
                indexes=bintools.owningStructure(aceGenes.starts[aceGeneIndex],aceGenes.strands[aceGeneIndex],otherGenes.starts,otherGenes.ends,otherGenes.strands,1)            
                ensemblGeneIndexes.extend(indexes)
                indexes=bintools.owningStructure(aceGenes.ends[aceGeneIndex],aceGenes.strands[aceGeneIndex],otherGenes.starts,otherGenes.ends,otherGenes.strands,1)        
                ensemblGeneIndexes.extend(indexes)
                '''add the ensembl gene indexes found in the previous for loop processing ensembl genes'''
                try:
                    ensemblGeneIndexes.extend(ace[aceGenes.IDs[aceGeneIndex]])
                except:
                    pass            
                ensemblGeneIndexes=set(ensemblGeneIndexes)
                '''recover the corresponding ensembl gene IDs'''
                ensemblGeneIDs=[]                          
                for i in range(len(ensemblGeneIndexes)):
                    ensemblGeneIDs.append(otherGenes.IDs[ensemblGeneIndexes.pop()])                        
                '''register the current ace gene'''
                ensemblGeneIDs.sort()
                otherGenesByGeneDB[aceGenes.IDs[aceGeneIndex]]=cPickle.dumps(ensemblGeneIDs,protocol=-1)
                geneNb=geneNb+1
                                    
    otherGenesByGeneDB.close() 
    otherGenesDB.close()
    aceGeneDB.close()
    t2=time.time()
    if log!=0:
        log.write('%s\t%s\t\t17\t%s_ensembl_genes_by_gene.bkdb\taceview\t%u\t%.2f\n'%(date.today(),species,species,geneNb,t2-t1))
    
   

if __name__=='__main__':        
#    os.chdir(os.environ['PYDATA'])
#    log=open('log.txt','a')       
#    t1=time.time()
#    genes_exons('mouse','Sep07',log=log)            
#    t2=time.time()
#    print 'processed in %u s'%int(round(t2-t1))        
#    print '%s\n'%('*'*5)            
    ensembl_genes('human')