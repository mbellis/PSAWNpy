#!/usr/bin/python
"""**Write ouptut files.**"""

import os,  bsddb, time, cPickle, sys
from datetime import date
import setenviron
setenviron.setenviron()


def write_probesets_by_gene(species,chipName,log=0):
    """
    **Write several probesets_by_gene files indexed on the number of targeting probes.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            log ------- handle of a log file for recording messages
        
        Input:
            %chip_probesets_by_gene.bkdb
            %chip_probe.bkdb
            
        Output:   
            ensembl_%chip_probesets_by_gene_%probeNb.txt 
            aceview_%chip_probesets_by_gene_%probeNb.txt
        
    """    
    
    for dbType in ['ensembl','aceview']:
        #open the probesets_by_gene DBs
        ok=1
        if dbType=='ensembl':
            path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_probesets_by_gene.bkdb"                
            psDB=bsddb.btopen(path,'r')
            path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_exons_by_gene.bkdb"                
            exonDB=bsddb.btopen(path,'r')
            ensFlag=1            
        elif dbType=='aceview':
            path=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probe.bkdb"
            if os.path.exists(path):
                path=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probesets_by_gene.bkdb"                    
                psDB=bsddb.btopen(path,'r')
                path=os.environ['PYDATA']+"/"+species+"/aceview/"+species+"_exons_by_gene.bkdb"                            
                exonDB=bsddb.btopen(path,'r')
                aceFlag=1            
            else:
                ok=0
        if ok:                                
            t1=time.time()
            #create txt file
            ps={}                    
            totalProbeNb=0                
            #process each gene in the current db
            for geneID in psDB:
                #load the exons of the current gene
                lastGroup=0
                lastExon=0
                try:
                    exon=cPickle.loads(exonDB[geneID])
                    lastGroup=max(exon.groups)
                    lastExon=max(exon.indexes)
                except:
                    'pass'                                    
                #load the probesets that target the current gene        
                probesets=cPickle.loads(psDB[geneID])
                #replace transcripts names by a integer
                #construct the list of all the existing transcripts
                #use a set to have unique id
                transcr=set() 
                #process each probe nb           
                for probeNb in probesets.probesetIndexes.keys():            
                    for T in probesets.targetedTranscripts[probeNb]:
                        #PATCH : must be a list                     
                        if isinstance(T,str):                        
                            transcr.add(T)
                        else:
                            for currT in T: 
                                transcr.add(currT)                            
                    for T in probesets.notTargetedTranscripts[probeNb]:                                                                        
                        if isinstance(T,str):
                            transcr.add(T)                        
                        else:
                            for currT in T:
                                transcr.add(currT)                        
                transcripts=[]            
                for i in range(len(transcr)):
                    transcripts.append(transcr.pop())                                
                transcripts.sort()                
                
                #write info    
                try:             
                    os.mkdir(os.environ['PYDATA']+"/"+species+"/txt")
                except:
                    pass                                
                for probeNb in probesets.probesetIndexes.keys():
                    if probeNb>=0:
                        #open if necessary a new txt file and keep its handle in a dictionnary
                        if probeNb not in ps.keys():
                            totalProbeNb=max(totalProbeNb,probeNb) 
                            ps[probeNb]=open(os.environ['PYDATA']+"/"+species+"/txt/"+'%s_%s_probesets_by_gene_%02u.txt'%(dbType,chipName,probeNb),'w')
                        #gene ID                    
                        ps[probeNb].write(geneID+'\t')                    
                        #probeset ID                            
                        ps[probeNb].write('{'+repr(probesets.probesetNames[probeNb])[1:-1]+'}'+'\t')
                        #probeset indexes       
                        ps[probeNb].write(repr(probesets.probesetIndexes[probeNb])+'\t')
                        #replace transcripts names by ranks            
                        for i in range(len(probesets.targetedTranscripts[probeNb])):                        
                            currT=probesets.targetedTranscripts[probeNb][i]
                            if isinstance(currT,str):                            
                                probesets.targetedTranscripts[probeNb][i]=[transcripts.index(currT)]
                            else:                                                                          
                                for j in range(len(currT)):                        
                                    probesets.targetedTranscripts[probeNb][i][j]=transcripts.index(probesets.targetedTranscripts[probeNb][i][j])+1
                        for i in range(len(probesets.notTargetedTranscripts[probeNb])):                        
                            currT=probesets.notTargetedTranscripts[probeNb][i]
                            if isinstance(currT,str):
                                try:                                
                                    probesets.notTargetedTranscripts[probeNb][i]=[transcripts.index(currT)]
                                except:
                                    print currT
                                    print probesets.notTargetedTranscripts[probeNb][i]
                                    probesets.notTargetedTranscripts[probeNb][i][0]=transcripts.index(currT)
                                    sys.exit()
                            else:
                                for j in range(len(currT)):                                
                                    try:
                                        probesets.notTargetedTranscripts[probeNb][i][j]=transcripts.index(probesets.notTargetedTranscripts[probeNb][i][j])+1
                                    except:                                
                                        print currT
                                        print 'transcripts: '+repr(transcripts)                                
                                        print probesets.notTargetedTranscripts[probeNb][i][j]
                                        print transcripts.index(probesets.notTargetedTranscripts[probeNb][i][j])
                                        probesets.notTargetedTranscripts[probeNb][i][j]=transcripts.index(probesets.notTargetedTranscripts[probeNb][i][j])+1
                                        sys.exit()
                        #targeted exons
                        exons=probesets.targetedExons[probeNb]
                        #exon ID                        
                        ps[probeNb].write('{')                    
                        for exonIndex in range(len(exons)):
                            ps[probeNb].write('{') 
                            if exons[exonIndex]!=[]:                            
                                for exonID in exons[exonIndex].keys():
                                    ps[probeNb].write("'"+exonID+"' ")                        
                            ps[probeNb].write('} ')    
                        ps[probeNb].write('}'+'\t')                        
                        #exon rank
                        if dbType=='ensembl':
                            ps[probeNb].write('{')                    
                            for exonIndex in range(len(exons)):
                                ps[probeNb].write('[') 
                                if exons[exonIndex]!=[]:                            
                                    for exonID in exons[exonIndex].keys():
                                        try:
                                            exonPos=exon.IDs.index(exonID)
                                            ps[probeNb].write(repr(exon.indexes[exonPos])+" ")
                                        except:
                                            'pass'                                                            
                                ps[probeNb].write('] ')    
                            ps[probeNb].write('}'+'\t')                                            
                        #probe nb
                        ps[probeNb].write('{')
                        for exonIndex in range(len(exons)):
                            ps[probeNb].write('[')
                            if exons[exonIndex]!=[]:                            
                                for exonID in exons[exonIndex].keys():
                                    ps[probeNb].write(repr(exons[exonIndex][exonID])+" ")
                            ps[probeNb].write('] ')    
                        ps[probeNb].write('}'+'\t')                    
                        #last exon
                        ps[probeNb].write(repr(lastExon)+'\t')
                        #last group
                        ps[probeNb].write(repr(lastGroup)+'\t')
                        #targeted groups
                        groups=probesets.targetedGroups[probeNb]
                        #group ID
                        ps[probeNb].write('{')
                        for groupIndex in range(len(groups)):
                            ps[probeNb].write('[')
                            if groups[groupIndex]!=[]:
                                for groupID in groups[groupIndex].keys():
                                    ps[probeNb].write(repr(groupID)+" ")
                            ps[probeNb].write('] ')    
                        ps[probeNb].write('}'+'\t')                    
                        #probe nb
                        ps[probeNb].write('{')
                        for groupIndex in range(len(groups)):
                            ps[probeNb].write('[')
                            if groups[groupIndex]!=[]:
                                for groupID in groups[groupIndex].keys():
                                    ps[probeNb].write(repr(groups[groupIndex][groupID])+" ")
                            ps[probeNb].write('] ')    
                        ps[probeNb].write('}'+'\t')                    
                        
                        #targeted transcripts                                                        
                        ps[probeNb].write('{'+repr(probesets.targetedTranscripts[probeNb])[1:-1]+'}'+'\t')
                        #not targeted transcripts            
                        ps[probeNb].write('{'+repr(probesets.notTargetedTranscripts[probeNb])[1:-1]+'}'+'\t')
                        #nb of probes not in exons (up + down + intron)
                        ps[probeNb].write('{'+repr(probesets.notInExonProbeNbs[probeNb])+'}'+'\t')
                        #nb of probes in the genes (up + down + exon + splice + intron)
                        ps[probeNb].write('{'+repr(probesets.inGeneProbeNbs[probeNb])+'}'+'\n')
                        
            for probeNb in range(1,totalProbeNb+1):
                #some probeNb could miss
                try:        
                    ps[probeNb].close()
                except:
                    'pass'        
            psDB.close()
            exonDB.close()
            t2=time.time()
            if log!=0:
                if ensFlag:
                    log.write('%s\t%s\t%s\t28\t%s_%s_probesets_by_gene.txt\tensembl\t\t%.2f\n'%(date.today(),species,chipName,dbType,chipName,t2-t1))
                elif aceFlag:                                              
                    log.write('%s\t%s\t%s\t28\t%s_%s_probesets_by_gene.txt\taceviewl\t\t%.2f\n'%(date.today(),species,chipName,chipName,dbType,t2-t1))
        
    

def write_probeset_list(species,chipName,log=0):
    """
    **Write several probesets_by_gene files indexed on the number of targeting probes.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            log ------- handle of a log file for recording messages
        
        Input:
            %chip_probeset.dump
            
        Output:   
            %chip_probesets_ensembl.txt         
            %chip_probesets_aceview.txt
    
    """
    
    os.chdir(os.environ['PYDATA']+"/"+species+"/")
    inFile=open(os.environ['PYDATA']+"/"+species+"/"+chipName+'_probeset.dump','r')
    probesets=cPickle.load(inFile)
    inFile.close()    
        
    
    outFile=open(os.environ['PYDATA']+"/"+species+"/txt/"+chipName+'_probesets_ensembl.txt','w')
    psPos=-1
    t1=time.time()    
    maxProbeNb=0
    for geneNbs in probesets.ensExonGeneNbs:
        maxProbeNb=max(maxProbeNb,len(geneNbs))
    for geneNbs in probesets.ensExonGeneNbs:
        psPos+=1
        output=[probesets.probesetIndex[psPos]];  
#        print repr(output)              
#        if len(geneNbs)>probeNb:                                
#            j=[sum(geneNbs[0:len(geneNbs)-probeNb+1])]
#            j.extend(geneNbs[len(geneNbs)-probeNb+1:])
#            geneNbs=j             
#        elif len(geneNbs)<probeNb:
#            geneNbs.extend([0]*(probeNb-len(geneNbs))) 
#        elif geneNbs is None:
#            geneNbs=[0]*probeNb
        if geneNbs is None:
            geneNbs=[0]*maxProbeNb
        else:        
            geneNbs.extend([0]*(maxProbeNb-len(geneNbs)))
        output.extend(geneNbs)                                                                                                                                                                                
        outFile.write(repr(output)+'\n')        
    outFile.close()    
    t2=time.time()
    if log!=0:
        log.write('%s\t%s\t%s\t29\t%s_probesets_ensembl.txt\tensembl\t\t%.2f\n'%(date.today(),species,chipName,chipName,t2-t1))
    
    path=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probe.bkdb"
    if os.path.exists(path):
        outFile=open(os.environ['PYDATA']+"/"+species+"/txt/"+chipName+'_probesets_aceview.txt','w')
        psPos=-1
        t1=time.time()
        maxProbeNb=0
        for geneNbs in probesets.aceExonGeneNbs:
            maxProbeNb=max(maxProbeNb,len(geneNbs))
        for geneNbs in probesets.aceExonGeneNbs:
            psPos+=1
            output=[probesets.probesetIndex[psPos]];        
#            if len(geneNbs)>probeNb:            
#                j=[sum(geneNbs[0:len(geneNbs)-probeNb+1])]
#                j.extend(geneNbs[len(geneNbs)-probeNb+1:])
#                geneNbs=j                                
#            elif len(geneNbs)<probeNb:
#                geneNbs.extend([0]*(probeNb-len(geneNbs)))
#            elif geneNbs is None:
#                geneNbs=[0]*probeNb
            if geneNbs is None:
                geneNbs=[0]*maxProbeNb
            else:        
                geneNbs.extend([0]*(maxProbeNb-len(geneNbs)))
            output.extend(geneNbs)    
            outFile.write(repr(output)+'\n')
        outFile.close()
        t2=time.time()
        if log!=0:
            log.write('%s\t%s\t%s\t29\t%s_probesets_aceview.txt\taceview\t\t%.2f\n'%(date.today(),species,chipName,chipName,t2-t1))
            
def write_ensembl_genes(species,log=0): 
    """
    **Write Ensembl gene list.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            log ------- handle of a log file for recording messages
        
        Input:
            %species_genes_by_region.bkdb
            %species_ensembl_genes_by_gene.bkdb
            
        Output:           
            %species_ensembl__genes.txt
        
    """
           
    t1=time.time()


    path=os.environ['PYDATA']+'/'+species+'/ensembl/'+species+'_genes_by_region.bkdb'
    byRegionDB=bsddb.btopen(path)    
    outFile=open(os.environ['PYDATA']+"/"+species+'/txt/'+species+'_ensembl_genes.txt','w')
        
                    
    ensGeneIDs=[]
    ensGeneNb=0
    for region in byRegionDB:
        ensGeneIDs=cPickle.loads(byRegionDB[region]).IDs
        ensGeneNb+=len(ensGeneIDs)
        for ensGeneID in ensGeneIDs:
            outFile.write('%s\n'%ensGeneID)
    outFile.close()        
                
    byRegionDB.close()       


    t2=time.time()
    if log!=0:
        log.write('%s\t%s\t\t30\t%s_ensembl_genes.txt\taceview\t%u\t%.2f\n'%(date.today(),species,species,ensGeneNb,t2-t1))            

def write_aceview_genes(species,log=0): 
    """
    **Write correspondance between AceView and Ensembl genes.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            log ------- handle of a log file for recording messages
        
        Input:
            %species_genes_by_region.bkdb
            %species_ensembl_genes_by_gene.bkdb
            
        Output:           
            %species_ens_by_ace_gene.txt
        
    """
           
#open ace genes by- ensembl region
    t1=time.time()


    path=os.environ['PYDATA']+'/'+species+'/ensembl/'+species+'_genes_by_region.bkdb'
    byRegionDB=bsddb.btopen(path)    
    path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_ensembl_genes_by_gene.bkdb'
    byAceGeneDB=bsddb.btopen(path)
    outFile=open(os.environ['PYDATA']+"/"+species+'/txt/'+species+'_ens_by_ace_gene.txt','w')
        
                    
    byAceGeneIDs=byAceGeneDB.keys()
    ensGeneIDs=[]
    for region in byRegionDB:
        ensGeneIDs.extend(cPickle.loads(byRegionDB[region]).IDs)        
    for byAceGeneID in byAceGeneIDs:
        ensGenes=cPickle.loads(byAceGeneDB[byAceGeneID])                        
        outFile.write(byAceGeneID+'\t')
        if ensGenes==[] or ensGenes==None:
            outFile.write('{''}\n')
        else:
            outFile.write('{')
            for i in range(len(ensGenes)):
                try:
                    #current ensemnl gene may have been already erased
                    ensGeneIDs.pop(ensGeneIDs.index(ensGenes[i]))
                except:
                    pass
                if i==0:
                    outFile.write("'"+ensGenes[i]+"'")
                else:
                    outFile.write(','+"'"+ensGenes[i]+"'") 
            outFile.write('}\n')
    for i in range(len(ensGeneIDs)):
        outFile.write('0\t%s\n'%ensGeneIDs[i])
    byRegionDB.close()       
    byAceGeneDB.close()
    outFile.close()
    t2=time.time()
    if log!=0:
        log.write('%s\t%s\t\t30\t%s_aceview_genes.txt\taceview\t%u\t%.2f\n'%(date.today(),species,species,len(byAceGeneIDs),t2-t1))
        
def write_ps_length(species,chipName):
    """Write the genomic distance between most 5' and most 3' probes
    for all probesets
    """
    outFile=open(os.environ['PYDATA']+'/'+chipName+'_probe_length.txt','w')
    
    path=os.environ['COMMONDB']+"/chip_probeset.bkdb"        
    chipDB=bsddb.btopen(path,'r')
    
    path=os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.bkdb"            
    psDB=bsddb.btopen(path)
    
    probesets=cPickle.loads(chipDB[chipName])   
    for probeset in probesets:
        ps=cPickle.loads(psDB[probeset])
        outFile.write('%s\t%u\n'%(probeset,ps.probesetTargetLength))
    outFile.close()    
                        
        
    
            
        
        
if __name__=='__main__':
    
    #write_aceview_genes('mouse')    
    write_ps_length('rat','m6')