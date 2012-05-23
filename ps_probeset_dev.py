#!/usr/bin/python
"""**Functions for assigning probe set to gene(s) according to the number of targeting probes.**"""

import bsddb, cPickle, os, time, sys
from bsddb import import db
import copy as copypy
import ps_class
from numpy import *
import seqcomp
import random
from datetime import date
import setenviron
setenviron.setenviron()

def position_mapping(species,chipName,dbType,region,log=0):
    """
    **Map positions to gene structures (exon, intron, upstream or downstream).**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            dbType ---- either 'ensemble' of 'aceview'
            region ---- Ensembl region id
            log -------- handle of a log file for recording messages
            
        Input:        
            %species_genes_by_region.bkdb or
            %species_genes_by_ensembl_region.bkdb (Ensembl or AceView database according to dbType argument)                
            %chip_positions_by_region.bkdb
            
        Output:    
            update %species_exons_by_gene.bkdb (Ensembl or AceView database according to dbType argument)
            update %chip_probe.bkdb (Ensembl or AceView database according to dbType argument)
            
        Log files:
            %chip_ensembl_position_mapping_log.txt or
            %chip_aceview_position_mapping_log.txt (according to dbType argument)    
        
    """
    
    t1=time.time()     
          
    if dbType=='ensembl':        
        log1=open(os.environ['PYDATA']+"/"+species+'/log/'+'%s_ensembl_position_mapping_log.txt'%chipName,'w')
        
        path=os.environ['PYDATA']+"/"+species+'/'+dbType            
        dbEnv = db.DBEnv()        
        dbEnv.open(path,db.DB_INIT_LOCK|db.DB_CREATE|db.DB_INIT_MPOOL)
        
        
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+'_genes_by_region.bkdb'        
        log1.write(path+'\n')
        genesByRegionDB=db.DB(None)
        genesByRegionDB=db.open(path)
        genesByRegionDB=bsddb.btopen(path,None,db.DB_RDONLY)        
        
        # uses only gene because probe position refers only to gene regions, not to predicted transcript region  
        # open database containing the exon positions grouped by gene
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_exons_by_gene.bkdb"
        #print path
        log1.write(path+'\n')
        exonsByGeneDB=db.DB(dbEnv)
        exonsByGeneDB=db.open(path,None,db.DB_BTREE,db.DB_CREATE)
                                    
        # open probe database to complete records
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+chipName+"_probe.bkdb"
        #print path
        log1.write(path+'\n')
        probeDB_name= chipName+'_probe.bkdb'
        probeDB=db.DB(dbEnv)
        probeDB=db.open(path,None,db.DB_BTREE,db.DB_CREATE)
        
    elif dbType=='aceview':        
                   
        log1=open(os.environ['PYDATA']+"/"+species+'/log/'+'%s_aceview_position_mapping_log.txt'%chipName,'w')
        
        path=os.environ['PYDATA']+"/"+species+'/'+dbType            
        dbEnv = db.DBEnv()        
        dbEnv.open(path,db.DB_INIT_LOCK|db.DB_CREATE|db.DB_INIT_MPOOL)    
        
        # open database containing the genes grouped by region
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_genes_by_ensembl_region.bkdb"
        print path
        log1.write(path+'\n')
        genesByRegionDB=db.DB(None)
        genesByRegionDB=db.open(path)
        genesByRegionDB=bsddb.btopen(path,None,db.DB_RDONLY)
        
        # uses only gene because probe position refers only to gene regions, not to predicted transcript region  
        # open database containing the exon positions grouped by gene
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_exons_by_gene.bkdb"
        print path
        log1.write(path+'\n')
        exonsByGeneDB=db.DB(dbEnv)
        exonsByGeneDB=db.open(path,None,db.DB_BTREE,db.DB_CREATE)
        
        # open probe database to complete records
        path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+chipName+"_probe.bkdb"
        print path
        log1.write(path+'\n')
        probeDB=db.DB(dbEnv)
        probeDB=db.open(path,None,db.DB_BTREE,db.DB_CREATE)
    
    dbEnv.set_lk_detect(db.DB_LOCK_YOUNGEST)
    exonsLockID = dbEnv.lock_id()
    probeLockID = dbEnv.lock_id()    
    # open database containing the probe positions grouped by region
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_positions_by_region.bkdb"    
    log1.write(path+'\n')
    positionsByRegionDB=bsddb.btopen(path)    

    t2=time.time()                            
    log1.write('\n'+'processing region %s'%region+'\n')            
    
    try:
        # get the list of positions in the current region
        positionList = cPickle.loads(positionsByRegionDB[region])                                                            
        # get the list of genes found in the current region
    except:
        positionList = None
        print 'no region %s in positionsBYRegionDB'%region
    try:        
        geneArrays=cPickle.loads(genesByRegionDB[region])   
        print 'processing region %s'%region                            
    except:       
        print 'no region %s for genesByRegionDB ',region        
        geneArrays=None                 
                                       
    if positionList is not None:
        for posL in range(len(positionList.medians)):
            position=ps_class.Position(probeID=positionList.probeIDs[posL],median=positionList.medians[posL],strand=positionList.strands[posL],mismatchNb=positionList.mismatches[posL])
            # test if there are genes in the current region
            foundPosition=0
            useCurrentProbe=0
            if geneArrays != None:
                # search if the current position is inside gene(s)
                geneIDs=[]
                geneStarts=[]
                geneEnds=[]                           
                geneIndexes=position.owningStructure(geneArrays,1)
                if geneIndexes!=[]:                        
                    for geneIndex in geneIndexes:
                        geneIDs.append(geneArrays.IDs[geneIndex])
                        geneStarts.append(geneArrays.starts[geneIndex])
                        geneEnds.append(geneArrays.ends[geneIndex])                                                                                                                                        
                    foundPosition=1
                    #load the probe corresponding to the current position
                    foundProbe=0
                    for probeID in position.probeID:
                        try:
                            probe=cPickle.loads(probeDB[probeID])                       
                            useCurrentProbe=1
                            foundProbe=1
                            break
                        except:
                            pass       
                    if foundProbe==0:
#                            print 'Genes: %s probe not in %s'%(position.probeID,probeDB_name)  
                        log1.write('Genes step : %s probe not in %s'%(position.probeID,probeDB_name)+'\n')                      
                    if useCurrentProbe==1:                        
                        #process all the genes that contains the current position                              
                        for geneIndex in range(len(geneIDs)):
                            geneID=geneIDs[geneIndex]                                
                            try:    
                                while(1):
                                    try:            
                                        lock = dbEnv.lock_get(exonsLockID, "exonid", db.DB_LOCK_WRITE,db.DB_LOCK_NOWAIT)            
                                        break
                                    except Exception,error:            
                                        pass             
                                exons=cPickle.loads(exonsByGeneDB.get[geneID])
                                processExons=1
                                if exons.IDs==[]:                                   
                                    for i in range(len(exons.indexes)):
                                        exons.IDs.insert(i,geneID+'.exon%u'%exons.indexes[i])
                                        
                                    exonsByGeneDB.put(geneID,cPickle.dumps(exons,protocol=-1))
                                    dbEnv.lock_put(lock)                                                                                
                            except:
                                processExons=0
#                                    print 'no entry for %s in exons_by_gene'%geneID
                                log1.write('no entry for %s in exons_by_gene'%geneID+'\n')
                            if processExons==1:                                
                                # test it the position is inside one or several exons
#                                    print exons.exonStarts,'\t',exons.exonEnds,'\t',exons.strands
#                                    print type(exons.exonStarts),'\t',type(exons.exonEnds),'\t',type(exons.strands)
                                                                                                                                                                    
                                [indexes,localisation]=position.exon(exons)                                
#                                    print 'exon: gene %s indexes %s localisation %s'%(geneID,repr(indexes),localisation)
                                if localisation=='e':                                                                                                                                
                                    try:    
                                        targets=probe.genes4exon[geneID]                                                                            
                                        for exonIndex in range(len(indexes)):                                                                                                                                                                                                                
                                                targets.repetition.append(1)
                                                targets.indexes.append(exons.indexes[indexes[exonIndex]])
                                                targets.groups.append(exons.groups[indexes[exonIndex]])                                            
                                                targets.mismatchNbs.append(position.mismatchNb)  
                                                targets.positions.append(position.median)
                                                targets.strands.append(position.strand)
                                                targets.IDs.append(exons.IDs[indexes[exonIndex]])                                 
                                    except:                                        
                                        targets=ps_class.Targets(repetition=[1],IDs=[exons.IDs[indexes[0]]],indexes=[exons.indexes[indexes[0]]],groups=[exons.groups[indexes[0]]],mismatchNbs=[position.mismatchNb],positions=[position.median],strands=[position.strand])                                        
                                        if len(indexes)>0:
                                            for exonIndex in range(1,len(indexes)):                                                                                        
                                                targets.repetition.append(1)
                                                targets.indexes.append(exons.indexes[indexes[exonIndex]])
                                                targets.groups.append(exons.groups[indexes[exonIndex]])                                            
                                                targets.mismatchNbs.append(position.mismatchNb)  
                                                targets.positions.append(position.median)
                                                targets.strands.append(position.strand)
                                                targets.IDs.append(exons.IDs[indexes[exonIndex]])                                                                                                             
                                    probe.genes4exon[geneID]=targets
                                # test it the position is inside one or several introns
#                                    print geneStarts[geneIndex],'\t',geneEnds[geneIndex]                                                                                
                                [indexes,localisation]=position.intron(exons,geneStarts[geneIndex],geneEnds[geneIndex])                                                                        
#                                    print 'intron: gene %s indexes %s localisation %s'%(geneID,repr(indexes),localisation)
                                if localisation=='i':                                                                                                                                        
                                    try:    
                                        targets=probe.genes4intron[geneID]
                                        for intronIndex in range(len(indexes)):
                                            targets.repetition.append(1)
                                            # introns are constructed from grouped exons
                                            # for example, if there are 20 exons, 5 of wich are overlapping
                                            # 15 groups are constructed giving 14 introns
                                            if position.strand==1:
                                                targets.indexes.append(indexes[intronIndex])
                                            else:
                                                targets.indexes.append(max(exons.groups)-indexes[intronIndex]-1)                                                                                                                                
                                            targets.mismatchNbs.append(position.mismatchNb)
                                            targets.positions.append(position.median)
                                            targets.strands.append(position.strand)                                   
                                    except:
                                        if position.strand==1:
                                            intron=indexes[0]
                                        else:
                                            intron=max(exons.groups)-indexes[0]-1
                                        targets=ps_class.Targets(repetition=[1],indexes=[intron],mismatchNbs=[position.mismatchNb],positions=[position.median],strands=[position.strand])
                                        if len(indexes)>0:
                                            for intronIndex in range(1,len(indexes)):
                                                targets.repetition.append(1)
                                                # introns are constructed from grouped exons
                                                # for example, if there are 20 exons, 5 of wich are overlapping
                                                # 15 groups are constructed giving 14 introns
                                                if position.strand==1:
                                                    targets.indexes.append(indexes[intronIndex])
                                                else:
                                                    targets.indexes.append(max(exons.groups)-indexes[intronIndex]-1)                                                                                                                                
                                                    targets.mismatchNbs.append(position.mismatchNb)
                                                    targets.positions.append(position.median)
                                                    targets.strands.append(position.strand)                                                                            
                                    probe.genes4intron[geneID]=targets                                    
                        
                #search for up or down positions   
                             
                [updownGeneIDs,updownGeneLocalisations]=position.up_or_down(geneArrays)  
                            
#                    print 'updown: genes %s localisation %s'%(updownGeneIDs,updownGeneLocalisations)
                if updownGeneIDs is not None:
                    foundPosition=1                    
                    if useCurrentProbe==0:
                        for probeID in position.probeID:
                            try:
                                probe=cPickle.loads(probeDB[probeID])                       
                                useCurrentProbe=1
                                foundProbe=1
                                break
                            except:
                                pass       
                        if foundProbe==0:                        
#                                print 'updown genes : %s probe not in %s'%(position.probeID,probeDB_name)    
                            log1.write('updown GENES : %s probe not in %s'%(position.probeID,probeDB_name)+'\n')
                    if useCurrentProbe==1:                                                                                  
                        for updownIndex in range(len(updownGeneIDs)):                            
                            if updownGeneLocalisations[updownIndex]=='u':                                                     
                                try:    
                                    targets=probe.genes4up[updownGeneIDs[updownIndex]]
                                    targets.repetition+=1                                    
                                    targets.mismatchNbs.append(position.mismatchNb)
                                    targets.positions.append(position.median)
                                    targets.strands.append(position.strand)                                   
                                except:
                                    targets=ps_class.Targets(repetition=[1],mismatchNbs=[position.mismatchNb],positions=[position.median],strands=[position.strand])                                                                                                                                    
                                probe.genes4up[updownGeneIDs[updownIndex]]=targets                                                                                                                    
                            elif updownGeneLocalisations[updownIndex]=='d':
                                try:    
                                    targets=probe.genes4down[updownGeneIDs[updownIndex]]
                                    targets.repetition+=1                                    
                                    targets.mismatchNbs.append(position.mismatchNb)
                                    targets.positions.append(position.median)
                                    targets.strands.append(position.strand)                                   
                                except:
                                    targets=ps_class.Targets(repetition=[1],mismatchNbs=[position.mismatchNb],positions=[position.median],strands=[position.strand])                            
                                probe.genes4down[updownGeneIDs[updownIndex]]=targets                                
                                                                                                                                                                                                                                                                
                if useCurrentProbe==1:
                    # update the probe data base  
                    while(1):
                        try:            
                            lock = dbEnv.lock_get(probeLockID, "probeid", db.DB_LOCK_WRITE,db.DB_LOCK_NOWAIT)            
                            break
                        except Exception,error:            
                            pass    
                    probeDB.put(position.probeID,cPickle.dumps(probe,protocol=-1))
                    dbEnv.lock_put(lock)                                                                                                                                                 
                                            
                                                     
            if foundPosition==0:
                for probeID in position.probeID:
                    try:
                        probe=cPickle.loads(probeDB[probeID])
                        probe.outOfGeneNb+=1
                        while(1):
                            try:            
                                lock = dbEnv.lock_get(probeLockID, "probeid", db.DB_LOCK_WRITE,db.DB_LOCK_NOWAIT)            
                                break
                            except Exception,error:            
                                pass    
                        probeDB.put(position.probeID,cPickle.dumps(probe,protocol=-1))
                        dbEnv.lock_put(lock)                                                                                                                 
                        foundProbe=1
                        break
                    except:
                        pass       
                if foundProbe==0:                                     
#                        print 'outOfGenes : %s probe not in %s'%(position.probeID,probeDB_name)                                                        
                    log1.write('outOfGenes : %s probe not in %s'%(position.probeID,probeDB_name))
    t3=time.time()    
    print 'region processed in %u s'%int(round(t3-t2))                
#                break        
#            break
                                    
    positionsByRegionDB.close()
    genesByRegionDB.close()
    exonsByGeneDB.close()
    probeDB.close()
    log1.close()
    t4=time.time()
    if log!=0:
        if dbType=='ensembl':
            log.write('%s\t%s\t\t19\t%s_exons_by_gene.bkdb\t\tensembl update\t%.2f\n'%(date.today(),species,species,t4-t1))
            log.write('%s\t%s\t%s\t19\t%s_probe.bkdb\t\tensembl update\t%.2f\n'%(date.today(),species,chipName,chipName,t4-t1))
        else:
            log.write('%s\t%s\t\t19\t%s_exons_by_gene.bkdb\t\taceview update\t%.2f\n'%(date.today(),species,species,t4-t1))
            log.write('%s\t%s\t%s\t19\t%s_probe.bkdb\t\taceview update\t%.2f\n'%(date.today(),species,chipName,chipName,t4-t1))
            
            
#    print 'process terminated'

def erase_gop(species,chipName):
    # open probe database to update records
    path=os.environ['PYDATA']+"/"+species+'/ensembl/'+chipName+"_probe.bkdb"
    print path        
    probeDB=bsddb.btopen(path)
    #erase target in GOP
    erasedNb=0
    for probeID in probeDB.keys():
        try:    
            currProbe=cPickle.loads(probeDB[probeID])
        except:
            pass        
        saveProbe=0
        for geneID in currProbe.genes4exon.keys():                        
            if geneID[0:3]=='GOP':             
                print geneID   
                saveProbe=1
                currProbe.genes4exon.pop(geneID)
                erasedNb+=1
        if saveProbe:
            probeDB[probeID]=cPickle.dumps(currProbe,protocol=-1)
    print '%u GOP erased'%erasedNb                
    probeDB.close()

def position_gopmapping(species,chipName,region,log=0):
    """
    **Map positions to group of probes (GOPs).**    
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            dbType ---- either 'ensemble' of 'aceview'
            region ---- Ensembl region id
            log -------- handle of a log file for recording messages
            
        Input:        
            %chip_gops__by_region.bkdb                    
            %chip_positions_by_region.bkdb
            
        Output:            
            update %chip_probe.bkdb (Ensembl or AceView database according to dbType argument)
            
        Log files:
            %species_ensembl_position_gopmapping_log.txt            
            
    """    
        
    #open log file specific for position_gopmapping    
    gop_log=open(os.environ['PYDATA']+"/"+species+'/log/'+species+'_ensembl_position_gopmapping_log.txt','w')
    gop_log.write('position_gopmapping\n')
    
    # open database containing the GOP by region for the current species
    path=os.environ['PYDATA']+"/"+species+'/ensembl/'+chipName+"_gops_by_region.bkdb"
    print path
    gop_log.write(path+'\n')
    gopDB=bsddb.btopen(path,'r')
    
    # open database containing the probe positions grouped by region for the current chip    
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_positions_by_region.bkdb"
    print path
    gop_log.write(path+'\n')
    positionDB=bsddb.btopen(path,'r')
        
    # open probe database to update records
    path=os.environ['PYDATA']+"/"+species+'/ensembl/'+chipName+"_probe.bkdb"
    print path
    gop_log.write(path+'\n')
    probeDB_name= chipName+"_probe.bkdb"
    probeDB=bsddb.btopen(path)        
    print '\n'        
        
        
    t1=time.time()            
    gop_log.write('\n'+'processing region %s'%region+'\n')            
    
    # get the list of positions in the current region of the current chip
    positionList = cPickle.loads(positionDB[region])                                                   
    # get the list of GOP found in the current region
    try:        
        gopArrays=cPickle.loads(gopDB[region])
        print 'processing region %s'%region                               
    except:       
        gopArrays=None                     
        print 'no region %s in gopDB'%region    
        
        
     
    if gopArrays is not None:
        #process each position of the current region of the current chip                                                            
        for posL in range(len(positionList.medians)):
            position=ps_class.Position(median=positionList.medians[posL],strand=positionList.strands[posL],mismatchNb=positionList.mismatches[posL])                                                                                                                                                                            
            useCurrentProbe=0                
            # search if the current position is contained inside one or several GOPs                                           
            gopIndexes=position.owningStructure(gopArrays,1)                
            if gopIndexes!=[]:                                
                #load the probe corresponding to the current position
                try: 
                    probe=cPickle.loads(probeDB[position.probeID])
                    useCurrentProbe=1
                    #process all the gop that contains the current position
                    for gopIndex in gopIndexes:
                        gopID =gopArrays.IDs[gopIndex]                      
                        #the same probe can be repeated inside a GOP                                                                                                                                                                                                                                                                   
                        try:                    
                            targets=probe.genes4exon[gopID]                                                                                                                                                                                                                                                                                                                
                            targets.repetition.append(1)
                            targets.indexes.append(0)
                            targets.groups.append(0)                                            
                            targets.mismatchNbs.append(position.mismatchNb)  
                            targets.positions.append(position.median)
                            targets.strands.append(position.strand)
                            targets.IDs.append('gop')                                                        
                        except:                                                           
                            targets=ps_class.Targets(repetition=[1],IDs=['gop'],indexes=[0],groups=[0],mismatchNbs=[position.mismatchNb],positions=[position.median],strands=[position.strand])                                                                                                                                                                                                                                                                                                                                         
                        probe.genes4exon[gopID]=targets
                                             
                except:
                    print '%s probe not in %s'%(position.probeID,probeDB_name)
                    gop_log.write('%s probe not in %s'%(position.probeID,probeDB_name)+'\n')
                    gop_log.close()
                    sys.exit()                                                                                                                                                                                                        
            if useCurrentProbe==1:#                                                                                                                               
                probeDB[position.probeID]=cPickle.dumps(probe,protocol=-1)                                                                
                                                                                                                                                                                                                                                            
    t2=time.time()    
    print 'region processed in %u s'%int(round(t2-t1))                
                                
    positionDB.close()
    gopDB.close()    
    probeDB.close()
    gop_log.close()
    print 'process terminated'
    

def assign_probeset(species,chipName,dbType,psRange=0,aceVersion='',log=log):
    """
    **Assign probe set to gene(s) or GOP(s) according to the number of targeting probes.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            dbType ---- either 'ensembl' or 'aceview'
            psRankge -- range of probe set to be processed
            aceVersion -- AceView version
            log -------- handle of a log file for recording messages
            
        Inputs:        
            %species_transcripts_by_gene.bkdb
            %species_transcript_sequences.bkdb
            %species_transcripts_by_exon.bkdb
            %species_ensembl_genes_by_gene.bkdb (if dbType is aceview)
            
        
        Outputs:
            update %chip_probeset.bkdb
            update %chip_probe.bkdb    
                    
        Log files:
            %chip_ensembl_assign_probeset_log.txt or
            %chip_aceview_assign_probeset_log.txt (according to dbTYpe)
                 
    """            
        
    logFile=open(os.environ['PYDATA']+"/"+species+'/log/'+'%s_%s_assign_probeset_log.txt'%(chipName,dbType),'w')
                     
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+chipName+"_probe.bkdb"
    print path
    probeDB=bsddb.btopen(path)
    
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_transcripts_by_gene.bkdb"
    print path
    transcriptsByGeneDB=bsddb.btopen(path,'r')
    
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_transcript_sequence.bkdb"
    print path
    transcriptSequencesDB=bsddb.btopen(path,'r')
    
    path=os.environ['PYDATA']+"/"+species+'/'+dbType+'/'+species+"_transcripts_by_exon.bkdb"
    print path
    transcriptsByExonDB=bsddb.btopen(path,'r')
    
    if dbType=='aceview':                
        path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_ensembl_genes_by_gene.bkdb'
        print path
        ensByGeneDB=bsddb.btopen(path,'r')        
                                
        
    # open the chip probeset data base 
    path=os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.bkdb"
    print path
    probesetDB=bsddb.btopen(path)
    probesetIDs=probesetDB.keys()
    if psRange==0:
        psRange=range(len(probesetIDs))   
    #recover maximum probe nb
    maxProbeNb=0
    for psRank in psRange:    
        probesetID=probesetIDs[psRank]                 
        try:        
            probeset=cPickle.loads(probesetDB[probesetID])
            probeIDs=probeset.probeIDs
            if probeIDs is not None:
                try:                
                    probeNb=len(probeIDs)
                    maxProbeNb=max(maxProbeNb,probeNb)
                except:
                    pass
        except:
            pass
    for psRank in psRange:    
        probesetID=probesetIDs[psRank]            
        # process each probeset        
        try:        
            probeset=cPickle.loads(probesetDB[probesetID])
        except:
            print "CPICKLE ERROR FOR PROBESET %s (rank %u)%(probesetID,psRank)"
            if log!=0:
                log.write('%s\t%s\t\t24\t%s_probeset.bkdb\t\tensembl update\t"CPICKLE ERROR FOR PROBESET %s (rank %u)\n'%(date.today(),species,species,probesetID,psRank))             
            sys.exit()                
        if dbType=='ensembl':
            probeset.ensExonGeneNbs=0
            probeset.ensIntronGeneNbs=0
            probeset.ensUpGeneNbs=0
            probeset.ensDownGeneNbs=0
            probeset.ensOutProbeNbs=0
            probeset.ensNisProbeNbs=0
            probeset.ourEnsGeneIDs=[]
            probeset.ourEnsGeneNb=0
            probeset.ourEnsProbeNb=0                
            probeset.ensGeneIDs4Ens=[]
            probeset.ourEnsGeneIDs4Ens=[]
            probeset.commonGeneIDs4Ens=[]
            probeset.ensGenesByProbeNb={}                
        elif dbType=='aceview':
            probeset.aceExonGeneNbs=0
            probeset.aceUpGeneNbs=0
            probeset.aceDownGeneNbs=0
            probeset.aceOutProbeNbs=0
            probeset.aceNisProbeNbs=0            
            probeset.ourAceGeneIDs=[]
            probeset.ourAceGeneNb=0
            probeset.ourAceProbeNb=0
            probeset.ourAceToEnsGeneNbs=[]
            probeset.ourAceToEnsGeneIDs=[]
            probeset.ensGeneIDs4Ace=[]
            probeset.ourEnsGeneIDs4Ace=[]
            probeset.commonGeneIDs4Ace=[]
            probeset.aceGenesByProbeNb={}
        #load all the corresponding probes
        probeIDs=probeset.probeIDs                
        if probeIDs is not None:
            try:                
                probeNb=len(probeIDs)        
            except:
                logFile.write('no probeIDs for probesetID: %s psRank: %u'%(probesetID,psRank))
                logFile.close()
                sys.exit()
                
            '''FIND GENES WITH PROBE(S) IN THEIR EXONS OR IN THEIR SPLICE JUNCTIONS
            INDICATE FOR EACH PROBE IF IT IS IN AN EXON OR IF IT IS NOT IN SEQUENCE (NIS)'''    
            #dictionary of probes     
            probes={}
            #list of targeted genes
            genes=[]        
            #gene IDs of targeted genes
            geneIDs=[]  
            #index of probes that are potentially not in sequence => verify if they are in a transcript sequence (overlapping two exons)  
            nisIndex=[0]*probeNb
            #index of probes that are at least in one exon 
            inExonIndex=[0]*probeNb
            #for each gene indicate which probe is in one of its exons
            inExonIndexes=[]        
            probeIndex=-1
                        
            # process each probe of the current probeset
            for probeIndex in range(probeNb):
                probeID=probeIDs[probeIndex]                                        
                probe=cPickle.loads(probeDB[probeID])                                      
                probes[probeID]=probe
                # recover index of probe with at least one position in an exon
                # and assemble the list of genes with a probe in one of their exon                                        
                if probe.genes4exon!={}:
                    inExonIndex[probeIndex]=1                    
                    for geneID in probe.genes4exon.keys():                                                                                                                    
                        if geneID not in geneIDs:                                            
                            # CREATE a TargetedGene INSTANCE
                            geneIDs.append(geneID)   
                            genes.append(ps_class.TargetedGene())
                            currIndex=[0]*probeNb
                            currIndex[probeIndex]=1
                            inExonIndexes.append(currIndex)
                        else:
                            # UPDATE the existing TargetedGene INSTANCE
                            geneIndex=geneIDs.index(geneID)
                            inExonIndexes[geneIndex][probeIndex]=1
                                                                                                                        
                # recover index of probes that are potentially not in sequence (nis)    
                if probe.genes4exon=={} and  probe.genes4intron=={} and  probe.genes4up=={} and  probe.genes4down=={} and probe.outOfGeneNb==0:
                    nisIndex[probeIndex]=1
                        
            # search potentially nis probes in existing transcripts
            # if no probe are in an exon, it is assumed that it is not worth doing the search 
            # on all the existing transcripts
            # therefore the search is restricted to the genes that are already targeted by at least one probe of
            # the current probeset         
            if geneIDs !=[]:                                  
                #recover transcript sequences corresponding to genes of the current probeset
                transcripts={}
                for geneID in geneIDs:               
                    geneIndex=geneIDs.index(geneID)                                                             
                    transcriptIDs=[]
                    if dbType=='ensembl':
                        try:               
                            transcriptIDs=cPickle.loads(transcriptsByGeneDB[geneID]).IDs
                        except:
                            pass
                    elif dbType=='aceview': 
                        try:                               
                            transcriptIDs=cPickle.loads(transcriptsByGeneDB[geneID])
                            for i in range(len(transcriptIDs)):
                                transcriptIDs[i]=geneID+'.'+transcriptIDs[i]+aceVersion
                        except:
                            pass
                    if transcriptIDs !=[]:                       
                        #recover transcript sequences   
                        for transcriptID in transcriptIDs:
                            transcript=None
                            try:                            
                                transcript=cPickle.loads(transcriptSequencesDB[transcriptID])
                            except:
                                if dbType=='aceview':
                                    try:
                                        transcriptID=transcriptID+'-unspliced'
                                        transcript=cPickle.loads(transcriptSequencesDB[transcriptID])
                                    except:                                            
                                        logFile.write('transcript %s of gene %s is not in cDNA file\n'%(transcriptID,geneID))
                                else:     
                                    logFile.write('transcript %s of gene %s is not in cDNA file\n'%(transcriptID,geneID))                       
                            if transcript is not None:        
                                sequence=transcript.sequence                                      
                                # search if not assigned probes match on of the current transcript  
                                # search if made for all the probe that were not assigned to an exon in the previous step
                                # in case that the probe targets a particular transcript, its index is set to 2 to 
                                # allow the test of other transcripts of the same gene                                                                         
                                for probeIndex in range(probeNb):                                    
                                    if inExonIndexes[geneIndex][probeIndex]==0 or inExonIndexes[geneIndex][probeIndex]==2:                                                                             
                                        [mismatchNbs,positions]= seqcomp.match_with_miss(pattern=probes[probeIDs[probeIndex]].sequence,seq=sequence,mismatchNb=1)                                                                                                                                                                                                                                                    
                                        if positions !=[]:                                       
                                            try:
                                                probes[probeIDs[probeIndex]].genes4splice[geneID].append(ps_class.Targets(repetition=[1]*len(positions),indexes=[-1]*len(positions),positions=positions,mismatchNbs=mismatchNbs,IDs=[transcriptID]*len(positions)))
                                            except:                                                  
                                                if probes[probeIDs[probeIndex]].genes4splice is None:                                                                                            
                                                    probes[probeIDs[probeIndex]].genes4splice={}
                                                probes[probeIDs[probeIndex]].genes4splice[geneID]=ps_class.Targets(repetition=[1]*len(positions),indexes=[-1]*len(positions),positions=positions,mismatchNbs=mismatchNbs,IDs=[transcriptID]*len(positions))                                                                                                                                                                                                                                                                                                                                
                                            
                                            inExonIndexes[geneIndex][probeIndex]=2
                                            inExonIndex[probeIndex]=2
                                            nisIndex[probeIndex]=0                                                                                                                        
                                            
            
            '''SCAN ALL THE PROBES OF THE CURRENT PROBE SET TO FIND ALL TARGETED GENES'''
            
            #VARIABLES CONTAINING PROBE LOCALISATION FOR EACH TARGETED GENE
            #initiate thes varaibles with the number of already targeted genes with probes in their exons of in their splice junctions            
            inExonProbeNb=[0]*len(geneIDs)        
            inSpliceProbeNb=[0]*len(geneIDs)
            inGeneProbeNb=[0]*len(geneIDs)
            inIntronProbeNb=[0]*len(geneIDs)
            upProbeNb=[0]*len(geneIDs)
            downProbeNb=[0]*len(geneIDs)        
            
            # INITIATE THE VARIABLES CONTAINING GENE INFORMATION FOR GENES WITH PROBES IN THEIR EXONS OR IN THEIR SPLICE JUNCTIONS
            #initiate TargtedGenes object attributes
            if geneIDs!=[]:        
                for i in range(len(geneIDs)):                  
                    genes[i].probeIndexes=range(probeNb)
                    genes[i].firstStructureIDs=[[]]*probeNb                                        
                    genes[i].firstStructureIndexes=[[]]*probeNb   
                    genes[i].sndStructureIndexes=[[]]*probeNb
                    genes[i].firstStructureGroups=[[]]*probeNb   
                    genes[i].sndStructureGroups=[[]]*probeNb                        
                    genes[i].probeLocalisations=[[]]*probeNb
                    #nis tag is set to all genes                                
                    genes[i].probePositions=[[]]*probeNb            
                    genes[i].probeRepetitionNbs=[[]]*probeNb            
                    genes[i].probeStrands=[[]]*probeNb            
                    genes[i].probeMismatchNbs=[[]]*probeNb
                if nisIndex.count(1)>0:    
                    for probeIndex in range(probeNb):
                        if nisIndex[probeIndex]==1:
                            for i in range(len(geneIDs)):
                                genes[i].probeLocalisations[probeIndex]=['n']    
            
            #DEFINE AND INITIATE THE TargetedGene OBJECT TO BE USED FOR GENES WITH NO PROBES IN THEIR EXONS OR IN THEIR SPLICE JUNCTIONS
            probeLocalisations=[[]]*probeNb
            if nisIndex.count(1)>0:
                for probeIndex in range(probeNb):
                        if nisIndex[probeIndex]==1:
                            for geneIndex in range(len(geneIDs)):
                                probeLocalisations[probeIndex]=['n']            
            newGene=ps_class.TargetedGene(probeIndexes=range(probeNb),\
                                         firstStructureIDs=[[]]*probeNb,\
                                         firstStructureIndexes=[[]]*probeNb,\
                                         sndStructureIndexes=[[]]*probeNb,\
                                         firstStructureGroups=[[]]*probeNb,\
                                         sndStructureGroups=[[]]*probeNb,\
                                         probeLocalisations=probeLocalisations,\
                                         probePositions=[[]]*probeNb,\
                                         probeRepetitionNbs=[[]]*probeNb,\
                                         probeStrands=[[]]*probeNb,\
                                         probeMismatchNbs=[[]]*probeNb)
            
                                          
                                                                
            #SCAN ALL THE PROBES TO FILL PROBE LOCALISATION AND GENE INFORMATION VARIABLES                    
            for probeIndex in range(probeNb):
                probe=probes[probeIDs[probeIndex]]    
                #update probe if in splice junction or nis information has been changed                                            
                if nisIndex[probeIndex]==1 or inExonIndex[probeIndex]==2:                                            
                    probeDB[probeIDs[probeIndex]]=cPickle.dumps(probe,protocol=-1)                                          
                if probe.genes4exon !={}:
                    # targetGenes have been already been created                            
                    for geneID in probe.genes4exon.keys():                    
                        geneIndex=geneIDs.index(geneID)                    
                        inExonProbeNb[geneIndex]+=1
                        inGeneProbeNb[geneIndex]+=1
    #                    if genes[geneIndex].firstStructureIndexes[probeIndex]==[]:                        
                        genes[geneIndex].firstStructureIDs[probeIndex]=[probe.genes4exon[geneID].IDs[0]]                        
                        genes[geneIndex].firstStructureIndexes[probeIndex]=[probe.genes4exon[geneID].indexes[0]]                                                                                                    
                        genes[geneIndex].firstStructureGroups[probeIndex]=[probe.genes4exon[geneID].groups[0]]
                        genes[geneIndex].probeLocalisations[probeIndex]=['e']                                                                                                                           
                        genes[geneIndex].probePositions[probeIndex]=[probe.genes4exon[geneID].positions[0]]
                        genes[geneIndex].probeRepetitionNbs[probeIndex]=[probe.genes4exon[geneID].repetition[0]]
                        genes[geneIndex].probeStrands[probeIndex]=[probe.genes4exon[geneID].strands[0]]
                        genes[geneIndex].probeMismatchNbs[probeIndex]=[probe.genes4exon[geneID].mismatchNbs[0]]
                        if len(probe.genes4exon[geneID].repetition)>1:                                                                   
                            for i in range(1,len(probe.genes4exon[geneID].repetition)):                                                                 
                                genes[geneIndex].firstStructureIDs[probeIndex].append(probe.genes4exon[geneID].IDs[i])                        
                                genes[geneIndex].firstStructureIndexes[probeIndex].append(probe.genes4exon[geneID].indexes[i])                                                                                                    
                                genes[geneIndex].firstStructureGroups[probeIndex].append(probe.genes4exon[geneID].groups[i])
                                genes[geneIndex].probeLocalisations[probeIndex].append('e')                                                                                                                            
                                genes[geneIndex].probePositions[probeIndex].append(probe.genes4exon[geneID].positions[i])
                                genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4exon[geneID].repetition[i])
                                genes[geneIndex].probeStrands[probeIndex].append(probe.genes4exon[geneID].strands[i])
                                genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4exon[geneID].mismatchNbs[i])
                                                                                                                                                                    
                if probe.genes4splice !={} and probe.genes4splice is not None:
                    # targetGenes have been already been created                                               
                    for geneID in probe.genes4splice.keys():
                        geneIndex=geneIDs.index(geneID)
                        inSpliceProbeNb[geneIndex]+=1
                        inGeneProbeNb[geneIndex]+=1                    
                        if genes[geneIndex].firstStructureIndexes[probeIndex]==[]:                        
                            genes[geneIndex].firstStructureIDs[probeIndex]=[probe.genes4splice[geneID].IDs[0]]                                        
                            genes[geneIndex].firstStructureIndexes[probeIndex]=[probe.genes4splice[geneID].indexes[0]]                     
                            genes[geneIndex].firstStructureGroups[probeIndex]=[-1]                    
                            genes[geneIndex].probeLocalisations[probeIndex]=['s']
                            genes[geneIndex].probePositions[probeIndex]=[-1]                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex]=[probe.genes4splice[geneID].repetition[0]]                                            
                            genes[geneIndex].probeMismatchNbs[probeIndex]=[probe.genes4splice[geneID].mismatchNbs[0]]
                        else:
                            genes[geneIndex].firstStructureIDs[probeIndex].append(probe.genes4splice[geneID].IDs[0])                                        
                            genes[geneIndex].firstStructureIndexes[probeIndex].append(probe.genes4splice[geneID].indexes[0])                     
                            genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                    
                            genes[geneIndex].probeLocalisations[probeIndex].append('s')
                            genes[geneIndex].probePositions[probeIndex].append(-1)                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4splice[geneID].repetition[0])                                            
                            genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4splice[geneID].mismatchNbs[0])
                        if len(probe.genes4splice[geneID].repetition)>1:                                                                   
                            for i in range(1,len(probe.genes4splice[geneID].repetition)):
                                genes[geneIndex].firstStructureIDs[probeIndex].append(probe.genes4splice[geneID].IDs[i])                                        
                                genes[geneIndex].firstStructureIndexes[probeIndex].append(probe.genes4splice[geneID].indexes[i])                     
                                genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                    
                                genes[geneIndex].probeLocalisations[probeIndex].append('s')
                                genes[geneIndex].probePositions[probeIndex].append(-1)                                        
                                genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4splice[geneID].repetition[i])                                            
                                genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4splice[geneID].mismatchNbs[i])
                                                                                                                                                                                                                                                                                                
                if probe.genes4intron !={}:                                   
                    for geneID in probe.genes4intron.keys():                    
                        try:
                            #test if the gene has been already processed
                            geneIndex=geneIDs.index(geneID)                                            
                        except:
                            #create a new TargetGene       
                            geneIndex=len(geneIDs)                             
                            geneIDs.append(geneID)
                            genes.append(copypy.deepcopy(newGene))
                            inExonProbeNb.append(0)                            
                            inSpliceProbeNb.append(0)
                            inGeneProbeNb.append(0)
                            inIntronProbeNb.append(0)
                            upProbeNb.append(0)
                            downProbeNb.append(0)
                        inIntronProbeNb[geneIndex]+=1                    
                        if genes[geneIndex].firstStructureIndexes[probeIndex]==[]:
                            genes[geneIndex].firstStructureIDs[probeIndex]=[-1]
                            genes[geneIndex].firstStructureIndexes[probeIndex]=[probe.genes4intron[geneID].indexes[0]]                                        
                            genes[geneIndex].firstStructureGroups[probeIndex]=[-1]
                            genes[geneIndex].probePositions[probeIndex]=[probe.genes4intron[geneID].positions[0]]                                                        
                            genes[geneIndex].probeLocalisations[probeIndex]=['i']
                            genes[geneIndex].probeRepetitionNbs[probeIndex]=[probe.genes4intron[geneID].repetition[0]]
                            genes[geneIndex].probeStrands[probeIndex]=[probe.genes4intron[geneID].strands[0]]
                            genes[geneIndex].probeMismatchNbs[probeIndex]=[probe.genes4intron[geneID].mismatchNbs[0]]
                        else:                        
                            genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                            genes[geneIndex].firstStructureIndexes[probeIndex].append(probe.genes4intron[geneID].indexes[0])                                        
                            genes[geneIndex].firstStructureGroups[probeIndex].append(-1)
                            genes[geneIndex].probePositions[probeIndex].append(probe.genes4intron[geneID].positions[0])                                                        
                            genes[geneIndex].probeLocalisations[probeIndex].append('i')
                            genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4intron[geneID].repetition[0])
                            genes[geneIndex].probeStrands[probeIndex].append(probe.genes4intron[geneID].strands[0])
                            genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4intron[geneID].mismatchNbs[0])
                        if len(probe.genes4intron[geneID].repetition)>1:                                                                   
                            for i in range(1,len(probe.genes4intron[geneID].repetition)):
                                genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                                genes[geneIndex].firstStructureIndexes[probeIndex].append(probe.genes4intron[geneID].indexes[i])                                        
                                genes[geneIndex].firstStructureGroups[probeIndex].append(-1)
                                genes[geneIndex].probePositions[probeIndex].append(probe.genes4intron[geneID].positions[i])                                                        
                                genes[geneIndex].probeLocalisations[probeIndex].append('i')
                                genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4intron[geneID].repetition[i])
                                genes[geneIndex].probeStrands[probeIndex].append(probe.genes4intron[geneID].strands[i])
                                genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4intron[geneID].mismatchNbs[i])                                             
                                                                                                                             
                if probe.genes4up !={}:                                                       
                    for geneID in probe.genes4up.keys():
                        try:
                            #test if the gene has been already processed
                            geneIndex=geneIDs.index(geneID)                                            
                        except:       
                            #create a new TargetGene
                            geneIndex=len(geneIDs)             
                            geneIDs.append(geneID)                                                                                                         
                            genes.append(copypy.deepcopy(newGene))                                                                                                        
                            inExonProbeNb.append(0)                            
                            inSpliceProbeNb.append(0)
                            inGeneProbeNb.append(0)
                            inIntronProbeNb.append(0)
                            upProbeNb.append(0)
                            downProbeNb.append(0)                    
                        upProbeNb[geneIndex]+=1                                   
                        if genes[geneIndex].firstStructureIndexes[probeIndex]==[]:                        
                            genes[geneIndex].firstStructureIDs[probeIndex]=[-1]
                            genes[geneIndex].firstStructureIndexes[probeIndex]=[-1]
                            genes[geneIndex].firstStructureGroups[probeIndex]=[-1]                                            
                            genes[geneIndex].probePositions[probeIndex]=[probe.genes4up[geneID].positions[0]]
                            genes[geneIndex].probeLocalisations[probeIndex]=['u']                                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex]=[probe.genes4up[geneID].repetition[0]]                    
                            genes[geneIndex].probeStrands[probeIndex]=[probe.genes4up[geneID].strands[0]]
                            genes[geneIndex].probeMismatchNbs[probeIndex]=[probe.genes4up[geneID].mismatchNbs[0]]
                        else:
                            genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                            genes[geneIndex].firstStructureIndexes[probeIndex].append(-1)
                            genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                                            
                            genes[geneIndex].probePositions[probeIndex].append(probe.genes4up[geneID].positions[0])
                            genes[geneIndex].probeLocalisations[probeIndex].append('u')                                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4up[geneID].repetition[0])                    
                            genes[geneIndex].probeStrands[probeIndex].append(probe.genes4up[geneID].strands[0])
                            genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4up[geneID].mismatchNbs[0])
                        if len(probe.genes4up[geneID].repetition)>1:                                                                   
                            for i in range(1,len(probe.genes4up[geneID].repetition)):
                                genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                                genes[geneIndex].firstStructureIndexes[probeIndex].append(-1)
                                genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                                            
                                genes[geneIndex].probePositions[probeIndex].append(probe.genes4up[geneID].positions[i])
                                genes[geneIndex].probeLocalisations[probeIndex].append('u')                                                        
                                genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4up[geneID].repetition[i])                    
                                genes[geneIndex].probeStrands[probeIndex].append(probe.genes4up[geneID].strands[i])
                                genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4up[geneID].mismatchNbs[i])
                                                               
                        
                if probe.genes4down !={}:#                                    
                    for geneID in probe.genes4down.keys():
                        try:
                            #test if the gene has been already processed                 
                            geneIndex=geneIDs.index(geneID)                                                
                        except:
                            #create a new TargetGene
                            geneIndex=len(geneIDs)                    
                            geneIDs.append(geneID)                                                    
                            genes.append(copypy.deepcopy(newGene))
                            inExonProbeNb.append(0)
                            inSpliceProbeNb.append(0)
                            inGeneProbeNb.append(0)
                            inIntronProbeNb.append(0)
                            upProbeNb.append(0)
                            downProbeNb.append(0)                                             
                        downProbeNb[geneIndex]+=1                                 
                        if genes[geneIndex].firstStructureIndexes[probeIndex]==[]:       
                            genes[geneIndex].firstStructureIDs[probeIndex]=[-1]
                            genes[geneIndex].firstStructureIndexes[probeIndex]=[-1]
                            genes[geneIndex].firstStructureGroups[probeIndex]=[-1]                                            
                            genes[geneIndex].probePositions[probeIndex]=[probe.genes4down[geneID].positions[0]]
                            genes[geneIndex].probeLocalisations[probeIndex]=['d']                                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex]=[probe.genes4down[geneID].repetition[0]]                    
                            genes[geneIndex].probeStrands[probeIndex]=[probe.genes4down[geneID].strands[0]]
                            genes[geneIndex].probeMismatchNbs[probeIndex]=[probe.genes4down[geneID].mismatchNbs[0]]
                        else:
                            genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                            genes[geneIndex].firstStructureIndexes[probeIndex].append(-1)
                            genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                                            
                            genes[geneIndex].probePositions[probeIndex].append(probe.genes4down[geneID].positions[0])
                            genes[geneIndex].probeLocalisations[probeIndex].append('d')                                                        
                            genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4down[geneID].repetition[0])                    
                            genes[geneIndex].probeStrands[probeIndex].append(probe.genes4down[geneID].strands[0])
                            genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4down[geneID].mismatchNbs[0])
                        if len(probe.genes4down[geneID].repetition)>1:                                                                   
                            for i in range(1,len(probe.genes4down[geneID].repetition)):
                                genes[geneIndex].firstStructureIDs[probeIndex].append(-1)
                                genes[geneIndex].firstStructureIndexes[probeIndex].append(-1)
                                genes[geneIndex].firstStructureGroups[probeIndex].append(-1)                                            
                                genes[geneIndex].probePositions[probeIndex].append(probe.genes4down[geneID].positions[i])
                                genes[geneIndex].probeLocalisations[probeIndex].append('d')                                                        
                                genes[geneIndex].probeRepetitionNbs[probeIndex].append(probe.genes4down[geneID].repetition[i])                    
                                genes[geneIndex].probeStrands[probeIndex].append(probe.genes4down[geneID].strands[i])
                                genes[geneIndex].probeMismatchNbs[probeIndex].append(probe.genes4down[geneID].mismatchNbs[i])
                                                                    
                                                                                    
            '''CALCULATE PROBE SET ATTRIBUTES'''
            #Count the number of genes that have i probes either in exons or splices junctions,
            #or in introns or in up or down regions
            inExonGeneNb=[]        
            inIntronGeneNb=[]
            upGeneNb=[]
            downGeneNb=[]                                            
            for i in range(probeNb,0,-1):
                inExonGeneNb.append(inGeneProbeNb.count(i))
                inIntronGeneNb.append(inIntronProbeNb.count(i))
                upGeneNb.append(upProbeNb.count(i))
                downGeneNb.append(downProbeNb.count(i))
            #Count for each probe the number of positions located out of gene
            # or not in sequence                                             
            outNb=[]
            nisNb=[]
            for probeIndex in range(probeNb):           
                outNb.append(probes[probeIDs[probeIndex]].outOfGeneNb)            
                if probes[probeIDs[probeIndex]].notInSequence is None:
                    nisNb.append(0)
                else:
                    nisNb.append(probes[probeIDs[probeIndex]].notInSequence)                    
            
            #CALCULATE TargetedGenes ATTRIBUTES             
            for geneIndex in range(len(geneIDs)):
                #fill probe localisation (exon(e),intron(i),up(u),down(d),out(o))                    
                genes[geneIndex].fill_info()    
                exonIDSet=None
                exonIDSet=set()   
                #process transcript information                     
                foundTranscript=0        
                if dbType=='ensembl':
                    genes[geneIndex].ensemblGenes=[geneIDs[geneIndex]]
                elif dbType=='aceview':                                               
                        genes[geneIndex].ensemblGenes=cPickle.loads(ensByGeneDB[geneIDs[geneIndex]])                                                                        
                for probeIndex in range(len(genes[geneIndex].probeLocalisations)):
                    for j in range(len(genes[geneIndex].probeLocalisations[probeIndex])):
                        currStructure=genes[geneIndex].probeLocalisations[probeIndex][j]                        
                        if ['e','i'].count(currStructure):
                            foundTranscript=1                                                
                            if genes[geneIndex].firstStructureIDs[probeIndex][j]!=-1:
                                exonIDSet.add(genes[geneIndex].firstStructureIDs[probeIndex][j])                                                                                                            
                if foundTranscript==1:
                    if geneIDs[geneIndex][:3]!='GOP':                                                
                        transcriptsIn=set()
                        transcriptsOut=set()    
                        transcripts=cPickle.loads(transcriptsByGeneDB[geneIDs[geneIndex]])
                        if dbType=='ensembl':
                            for transcriptID in transcripts.IDs:
                                transcriptsOut.add(transcriptID)
                            for exonID in exonIDSet:                                     
                                transcripts=cPickle.loads(transcriptsByExonDB[exonID])                    
                                for transcriptID in transcripts.IDs:
                                    transcriptsIn.add(transcriptID)                        
                        elif dbType=='aceview':    
                            for transcriptID in transcripts:
                                transcriptsOut.add(transcriptID)                            
                            for exonID in exonIDSet:                                     
                                transcripts=cPickle.loads(transcriptsByExonDB[exonID])                    
                                for transcriptID in transcripts:
                                    transcriptsIn.add(transcriptID)            
                        TList=list(transcriptsIn)
                        TList.sort()        
                        genes[geneIndex].targetedTranscripts=TList
                        TList=list(transcriptsOut-transcriptsIn)
                        TList.sort()
                        genes[geneIndex].notTargetedTranscripts=TList
                    
                            
                                                                            
            '''SAVE PROBE SETS'''     
            #complete with 0 up to MmaxProbeNb to have vectors of same size
            addedZ=[]
            if probeNb<maxProbeNb:
                addedZ=[0]*(maxProbeNb-probeNb)                   
            if dbType=='ensembl':               
                probeset.ensExonGeneNbs=addedZ+inExonGeneNb
                probeset.ensIntronGeneNbs=addedZ+inIntronGeneNb
                probeset.ensUpGeneNbs=addedZ+upGeneNb
                probeset.ensDownGeneNbs=addedZ+downGeneNb
                probeset.ensOutProbeNbs=addedZ+outNb
                probeset.ensNisProbeNbs=addedZ+nisNb
                probeset.ensOurEnsGeneIDs=[]
            elif dbType=='aceview':
                probeset.aceExonGeneNbs=addedZ+inExonGeneNb
                probeset.aceIntronGeneNbs=addedZ+inIntronGeneNb
                probeset.aceUpGeneNbs=addedZ+upGeneNb
                probeset.aceDownGeneNbs=addedZ+downGeneNb
                probeset.aceOutProbeNbs=addedZ+outNb
                probeset.aceNisProbeNbs=addedZ+nisNb
                probeset.aceOurEnsGeneIDs=[]
            
            #ProbeNb and GeneNb variable are indexed on different probe number order
            probeNbList=range(probeNb,0,-1)            
            #search the genes with the highest nb of probes in exon or splice            
            #find the gene IDs of genes with probe(s) in exon or in splice junction
            #for i in range(probeNb)
            #start with the smallest number of probes (probeNbList(probeNb-1)=1)
            #and ends with the largest number of probes (probeNbList(0)=probeNb)
            for i in range(probeNb):                            
                if inExonGeneNb[i]>0:                
#                    for geneIndex in range(len(inExonProbeNb)):
#                        if inExonProbeNb[geneIndex]+inSpliceProbeNb[geneIndex]==probeNbList[i]:                        
                    for geneIndex in range(len(inGeneProbeNb)):
                        if inGeneProbeNb[geneIndex]==probeNbList[i]:                            
                            if dbType=='ensembl':
                                probeset.ourEnsGeneIDs.append(geneIDs[geneIndex])
                            elif dbType=='aceview':                        
                                probeset.ourAceGeneIDs.append(geneIDs[geneIndex])
                    if dbType=='ensembl':            
                        probeset.ourEnsGeneNb=inExonGeneNb[i]
                        probeset.ourEnsProbeNb=probeNbList[i]
                    elif dbType=='aceview':        
                        probeset.ourAceGeneNb=inExonGeneNb[i]
                        probeset.ourAceProbeNb=probeNbList[i]                    
                    break                
            if dbType=='ensembl':
                probeset.ourEnsGeneIDs.sort()                
            elif dbType=='aceview':                        
                probeset.ourAceGeneIDs.sort()
            #find the corespondant Ensembl gene for AceView gene                    
            if dbType=='aceview':
                if probeset.ourAceGeneIDs !=[]:
                    ourAceToEnsGeneIDs=[]
                    ourAceToEnsGeneNbs=[]                    
                    for aceGene in probeset.ourAceGeneIDs:
                        try:                                        
                            ourAceToEnsGeneIDs.extend(cPickle.loads(ensByGeneDB[aceGene]))                                                    
                            ourAceToEnsGeneNbs.append(len(cPickle.loads(ensByGeneDB[aceGene])))
                        except: 
                            pass                                                                       
                    ourAceToEnsGeneIDs.sort()                    
                    probeset.ourAceToEnsGeneIDs=ourAceToEnsGeneIDs
                    probeset.ourAceToEnsGeneNbs=ourAceToEnsGeneNbs  
            
            ''' Put in a dictionary indexed on the number of probes (from 1 to probeNb) the genes that have the
            corresponding number of probe in exons or in splice junctions
            Put all genes that have also probes in introns, up and down regions in the index 0'''                  
            if dbType=='ensembl':
                probeset.ensIDbyProbeNb={}
                for i in range(probeNb):            
                    if inExonGeneNb[i]>0:     
                        for geneIndex in range(len(inGeneProbeNb)):
                            if inGeneProbeNb[geneIndex]==probeNbList[i]:
                                try:                                
                                    probeset.ensGenesByProbeNb[probeNbList[i]][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.ensGenesByProbeNb is None:
                                        probeset.ensGenesByProbeNb={}                                                             
                                    probeset.ensGenesByProbeNb[probeNbList[i]]={}
                                    probeset.ensGenesByProbeNb[probeNbList[i]][geneIDs[geneIndex]]=genes[geneIndex]
                    if inIntronGeneNb[i]>0:
                        for geneIndex in range(len(inIntronProbeNb)):
                            if inIntronProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.ensGenesByProbeNb is None:
                                        probeset.ensGenesByProbeNb={}
                                    probeset.ensGenesByProbeNb[0]={}    
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]                                
                    if upGeneNb[i]>0:
                        for geneIndex in range(len(upProbeNb)):
                            if upProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.ensGenesByProbeNb is None:
                                        probeset.ensGenesByProbeNb={}
                                    probeset.ensGenesByProbeNb[0]={}
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]                            
                    if downGeneNb[i]>0:
                        for geneIndex in range(len(downProbeNb)):
                            if downProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.ensGenesByProbeNb is None:
                                        probeset.ensGenesByProbeNb={}
                                    probeset.ensGenesByProbeNb[0]={}
                                    probeset.ensGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]         
                if probeset.ensGeneIDs is not None:
                    ensemblGenesSet=set(probeset.ensGeneIDs)
                else:
                    ensemblGenesSet=set([])
                if probeset.ourEnsGeneIDs is not None:   
                    ourGenesSet=set(probeset.ourEnsGeneIDs)
                else:
                    ourGenesSet=set([])
                GList=list(ensemblGenesSet-ourGenesSet)
                GList.sort()    
                probeset.ensGeneIDs4Ens=GList
                GList=list(ourGenesSet-ensemblGenesSet)
                GList.sort()
                probeset.ourEnsGeneIDs4Ens=GList
                GList=list(ourGenesSet.intersection(ensemblGenesSet))
                GList.sort()
                probeset.commonGeneIDs4Ens=GList
            elif dbType=='aceview':
                probeset.aceIDbyProbeNb={}
                #for i in range(probeNb):            
                for i in range(probeNb,0,-1):
                    if inExonGeneNb[i]>0:
                        for geneIndex in range(len(inGeneProbeNb)):
                            if inGeneProbeNb[geneIndex]==probeNbList[i]:        
                                try:                                
                                    probeset.aceGenesByProbeNb[probeNbList[i]][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.aceGenesByProbeNb is None:
                                        probeset.aceGenesByProbeNb={}                                                             
                                    probeset.aceGenesByProbeNb[probeNbList[i]]={}
                                    probeset.aceGenesByProbeNb[probeNbList[i]][geneIDs[geneIndex]]=genes[geneIndex]
                    if inIntronGeneNb[i]>0:
                        for geneIndex in range(len(inIntronProbeNb)):
                            if inIntronProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.aceGenesByProbeNb is None:
                                        probeset.aceGenesByProbeNb={}
                                    probeset.aceGenesByProbeNb[0]={}    
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]                                
                    if upGeneNb[i]>0:
                        for geneIndex in range(len(upProbeNb)):
                            if upProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.aceGenesByProbeNb is None:
                                        probeset.aceGenesByProbeNb={}
                                    probeset.aceGenesByProbeNb[0]={}
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]                            
                    if downGeneNb[i]>0:
                        for geneIndex in range(len(downProbeNb)):
                            if downProbeNb[geneIndex]>0:
                                try:                                
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                                except:
                                    if probeset.aceGenesByProbeNb is None:
                                        probeset.aceGenesByProbeNb={}
                                    probeset.aceGenesByProbeNb[0]={}
                                    probeset.aceGenesByProbeNb[0][geneIDs[geneIndex]]=genes[geneIndex]
                if probeset.ensGeneIDs is not None:
                    ensemblGenesSet=set(probeset.ensGeneIDs)
                else:
                    ensemblGenesSet=set([])
                if probeset.ourAceToEnsGeneIDs is not None:   
                    ourGenesSet=set(probeset.ourAceToEnsGeneIDs)
                else:
                    ourGenesSet=set([])
                GList=list(ensemblGenesSet-ourGenesSet)
                GList.sort()    
                probeset.ensGeneIDs4Ace=GList
                GList=list(ourGenesSet-ensemblGenesSet)
                GList.sort()
                probeset.ourEnsGeneIDs4Ace=GList
                GList=list(ourGenesSet.intersection(ensemblGenesSet))
                GList.sort()
                probeset.commonGeneIDs4Ace=GList   
        probesetDB[probesetID]=cPickle.dumps(probeset,protocol=-1)            
    logFile.close()
    probeDB.close()
    probesetDB.close()
    transcriptsByGeneDB.close()
    transcriptSequencesDB.close()
    transcriptsByExonDB.close()
    if dbType=='aceview':
        ensByGeneDB.close()
    

        
                
def make_probeset_list(species,chipName): 
    """
    **Make a dump of probesets.**
    ::
            
        Atributes:
            species --- species name
            chipName -- chip name
            
        Input:    
            %chip_probeset.bkdb
        
        Output:
            %chip_probeset.dump
        
    """
    
    # open the chip probeset data base 
    path=os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.bkdb"
    print path
    probesetDB=bsddb.btopen(path)
    
    
    
    # empty variables
    probesetIDs=[]
    probesetIndex=[]
    probeNb=[]
    ensGeneIDnb=[]
    affyGeneIDnb=[]                
    ensExonGeneNbs=[]
    ensIntronGeneNbs=[]
    ensUpGeneNbs=[]
    ensDownGeneNbs=[]
    ensOutProbeNbs=[]
    ensNisProbeNbs=[]        
    ourEnsGeneNb=[]
    ourEnsProbeNb=[]
    ensGeneIDs4EnsNb=[]
    ourEnsGeneIDs4EnsNb=[]    
    commonGeneIDs4EnsNb=[]
    ourEnsTargetedExons=[]
    ourEnsTargetedGroups=[]
    ourEnsTargetedTranscripts=[]
    ourEnsNotTargetedTranscripts=[]
      
    aceExonGeneNbs=[]
    aceIntronGeneNbs=[]
    aceUpGeneNbs=[]
    aceDownGeneNbs=[]
    aceOutProbeNbs=[]
    aceNisProbeNbs=[]              
    ourAceGeneNb=[]
    ourAceProbeNb=[]
                
    ensGeneIDs4AceNb=[]
    ourEnsGeneIDs4AceNb=[]
    commonGeneIDs4AceNb=[]
    ourAceTargetedExons=[]
    ourAceTargetedGroups=[]
    ourAceTargetedTranscripts=[]
    ourAceNotTargetedTranscripts=[]
    
    rank=-1
#    print 'processing data ...'
    for probesetID in probesetDB.keys():
#    for probesetID in ['92573_at']:
        rank +=1
#        print rank,': ',probesetID
        probeset=cPickle.loads(probesetDB[probesetID])
        probesetIDs.append(probeset.probesetID)
        probesetIndex.append(probeset.probesetIndex)
        probeNb.append(probeset.probeNb)
        # ENSEMBL INFO  
        if probeset.ensGeneIDs is not None:      
            ensGeneIDnb.append(len(probeset.ensGeneIDs))
        else:
            ensGeneIDnb.append(0)
        if probeset.affyGeneIDs is not None:        
            affyGeneIDnb.append(len(probeset.affyGeneIDs))
        else:
            affyGeneIDnb.append(0)                                
        ensExonGeneNbs.append(probeset.ensExonGeneNbs)
        ensIntronGeneNbs.append(probeset.ensIntronGeneNbs)
        ensUpGeneNbs.append(probeset.ensUpGeneNbs)
        ensDownGeneNbs.append(probeset.ensDownGeneNbs)
        ensOutProbeNbs.append(probeset.ensOutProbeNbs)
        ensNisProbeNbs.append(probeset.ensNisProbeNbs)        
        ourEnsGeneNb.append(probeset.ourEnsGeneNb)
        ourEnsProbeNb.append(probeset.ourEnsProbeNb)
        if probeset.ensGeneIDs4Ens is not None:
            ensGeneIDs4EnsNb.append(len(probeset.ensGeneIDs4Ens))
        else:
            ensGeneIDs4EnsNb.append(0)
        if probeset.ourEnsGeneIDs4Ens is not None:    
            ourEnsGeneIDs4EnsNb.append(len(probeset.ourEnsGeneIDs4Ens))
        else:
            ourEnsGeneIDs4EnsNb.append(0)
        if probeset.commonGeneIDs4Ens is not None:    
            commonGeneIDs4EnsNb.append(len(probeset.commonGeneIDs4Ens))
        else:
            commonGeneIDs4EnsNb.append(0)
        currExons={};        
        if probeset.ourEnsProbeNb>0:            
            genes=probeset.ensGenesByProbeNb[probeset.ourEnsProbeNb]            
            for gene in genes.keys():                          
                exons={};
                for probeIndexe in range(len(genes[gene].probeIndexes)):
                    if genes[gene].firstStructureGroups[probeIndexe]!=[]:                            
                        foundExons=[];
                        for item in range(len(genes[gene].firstStructureGroups[probeIndexe])):
                            if genes[gene].probeLocalisations[probeIndexe][item]=='e':
                                currExon=genes[gene].firstStructureIDs[probeIndexe][item]
                                if foundExons.count(currExon)==0:
                                    foundExons.append(currExon)
                                    if currExon in exons.keys():
                                        exons[currExon]=exons[currExon]+1
                                    else:
                                        exons[currExon]=1                                     
                currExons[gene]=exons                                                                      
        ourEnsTargetedExons.append(currExons)               
        currGroups={};        
        if probeset.ourEnsProbeNb>0:            
            genes=probeset.ensGenesByProbeNb[probeset.ourEnsProbeNb]            
            for gene in genes.keys():                          
                groups={};
                for probeIndexe in range(len(genes[gene].probeIndexes)):
                    if genes[gene].firstStructureGroups[probeIndexe]!=[]:                            
                        foundGroup=[];
                        for item in range(len(genes[gene].firstStructureGroups[probeIndexe])):
                            if genes[gene].probeLocalisations[probeIndexe][item]=='e':
                                currGroup=genes[gene].firstStructureGroups[probeIndexe][item]
                                if foundGroup.count(currGroup)==0:
                                    foundGroup.append(currGroup)
                                    if currGroup in groups.keys():
                                        groups[currGroup]=groups[currGroup]+1
                                    else:
                                        groups[currGroup]=1                                     
                currGroups[gene]=groups                                                                      
        ourEnsTargetedGroups.append(currGroups)           
        transcripts=[];
        notTranscripts=[];
        if probeset.ourEnsProbeNb>0:            
            genes=probeset.ensGenesByProbeNb[probeset.ourEnsProbeNb]
            transcripts=[];
            for gene in genes.keys():               
                if genes[gene].targetedTranscripts is not None:
                    transcripts.extend(genes[gene].targetedTranscripts)                                                                                                                                            
                if genes[gene].notTargetedTranscripts is not None:    
                    notTranscripts.extend(genes[gene].notTargetedTranscripts)                
        ourEnsTargetedTranscripts.append(transcripts)       
        ourEnsNotTargetedTranscripts.append(notTranscripts)
            
        
        #ACEVIEW INFO     
        aceExonGeneNbs.append(probeset.aceExonGeneNbs)
        aceIntronGeneNbs.append(probeset.aceIntronGeneNbs)
        aceUpGeneNbs.append(probeset.aceUpGeneNbs)
        aceDownGeneNbs.append(probeset.aceDownGeneNbs)
        aceOutProbeNbs.append(probeset.aceOutProbeNbs)
        aceNisProbeNbs.append(probeset.aceNisProbeNbs)        
        ourAceGeneNb.append(probeset.ourAceGeneNb)
        ourAceProbeNb.append(probeset.ourAceProbeNb)
                     
        if probeset.ensGeneIDs4Ace is not None:                        
            ensGeneIDs4AceNb.append(len(probeset.ensGeneIDs4Ace))
        else:
            ensGeneIDs4AceNb.append(0)
        if probeset.ourEnsGeneIDs4Ace is not None:    
            ourEnsGeneIDs4AceNb.append(len(probeset.ourEnsGeneIDs4Ace))
        else:
            ourEnsGeneIDs4AceNb.append(0)
        if probeset.commonGeneIDs4Ace is not None:    
            commonGeneIDs4AceNb.append(len(probeset.commonGeneIDs4Ace))
        else:
            commonGeneIDs4AceNb.append(0)            
            
        currExons={};        
        if probeset.ourAceProbeNb>0:                
            genes=probeset.aceGenesByProbeNb[probeset.ourAceProbeNb]                        
            for gene in genes.keys():                          
                exons={};
                for probeIndexe in range(len(genes[gene].probeIndexes)):
                    if genes[gene].firstStructureGroups[probeIndexe]!=[]:                            
                        foundExons=[];
                        for item in range(len(genes[gene].firstStructureGroups[probeIndexe])):
                            if genes[gene].probeLocalisations[probeIndexe][item]=='e':
                                currExon=genes[gene].firstStructureIDs[probeIndexe][item]
                                if foundExons.count(currExon)==0:
                                    foundExons.append(currExon)
                                    if currExon in exons.keys():
                                        exons[currExon]=exons[currExon]+1
                                    else:
                                        exons[currExon]=1                                     
                currExons[gene]=exons                                                                      
        ourAceTargetedExons.append(currExons)               
        
        currGroups={};        
        if probeset.ourAceProbeNb>0:            
            genes=probeset.aceGenesByProbeNb[probeset.ourAceProbeNb]            
            for gene in genes.keys():                          
                groups={};
                for probeIndexe in range(len(genes[gene].probeIndexes)):
                    if genes[gene].firstStructureGroups[probeIndexe]!=[]:                            
                        foundGroup=[];
                        for item in range(len(genes[gene].firstStructureGroups[probeIndexe])):
                            if genes[gene].probeLocalisations[probeIndexe][item]=='e':
                                currGroup=genes[gene].firstStructureGroups[probeIndexe][item]
                                if foundGroup.count(currGroup)==0:
                                    foundGroup.append(currGroup)
                                    if currGroup in groups.keys():
                                        groups[currGroup]=groups[currGroup]+1
                                    else:
                                        groups[currGroup]=1                                     
                currGroups[gene]=groups                                                                      
        ourAceTargetedGroups.append(currGroups)       
                                
        transcripts=[];
        notTranscripts=[];
        if probeset.ourAceProbeNb>0:            
            genes=probeset.aceGenesByProbeNb[probeset.ourAceProbeNb]
            transcripts=[];
            for gene in genes.keys():                          
                if genes[gene].targetedTranscripts is not None:
                    if genes[gene].targetedTranscripts!=[]:                    
                        for letter in genes[gene].targetedTranscripts:
                            transcripts.append(gene+'.'+letter)                                                                                                                                                                                        
                if genes[gene].notTargetedTranscripts is not None:
                    if genes[gene].notTargetedTranscripts!=[]:
                        for letter in genes[gene].notTargetedTranscripts:
                            notTranscripts.append(gene+'.'+letter)                                    
        ourAceTargetedTranscripts.append(transcripts)       
        ourAceNotTargetedTranscripts.append(notTranscripts)
                                    
        
#        break
        
    probesetList=ps_class.ProbesetList(probesetID=probesetIDs,probesetIndex=probesetIndex,probeNb=probeNb,ensGeneIDnb=ensGeneIDnb,affyGeneIDnb=affyGeneIDnb,\
                 ensExonGeneNbs=ensExonGeneNbs,ensIntronGeneNbs=ensIntronGeneNbs,ensUpGeneNbs=ensUpGeneNbs,ensDownGeneNbs=ensDownGeneNbs,\
                 ensOutProbeNbs=ensOutProbeNbs,ensNisProbeNbs=ensNisProbeNbs,\
                 aceExonGeneNbs=aceExonGeneNbs,aceIntronGeneNbs=aceIntronGeneNbs,aceUpGeneNbs=aceUpGeneNbs,aceDownGeneNbs=aceDownGeneNbs,\
                 aceOutProbeNbs=aceOutProbeNbs,aceNisProbeNbs=aceNisProbeNbs,\
                 ensGeneIDs4EnsNb=ensGeneIDs4EnsNb, ourEnsGeneNb=ourEnsGeneNb,ourEnsProbeNb=ourEnsProbeNb,\
                 ensGeneIDs4AceNb=ensGeneIDs4AceNb, ourAceGeneNb=ourAceGeneNb,ourAceProbeNb=ourAceProbeNb,\
                 ourEnsGeneIDs4EnsNb=ourEnsGeneIDs4EnsNb,commonGeneIDs4EnsNb= commonGeneIDs4EnsNb,\
                 ourEnsGeneIDs4AceNB=ourEnsGeneIDs4AceNb,commonGeneIDs4AceNb=commonGeneIDs4AceNb,\
                 ourEnsTargetedExons=ourEnsTargetedExons,\
                 ourEnsTargetedGroups=ourEnsTargetedGroups,\
                 ourEnsTargetedTranscripts=ourEnsTargetedTranscripts, ourEnsNotTargetedTranscripts=ourEnsNotTargetedTranscripts,\
                 ourAceTargetedExons=ourAceTargetedExons,\
                 ourAceTargetedGroups=ourAceTargetedGroups,\
                 ourAceTargetedTranscripts=ourAceTargetedTranscripts,ourAceNotTargetedTranscripts=ourAceNotTargetedTranscripts)    
    probesetDB.close()        
    os.chdir(os.environ['PYDATA']+"/" + species+"/")    
    outFile=open(os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.dump",'w')
#    print 'data processd'
    print os.getcwd()
    cPickle.dump(probesetList,outFile)
    outFile.close()

def probesets_by_gene(species,chipName,log=0):
    """
    **Find probe sets targeting genes.**
    ::
    
        Arguments:
            species --- species name
            chipName -- chip name
            
        Inputs:
            %chip_probeset_bkdb
            
        Output:    
            create %chip_probesets_by_gene.bkdb (for Ensembl)
            create %chip_probesets_by_gene.bkdb (for AceView)
        
    """
    
    t1=time.time()
    # open the chip probeset data base 
    path=os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.bkdb"
    
    print path
    probesetDB=bsddb.btopen(path,'r')
    
    #create the probesets_by_gene DBs
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_probesets_by_gene.bkdb"
    print path
    if os.path.exists(path):
        os.remove(path)
    ensemblDB=bsddb.btopen(path,'w')
    
    
    path=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probe.bkdb"
    if os.path.exists(path):        
        path=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probesets_by_gene.bkdb"
        if os.path.exists(path):
            os.remove(path)            
        aceviewDB=bsddb.btopen(path,'w')
        aceFlag=1
    else:
        aceFlag=0
    
    
    # empty dictionnaries whose keys are gene ID (ensembl or aceview)
    ensProbesetIndexes={}
    aceProbesetIndexes={}
    ensProbesetIDs={}
    aceProbesetIDs={}
    ensTargetedExons={}
    ensTargetedGroups={}
    ensTargetedTranscripts={}
    aceTargetedExons={}
    aceTargetedGroups={}
    aceTargetedTranscripts={}
    ensNotTargetedTranscripts={}
    aceNotTargetedTranscripts={}
    ensGeneIDs=[];
    aceGeneIDs=[];
    ensInGeneProbeNbs={};    
    aceInGeneProbeNbs={};
    ensNotInExonProbeNbs={};    
    aceNotInExonProbeNbs={};
    probesetIDs=probesetDB.keys()
    for probesetID in probesetIDs:
        # process the current probeset
        probeset=cPickle.loads(probesetDB[probesetID])
        #process Ensembl genes         
        for probeNb in probeset.ensGenesByProbeNb.keys():            
            #process the current gene
            for targetedGeneID in probeset.ensGenesByProbeNb[probeNb]:                                 
                targetedGene=probeset.ensGenesByProbeNb[probeNb][targetedGeneID]                
                probeIndexes=targetedGene.probeIndexes
                firstStructureIDs=targetedGene.firstStructureIDs
                firstStructureGroups=targetedGene.firstStructureGroups
                probeLocalisations=targetedGene.probeLocalisations
                inGeneProbeNb=array(targetedGene.inExonProbeNb)+array(targetedGene.inSpliceProbeNb)+array(targetedGene.inIntronProbeNb)+array(targetedGene.upProbeNb)+array(targetedGene.downProbeNb)                                
                inTranscriptProbeNb=array(targetedGene.inExonProbeNb)+array(targetedGene.inSpliceProbeNb)                
                notInExonProbeNb=len(set(inGeneProbeNb.nonzero()[0]).difference(set(inTranscriptProbeNb.nonzero()[0])))
                inGeneProbeNb=len(inGeneProbeNb.nonzero()[0])                
                
                #new gene                            
                if ensGeneIDs.count(targetedGeneID)==0:                    
                    ensGeneIDs.append(targetedGeneID)
                    ensProbesetIDs[targetedGeneID]={}
                    ensProbesetIDs[targetedGeneID][probeNb]=[probeset.probesetID]
                    ensProbesetIndexes[targetedGeneID]={}
                    ensProbesetIndexes[targetedGeneID][probeNb]=[probeset.probesetIndex]
                    ensTargetedExons[targetedGeneID]={}
                    ensTargetedGroups[targetedGeneID]={}
                    ensTargetedTranscripts[targetedGeneID]={}
                    ensNotTargetedTranscripts[targetedGeneID]={}
                    if probeNb==0:
#                        if targetedGeneID=='ENSMUSG00000000131':                            
#                            print targetedGene.targetedTranscripts
#                            print targetedGene.notTargetedTranscripts
                        ensTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                        ensTargetedExons[targetedGeneID][probeNb]=[[]]
                        ensTargetedGroups[targetedGeneID][probeNb]=[[]]
                        if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                            allTranscripts=[]
                        else:  
                            allTranscripts=targetedGene.targetedTranscripts[:]
                        if targetedGene.notTargetedTranscripts is not None:    
                            allTranscripts.extend(targetedGene.notTargetedTranscripts)
                        if allTranscripts==[]:
                            ensNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                        else:                                                                                                                    
                            ensNotTargetedTranscripts[targetedGeneID][probeNb]=allTranscripts
#                        if targetedGeneID=='ENSMUSG00000000131':
#                            print allTranscripts
#                            print ensNotTargetedTranscripts[targetedGeneID][probeNb]                            
                    else:
                        if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                            ensTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                        else:
                            ensTargetedTranscripts[targetedGeneID][probeNb]=[targetedGene.targetedTranscripts]
                        if targetedGene.notTargetedTranscripts is None or targetedGene.notTargetedTranscripts==[]:
                            ensNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                        else:
                            ensNotTargetedTranscripts[targetedGeneID][probeNb]=[targetedGene.notTargetedTranscripts]
                    ensInGeneProbeNbs[targetedGeneID]={}                    
                    ensInGeneProbeNbs[targetedGeneID][probeNb]= [inGeneProbeNb]
                    ensNotInExonProbeNbs[targetedGeneID]={}                    
                    ensNotInExonProbeNbs[targetedGeneID][probeNb]= [notInExonProbeNb]
                        
                    exons={};
                    for probeIndexe in range(len(probeIndexes)):
                        if firstStructureGroups[probeIndexe]!=[]:                            
                            foundExons=[];
                            for item in range(len(firstStructureGroups[probeIndexe])):
                                if probeLocalisations[probeIndexe][item]=='e':
                                    currExon=firstStructureIDs[probeIndexe][item]
                                    if foundExons.count(currExon)==0:
                                        foundExons.append(currExon)
                                        if currExon in exons.keys():
                                            exons[currExon]=exons[currExon]+1
                                        else:
                                            exons[currExon]=1                                     
                    ensTargetedExons[targetedGeneID][probeNb]=[exons]
                                        
                    groups={};
                    for probeIndexe in range(len(probeIndexes)):
                        if firstStructureGroups[probeIndexe]!=[]:                            
                            foundGroup=[];
                            for item in range(len(firstStructureGroups[probeIndexe])):
                                if probeLocalisations[probeIndexe][item]=='e':
                                    currGroup=firstStructureGroups[probeIndexe][item]
                                    if foundGroup.count(currGroup)==0:
                                        foundGroup.append(currGroup)
                                        if currGroup in groups.keys():
                                            groups[currGroup]=groups[currGroup]+1
                                        else:
                                            groups[currGroup]=1                                     
                    ensTargetedGroups[targetedGeneID][probeNb]=[groups]
                    
                    
                #existing gene                                            
                else:
                    #the gene exist but the key probeNb is may be not yet used                    
                    try:                                        
                        ensProbesetIDs[targetedGeneID][probeNb].append(probeset.probesetID)                        
                        ensProbesetIndexes[targetedGeneID][probeNb].append(probeset.probesetIndex)                        
                        if probeNb==0:                             
#                            if targetedGeneID=='ENSMUSG00000000131':                            
#                                print targetedGene.targetedTranscripts
#                                print targetedGene.notTargetedTranscripts
                            ensTargetedExons[targetedGeneID][probeNb].append([])                            
                            ensTargetedGroups[targetedGeneID][probeNb].append([])                            
                            ensTargetedTranscripts[targetedGeneID][probeNb].append([])                            
                            if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                                allTranscripts=[]
                            else:                                                                
                                allTranscripts=targetedGene.targetedTranscripts[:]
                            if targetedGene.notTargetedTranscripts is not None:     
                                allTranscripts.extend(targetedGene.notTargetedTranscripts)
                            ensNotTargetedTranscripts[targetedGeneID][probeNb].append(allTranscripts)
#                            if targetedGeneID=='ENSMUSG00000000131':
#                                print allTranscripts
#                                print ensNotTargetedTranscripts[targetedGeneID][probeNb]
                                                                                
                        else:                                                
                            if targetedGene.targetedTranscripts is not None:                          
                                ensTargetedTranscripts[targetedGeneID][probeNb].append(targetedGene.targetedTranscripts)                                
                            if targetedGene.notTargetedTranscripts is not None:    
                                ensNotTargetedTranscripts[targetedGeneID][probeNb].append(targetedGene.notTargetedTranscripts)                                                            
                            exons={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundExons=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currExon=firstStructureIDs[probeIndexe][item]
                                            if foundExons.count(currExon)==0:
                                                foundExons.append(currExon)
                                                if currExon in exons.keys():
                                                    exons[currExon]=exons[currExon]+1
                                                else:
                                                    exons[currExon]=1                                     
                            ensTargetedExons[targetedGeneID][probeNb].append(exons)    
                            
                            groups={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundGroup=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currGroup=firstStructureGroups[probeIndexe][item]
                                            if foundGroup.count(currGroup)==0:
                                                foundGroup.append(currGroup)
                                                if currGroup in groups.keys():
                                                    groups[currGroup]=groups[currGroup]+1
                                                else:
                                                    groups[currGroup]=1                                     
                            ensTargetedGroups[targetedGeneID][probeNb].append(groups)                                
                        ensInGeneProbeNbs[targetedGeneID][probeNb].append(inGeneProbeNb)                        
                        ensNotInExonProbeNbs[targetedGeneID][probeNb].append(notInExonProbeNb)
                    except:
                        #existing gene but new key
                        #print 'except for '+targetedGeneID+' '+probeset.probesetID 
                        try:
                            #
                            if ensProbesetIDs[targetedGeneID][probeNb]!=[]:                                
                                print 'error ensembl: '+repr(ensProbesetIDs[targetedGeneID][probeNb])
                                sys.exit()                            
                        except:
                            pass                                                                                                                                        
                        ensProbesetIDs[targetedGeneID][probeNb]=[probeset.probesetID]
                        ensProbesetIndexes[targetedGeneID][probeNb]=[probeset.probesetIndex]
                        if probeNb==0:
                            ensTargetedExons[targetedGeneID][probeNb]=[[]]
                            ensTargetedGroups[targetedGeneID][probeNb]=[[]]
                            ensTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                                allTranscripts=[]
                            else:        
                                allTranscripts=targetedGene.targetedTranscripts[:]
                            if targetedGene.notTargetedTranscripts is not None:    
                                allTranscripts.extend(targetedGene.notTargetedTranscripts)
                            if allTranscripts==[]:
                                ensNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            else:        
                                ensNotTargetedTranscripts[targetedGeneID][probeNb]=allTranscripts
                        else:                            
                            if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:    
                                ensTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            else:
                                ensTargetedTranscripts[targetedGeneID][probeNb]=[targetedGene.targetedTranscripts]
                            if targetedGene.notTargetedTranscripts is None or targetedGene.notTargetedTranscripts==[]:
                                ensNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            else:        
                                ensNotTargetedTranscripts[targetedGeneID][probeNb]=[targetedGene.notTargetedTranscripts]                            
                                
                            exons={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundExons=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currExon=firstStructureIDs[probeIndexe][item]
                                            if foundExons.count(currExon)==0:
                                                foundExons.append(currExon)
                                                if currExon in exons.keys():
                                                    exons[currExon]=exons[currExon]+1
                                                else:
                                                    exons[currExon]=1                                     
                            ensTargetedExons[targetedGeneID][probeNb]=[exons]    
                            
                            groups={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundGroup=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currGroup=firstStructureGroups[probeIndexe][item]
                                            if foundGroup.count(currGroup)==0:
                                                foundGroup.append(currGroup)
                                                if currGroup in groups.keys():
                                                    groups[currGroup]=groups[currGroup]+1
                                                else:
                                                    groups[currGroup]=1                                     
                            ensTargetedGroups[targetedGeneID][probeNb]=[groups]                                    
                        ensInGeneProbeNbs[targetedGeneID][probeNb]=[inGeneProbeNb]                    
                        ensNotInExonProbeNbs[targetedGeneID][probeNb]=[notInExonProbeNb]
                            
        #process AceView genes
        if aceFlag:                    
            for probeNb in probeset.aceGenesByProbeNb.keys():           
                for targetedGeneID in probeset.aceGenesByProbeNb[probeNb]:                
                    # process the current gene
                    targetedGene=probeset.aceGenesByProbeNb[probeNb][targetedGeneID]                
                    probeIndexes=targetedGene.probeIndexes
                    firstStructureIDs=targetedGene.firstStructureIDs
                    firstStructureGroups=targetedGene.firstStructureGroups
                    probeLocalisations=targetedGene.probeLocalisations
                    inGeneProbeNb=array(targetedGene.inExonProbeNb)+array(targetedGene.inSpliceProbeNb)+array(targetedGene.inIntronProbeNb)+array(targetedGene.upProbeNb)+array(targetedGene.downProbeNb)                                
                    inTranscriptProbeNb=array(targetedGene.inExonProbeNb)+array(targetedGene.inSpliceProbeNb)                
                    notInExonProbeNb=len(set(inGeneProbeNb.nonzero()[0]).difference(set(inTranscriptProbeNb.nonzero()[0])))
                    inGeneProbeNb=len(inGeneProbeNb.nonzero()[0])                
                    #new gene
                    if aceGeneIDs.count(targetedGeneID)==0:
                        aceGeneIDs.append(targetedGeneID)
                        aceProbesetIDs[targetedGeneID]={}
                        aceProbesetIDs[targetedGeneID][probeNb]=[probeset.probesetID]
                        aceProbesetIndexes[targetedGeneID]={}
                        aceProbesetIndexes[targetedGeneID][probeNb]=[probeset.probesetIndex]
                        aceTargetedExons[targetedGeneID]={}
                        aceTargetedGroups[targetedGeneID]={}
                        aceTargetedTranscripts[targetedGeneID]={}
                        aceNotTargetedTranscripts[targetedGeneID]={}
                        if probeNb==0:
                            aceTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            aceTargetedExons[targetedGeneID][probeNb]=[[]]
                            aceTargetedGroups[targetedGeneID][probeNb]=[[]]
                            if targetedGene.targetedTranscripts is None:
                                allTranscripts=[]
                            else:
                                allTranscripts=targetedGene.targetedTranscripts[:]
                            if targetedGene.notTargetedTranscripts is not None:
                                allTranscripts.extend(targetedGene.notTargetedTranscripts)
                            if allTranscripts==[]:
                                aceNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            else:                                                                                                                    
                                aceNotTargetedTranscripts[targetedGeneID][probeNb]=allTranscripts
                        else:
                            if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                                aceTargetedTranscripts[targetedGeneID][probeNb]=[[]];
                            else:                                              
                                aceTargetedTranscripts[targetedGeneID][probeNb]=[list(targetedGene.targetedTranscripts)]                        
                            if targetedGene.notTargetedTranscripts is None:
                                aceNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                            else:                                                                                                                                                
                                aceNotTargetedTranscripts[targetedGeneID][probeNb]=[list(targetedGene.notTargetedTranscripts)]                    
                                
                            exons={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundExons=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currExon=firstStructureIDs[probeIndexe][item]
                                            if foundExons.count(currExon)==0:
                                                foundExons.append(currExon)
                                                if currExon in exons.keys():
                                                    exons[currExon]=exons[currExon]+1
                                                else:
                                                    exons[currExon]=1                                     
                            aceTargetedExons[targetedGeneID][probeNb]=[exons]        
                            
                            groups={};
                            for probeIndexe in range(len(probeIndexes)):
                                if firstStructureGroups[probeIndexe]!=[]:                            
                                    foundGroup=[];
                                    for item in range(len(firstStructureGroups[probeIndexe])):
                                        if probeLocalisations[probeIndexe][item]=='e':
                                            currGroup=firstStructureGroups[probeIndexe][item]
                                            if foundGroup.count(currGroup)==0:
                                                foundGroup.append(currGroup)
                                                if currGroup in groups.keys():
                                                    groups[currGroup]=groups[currGroup]+1
                                                else:
                                                    groups[currGroup]=1                                     
                            aceTargetedGroups[targetedGeneID][probeNb]=[groups]                                    
                                
                        aceInGeneProbeNbs[targetedGeneID]={}                    
                        aceInGeneProbeNbs[targetedGeneID][probeNb]=[inGeneProbeNb]                        
                        aceNotInExonProbeNbs[targetedGeneID]={}                    
                        aceNotInExonProbeNbs[targetedGeneID][probeNb]=[notInExonProbeNb]
                    #existing gene    
                    else:         
                        #the gene exist but the key probeNb is may be not yet used                                               
                        try:                     
                            aceProbesetIDs[targetedGeneID][probeNb].append(probeset.probesetID)                        
                            aceProbesetIndexes[targetedGeneID][probeNb].append(probeset.probesetIndex)
                            if probeNb==0:
                                aceTargetedExons[targetedGeneID][probeNb].append([])
                                aceTargetedGroups[targetedGeneID][probeNb].append([])
                                aceTargetedTranscripts[targetedGeneID][probeNb].append([])
                                if targetedGene.targetedTranscripts is None:
                                    allTranscripts=[]
                                else:
                                    allTranscripts=targetedGene.targetedTranscripts[:]
                                if targetedGene.notTargetedTranscripts is not None:
                                    allTranscripts.extend(targetedGene.notTargetedTranscripts)
                                aceNotTargetedTranscripts[targetedGeneID][probeNb].append(allTranscripts)
                            else:                        
                                if targetedGene.targetedTranscripts is not None:
                                    aceTargetedTranscripts[targetedGeneID][probeNb].append(list(targetedGene.targetedTranscripts))                            
                                if targetedGene.notTargetedTranscripts is not None:    
                                    aceNotTargetedTranscripts[targetedGeneID][probeNb].append(list(targetedGene.notTargetedTranscripts))                                                                            
                                    
                                exons={};
                                for probeIndexe in range(len(probeIndexes)):
                                    if firstStructureGroups[probeIndexe]!=[]:                            
                                        foundExons=[];
                                        for item in range(len(firstStructureGroups[probeIndexe])):
                                            if probeLocalisations[probeIndexe][item]=='e':
                                                currExon=firstStructureIDs[probeIndexe][item]
                                                if foundExons.count(currExon)==0:
                                                    foundExons.append(currExon)
                                                    if currExon in exons.keys():
                                                        exons[currExon]=exons[currExon]+1
                                                    else:
                                                        exons[currExon]=1                                     
                                aceTargetedExons[targetedGeneID][probeNb].append(exons)        
                                    
                                groups={};
                                for probeIndexe in range(len(probeIndexes)):
                                    if firstStructureGroups[probeIndexe]!=[]:                            
                                        foundGroup=[];
                                        for item in range(len(firstStructureGroups[probeIndexe])):
                                            if probeLocalisations[probeIndexe][item]=='e':
                                                currGroup=firstStructureGroups[probeIndexe][item]
                                                if foundGroup.count(currGroup)==0:
                                                    foundGroup.append(currGroup)
                                                    if currGroup in groups.keys():
                                                        groups[currGroup]=groups[currGroup]+1
                                                    else:
                                                        groups[currGroup]=1                                     
                                aceTargetedGroups[targetedGeneID][probeNb].append(groups)                                    
                            aceInGeneProbeNbs[targetedGeneID][probeNb].append(inGeneProbeNb)                        
                            aceNotInExonProbeNbs[targetedGeneID][probeNb].append(notInExonProbeNb)
                        except:
                            #existing gene but new probe nb
                            try:
                                if aceProbesetIDs[targetedGeneID][probeNb]!=[]:
                                    print 'error aceview: '+repr(aceProbesetIDs[targetedGeneID][probeNb])
                                    sys.exit()                            
                            except:
                                pass                        
                            aceProbesetIDs[targetedGeneID][probeNb]=[probeset.probesetID]
                            aceProbesetIndexes[targetedGeneID][probeNb]=[probeset.probesetIndex]
                            if probeNb==0:
                                aceTargetedExons[targetedGeneID][probeNb]=[[]]
                                aceTargetedExons[targetedGeneID][probeNb]=[[]]
                                aceTargetedGroups[targetedGeneID][probeNb]=[[]]
                                aceTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                                if targetedGene.targetedTranscripts is None:
                                    allTranscripts=[]
                                else:
                                    allTranscripts=targetedGene.targetedTranscripts[:]
                                if targetedGene.notTargetedTranscripts is not None:
                                    allTranscripts.extend(targetedGene.notTargetedTranscripts)
                                if allTranscripts==[]:
                                    aceNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                                else:                                                                                                                    
                                    aceNotTargetedTranscripts[targetedGeneID][probeNb]=allTranscripts                                                                              
                            else:
                                if targetedGene.targetedTranscripts is None or targetedGene.targetedTranscripts==[]:
                                        aceTargetedTranscripts[targetedGeneID][probeNb]=[[]]
                                else:
                                    aceTargetedTranscripts[targetedGeneID][probeNb]=[list(targetedGene.targetedTranscripts)]                            
                                
                                if targetedGene.notTargetedTranscripts is None or targetedGene.notTargetedTranscripts==[]:
                                    aceNotTargetedTranscripts[targetedGeneID][probeNb]=[[]]    
                                else:
                                    aceNotTargetedTranscripts[targetedGeneID][probeNb]=[list(targetedGene.notTargetedTranscripts)]
                                    
                                exons={};
                                for probeIndexe in range(len(probeIndexes)):
                                    if firstStructureGroups[probeIndexe]!=[]:                            
                                        foundExons=[];
                                        for item in range(len(firstStructureGroups[probeIndexe])):
                                            if probeLocalisations[probeIndexe][item]=='e':
                                                currExon=firstStructureIDs[probeIndexe][item]
                                                if foundExons.count(currExon)==0:
                                                    foundExons.append(currExon)
                                                    if currExon in exons.keys():
                                                        exons[currExon]=exons[currExon]+1
                                                    else:
                                                        exons[currExon]=1                                     
                                aceTargetedExons[targetedGeneID][probeNb]=[exons]        
                                    
                                groups={};
                                for probeIndexe in range(len(probeIndexes)):
                                    if firstStructureGroups[probeIndexe]!=[]:                            
                                        foundGroup=[];
                                        for item in range(len(firstStructureGroups[probeIndexe])):
                                            if probeLocalisations[probeIndexe][item]=='e':
                                                currGroup=firstStructureGroups[probeIndexe][item]
                                                if foundGroup.count(currGroup)==0:
                                                    foundGroup.append(currGroup)
                                                    if currGroup in groups.keys():
                                                        groups[currGroup]=groups[currGroup]+1
                                                    else:
                                                        groups[currGroup]=1                                     
                                aceTargetedGroups[targetedGeneID][probeNb]=[groups]                                                                                            
                            aceInGeneProbeNbs[targetedGeneID][probeNb]=[inGeneProbeNb]                                                        
                            aceNotInExonProbeNbs[targetedGeneID][probeNb]=[notInExonProbeNb]
    print 'saving data ....'
    t2=time.time()
    ensGeneIDs.sort()                            
    for GeneID in ensGeneIDs:        
        probesets=ps_class.TargetingProbeset(probesetNames=ensProbesetIDs[GeneID],probesetIndexes=ensProbesetIndexes[GeneID],targetedExons=ensTargetedExons[GeneID],targetedGroups=ensTargetedGroups[GeneID],targetedTranscripts=ensTargetedTranscripts[GeneID],notTargetedTranscripts=ensNotTargetedTranscripts[GeneID],inGeneProbeNbs=ensInGeneProbeNbs[GeneID],notInExonProbeNbs=ensNotInExonProbeNbs[GeneID])
        ensemblDB[GeneID]=cPickle.dumps(probesets,protocol=-1)
    ensemblDB.close() 
    probesetDB.close()                   
    if log!=0:
        log.write('%s\t%s\t%s\t27\t%s_probesets_by_gene.bkdb\tensembl\t%u\t%.2f\n'%(date.today(),species,chipName,chipName,len(ensGeneIDs),t2-t1))        
    if aceFlag:    
        aceGeneIDs.sort()                            
        for GeneID in aceGeneIDs:        
            probesets=ps_class.TargetingProbeset(probesetNames=aceProbesetIDs[GeneID],probesetIndexes=aceProbesetIndexes[GeneID],targetedExons=aceTargetedExons[GeneID],targetedGroups=aceTargetedGroups[GeneID],targetedTranscripts=aceTargetedTranscripts[GeneID],notTargetedTranscripts=aceNotTargetedTranscripts[GeneID],inGeneProbeNbs=aceInGeneProbeNbs[GeneID],notInExonProbeNbs=aceNotInExonProbeNbs[GeneID])
            aceviewDB[GeneID]=cPickle.dumps(probesets,protocol=-1)        
        aceviewDB.close()
        if log!=0:
            log.write('%s\t%s\t%s\t27\t%s_probesets_by_gene.bkdb\taceviewl\t%u\t%.2f\n'%(date.today(),species,chipName,chipName,len(aceGeneIDs),t2-t1))       
    
            

   

    

   
if __name__=='__main__':        
    erase_gop('mouse','m8')
#    position_gopmapping('mouse','m8','1')
    pass
#    os.chdir(os.environ['PYDATA'])
#    log=open('log.txt','a')       
#    t1=time.time()        
#    assign_probeset('mouse','m5','ensembl',log=log)
#    t2=time.time()
#    print 'processed in %u s'%int(round(t2-t1))        
#    print '%s\n'%('*'*5)
#    log.write('%s\t%s\t\t23\t%s_probeset.bkdb\t\tensembl update\t%.2f\n'%(date.today(),'mouse','mouse',t2-t1))

#    write_probeset_list('mouse','m8')

#    assign_probeset('mouse','m8','ensembl',psRange=[9])
    
#    print '*'*5+' 28: WRITE PROBESETS BY GENE for '+'m8' 
#    print '\n'        
#    write_probesets_by_gene('mouse','m8',log=log)
#    print '%s\n'%('*'*5)
#    
#
#    print '*'*5+' 29: WRITE PROBESET LIST for '+'m8'
#    print '\n'    
#    write_probeset_list('mouse','m8',11,log=log)
#    print '*'*5 
#    print '%s\n'%('*'*5)
#
#
#    print '*'*5+' 30: WRITE ACEVIEW GENES for '+'mouse'
#    print '\n'            
#    write_aceview_genes('mouse',log=log)        
#    print '%s\n'%('*'*5)        