#!/usr/bin/python
"""**Functions for processing positions.**"""
import os, sys, cPickle, time
from datetime import date
from numpy import array,int8,uint8,uint32,append,zeros, where
import bsddb
from bsddb import db
try:
    import MySQLdb
except:
    pass
from dbtools import makebt
import ps_class

def positions_by_region(species,myChipList,ensVersion,ensTable,ensChipList,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find positions by Ensembl region for a list of chips.**    
    ::
    
        Arguments:
            species ----- species name
            myChipList -- list of chip names
            ensVersion -- Ensembl version
            ensTable ---- Ensembl table to be interrogated by MySQL                
            host ------- Ensembl database address    
            port ------- Ensembl port
            user ------- user name
            pswd ------- user pasword
            log -------- handle of a log file for recording messages
                
        Output:    
            Create %chip_positions_by_region.bkdb
                key ==== Ensembl region id
                value == PositionList object
                                    
    """
                     
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_chromosome_by_region.bkdb"            
    regionDB=bsddb.btopen(path)
    regions=regionDB.keys()
    regionDB.close()       
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 18 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 18 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)
    
    # Create a mySQL cursor         
    cursor = connection.cursor()    
    failed=0     
    for ensChipName in ensChipList:
        print 'processing chip %s'%ensChipName
        if failed==0:
            t1=time.time()
            chipRank=ensChipList.index(ensChipName)                                
            myChipName=myChipList[chipRank]
            
            # Initialize the object dictionary used to create the berkeley database
            dict={}
            # Put each resulting lines in dictionary
            region=''
#            dict={}    
#            positions=[]
            positionNb=0
            
            #information arranged in array
            arrayDict={};                        
        
            for region in regions:
                if failed==0:
                    print 'processing region %s'%region
                    #find seq_region_id
                    if int(ensVersion.split("_")[0])>=48:            
                        query="""
                        SELECT seq_region_id
                        FROM seq_region
                        WHERE core_seq_region_id="%s"
                        AND schema_build="%s"
                        """%(region,ensVersion)
                        cursor.execute (query)
                        lines=cursor.fetchall()
                        try:
                            seqRegionId=str(lines[0][0])
                            isRegion=1;
                        except:
                            isRegion=0
                    if isRegion:
                        t2=time.time()
                        if int(ensVersion.split("_")[0])<48:
                            query = """
                            SELECT DISTINCT  
                            oligo_probe.name,                                                                     
                            oligo_feature.seq_region_start,
                            oligo_feature.seq_region_end, 
                            oligo_feature.seq_region_strand,
                            oligo_feature.mismatches 
                            FROM 
                            oligo_probe, oligo_array, oligo_feature 
                            WHERE
                            oligo_array.name='%s'         
                            AND oligo_feature.seq_region_id='%s'
                            AND oligo_probe.oligo_probe_id=oligo_feature.oligo_probe_id                
                            AND oligo_probe.oligo_array_id=oligo_array.oligo_array_id                   
             
                            ORDER BY 
                            oligo_probe.oligo_probe_id, oligo_feature.seq_region_start
                            """ %(region,ensChipName)                        
                            t1=time.time()                            
                            tryQuery=1                                      
                            while tryQuery:
                                try:    
                                    cursor.execute (query)
                                    tryQuery=0
                                    failed=0
                                except:
                                    print 'lost MySQL connection for region %s'%region
                                    tryQuery=raw_input('would you like to continue (enter 1) or to cancel (enter 0)')                                   
                                    if tryQuery=='0':
                                        tryQuery=0                     
                                        failed=1
                                        if log !=0:
                                            log.write('%s\t%s\t\t18\t%s_position_by_region.bkdb\tens_%s\tCANCELED FROM REGION %s\t\n'%(date.today(),species,myChipName,ensVersion,region))
                                        sys.exit()                        
                                    elif tryQuery!='1':
                                        print 'enter 1 or 0'
                            cursor.execute (query)
                        else:
                            query = """
                            SELECT DISTINCT    
                            probe.name,                                                                             
                            probe_feature.seq_region_start,
                            probe_feature.seq_region_end, 
                            probe_feature.seq_region_strand,
                            probe_feature.mismatches 
                            FROM 
                            probe, array_chip, probe_feature 
                            WHERE         
                            array_chip.name='%s'
                            AND probe_feature.seq_region_id="%s"                
                            AND probe.array_chip_id=array_chip.array_chip_id                                                
                            AND probe_feature.probe_id =probe.probe_id
                            """ %(ensChipName,seqRegionId)
                            #print query                
                            t1=time.time()  
                            tryQuery=1                                      
                            while tryQuery:
                                try:    
                                    cursor.execute (query)
                                    tryQuery=0
                                    failed=0
                                except:
                                    print 'lost MySQL connection for region %s'%region
                                    tryQuery=raw_input('would you like to continue (enter 1) or to cancel (enter 0)')                                   
                                    if tryQuery=='0':
                                        tryQuery=0                     
                                        failed=1
                                        if log !=0:
                                            log.write('%s\t%s\t\t18\t%s_position_by_region.bkdb\tens_%s\tCANCELED FROM REGION %s\t\n'%(date.today(),species,myChipName,ensVersion,region))
                                        sys.exit()                        
                                    elif tryQuery!='1':
                                        print 'enter 1 or 0'    
                        if tryQuery==0 and failed==0:    
                            t3=time.time()
                            lines = cursor.fetchall()
                            print "mySQL database %s delivered %u results in %.2f s for region %s for chip %s\n"%("ensembl_"+ensTable,len(lines),t3-t2,region,myChipName)                                                           
                            positionNb=positionNb+len(lines)    
                            
                                
                          
                            probeIDArray=array((),dtype='str');
                            medianArray=array((),dtype=uint32);
                            strandArray=array((),dtype=int8);
                            mismatchNbArray=array((),dtype=uint8);
    #                        positions=[]
                                                    
                            
                            for n in lines:                                                    
                                probeID=n[0]
                                if probeID[-1]==';':
                                    probeID=probeID[:-1]
                                if ':' in probeID:    
                                    probeIDs=probeID.split(":")
                                    probeID="%04u:%04u"%(int(probeIDs[0]),int(probeIDs[1]))                                                                                                                                                                                                                                                                                                
                                probeIDArray=append(probeIDArray,probeID)
                                medianArray=append(medianArray,uint32(n[2]-(n[2]-n[1])/2))
                                strandArray=append(strandArray,int8(n[3]))
                                mismatchNbArray=append(mismatchNbArray,uint8(n[4]))
                                                    
                            
                            #save last region                      
                            #sort data on median position
                            sortIndex=medianArray.argsort() 
                            probeIDArray=probeIDArray[sortIndex];
                            medianArray= medianArray[sortIndex];
                            strandArray=strandArray[sortIndex];
                            mismatchNbArray=mismatchNbArray[sortIndex];        
                            arrayDict[region]=ps_class.PositionList(probeIDList=list(probeIDArray),medianArray=medianArray,strandArray=strandArray,mismatchNbArray=mismatchNbArray)                                                        
        # create berkeley database
        if failed==0:
            print 'TOTAL NB OF POSITIONS: %u\n'%positionNb
            if log!=0:
                t4=time.time() 
                log.write('%s\t%s\t%s\t18\t%s_position_by_region.bkdb\tens_%s\t%u\t%.2f\n'%(date.today(),species,myChipName,myChipName,ensVersion,positionNb,t4-t1))                                                   
            path=os.environ['PYDATA']+"/"+species+"/ensembl/"+myChipName+"_positions_by_region.bkdb"
            if os.path.exists(path):
                os.remove(path)
            makebt(dict=arrayDict,path=path,newFlag=1,objFlag=1)                  
    cursor.close ()
    connection.close ()
    
                      
def merge_positions(species,chipList,region,log=0):    
    """
    **Merge positions by ensembl region for a species.**
    ::
    
        creates one database merging the probe position of all the chip sets of a species
        +1 (Plus) and -1 (Minus) strands position are separated in order to have position ordered inside
        each category
            
        Arguments:
            species --- species name
            chipList -- list of chips
            log -------- handle of a log file for recording messages
            
        Input:    
            %chip_positions_by_region.bkdb
            
        Output:
            Create %species_positions_by_region.bkdb
                key ==== Ensembl region id
                value == PositionList object
            
    """
    
    t1=time.time()
    #path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_positions_by_region.bkdb"
    #WARNING : process all the chips at once (this version does not allow incremental feeding of the position data base)
    #if os.path.exists(path):
    #    os.remove(path)            
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+'%s_positions_by_region.bkdb'%species
    mergedPosDB=bsddb.btopen(path)  
    try:     
        a=cPickle.loads(mergedPosDB[region])
        print '\nLOADED REGION %s\n'%region                    
        mergedPosDB.close()        
    except:
        mergedPosDB.close()
        print '\nPROCESSING REGION %s\n'%region                                                
        #calculate the total number of referenced positions
        positionNb=0
        posByChip={}
        keptByChip={}        
        for chipName in chipList:    
            path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_positions_by_region.bkdb"        
            positionDB=bsddb.btopen(path)        
            posByChip[chipName]=0
            keptByChip[chipName]=0
            #for region in positionDB.keys():
            try:
                positions=cPickle.loads(positionDB[region])
                posNb=len(positions.medians)            
                posByChip[chipName]=posByChip[chipName]+posNb                
                positionNb+=posNb
            except:
                pass
    #            if region in positionNb.keys():            
    #                positionNb[region]=positionNb[region]+posNb
    #            else:
    #                positionNb[region]=posNb
            positionDB.close()
            print 'chip %s: %u positions to be scanned'%(chipName,posByChip[chipName])
        print '\nTOTAL NB OF POSITIONS TO BE SCANNED: %u\n'%positionNb    
        
        if positionNb>0:
                #process current region
                #for region in positionNb.keys():
                print 'processing region %s'%region
                #create empty variable for positions on plus strand                                     
                chipsPlus=[]  
                probeIDsPlus=[]                             
                mediansPlus=zeros(positionNb,dtype=uint32)
                strandsPlus=zeros(positionNb,dtype=int8)
                mismatchesPlus=zeros(positionNb,dtype=uint8)   
                
                #create empty variable for positions on minus strand        
                chipsMinus=[]   
                probeIDsMinus=[]
                mediansMinus=zeros(positionNb,dtype=uint32)
                strandsMinus=zeros(positionNb,dtype=int8)
                mismatchesMinus=zeros(positionNb,dtype=uint8)
                    
                #fill the variables        
                plusIndex=0
                minusIndex=0
                for chipName in chipList:  
                    chipRank= int(chipName[1:]) 
                    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_positions_by_region.bkdb"        
                    positionDB=bsddb.btopen(path)
                    positions=cPickle.loads(positionDB[region])                                                                    
                    for posRank in range(len(positions.medians)):                                
                        if positions.strands[posRank]==1: 
                            if positions.medians[posRank] not in mediansPlus:
                                keptByChip[chipName]+=1                                        
                                chipsPlus.append(set([chipRank]))                                               
                                probeIDsPlus.append(set([positions.probeIDs[posRank]]))                                                            
                                mediansPlus[plusIndex]=positions.medians[posRank]                                
                                strandsPlus[plusIndex]=positions.strands[posRank]
                                mismatchesPlus[plusIndex]=positions.mismatches[posRank]
                                plusIndex+=1
                            else:                        
                                plusPos=where(mediansPlus==positions.medians[posRank])[0][0]
                                chipsPlus[plusPos].add(chipRank)                         
                                probeIDsPlus[plusPos].add(positions.probeIDs[posRank])
                        else:
                            if positions.medians[posRank] not in mediansMinus:
                                keptByChip[chipName]+=1                                        
                                chipsMinus.append(set([chipRank]))                 
                                probeIDsMinus.append(set([positions.probeIDs[posRank]]))                              
                                mediansMinus[minusIndex]=positions.medians[posRank]
                                strandsMinus[minusIndex]=positions.strands[posRank]
                                mismatchesMinus[minusIndex]=positions.mismatches[posRank]
                                minusIndex+=1
                            else:
                                minusPos=where(mediansMinus==positions.medians[posRank])[0][0]                    
                                chipsMinus[minusPos].add(chipRank)                                                
                                probeIDsMinus[minusPos].add(positions.probeIDs[posRank])              
                    positionDB.close()            
                    
                #process plus strand positions
                #remove zeros positions
                                
                mediansPlus=mediansPlus[0:plusIndex]
                strandsPlus=strandsPlus[0:plusIndex]
                mismatchesPlus=mismatchesPlus[0:plusIndex]       
                #reorder data on median position
                sortIndex=mediansPlus.argsort()
                chipsPlusList=[]
                probeIDPlusList=[]
                for i in sortIndex:        
                    chipsPlusList.append(chipsPlus[i])        
                    probeIDPlusList.append(probeIDsPlus[i])
                chipsPlus=[]                   
                mediansPlus=mediansPlus[sortIndex]
                strandsPlus=strandsPlus[sortIndex]
                mismatchesPlus=mismatchesPlus[sortIndex]        
                
                #process minus strand positions
                #remove zeros positions                    
                mediansMinus=mediansMinus[0:minusIndex]
                strandsMinus=strandsMinus[0:minusIndex]
                mismatchesMinus=mismatchesMinus[0:minusIndex]        
                #reorder data on median position
                sortIndex=mediansMinus.argsort()
                chipsMinusList=[]
                probeIDMinusList=[]
                for i in sortIndex: 
                    chipsMinusList.append(chipsMinus[i])
                    probeIDMinusList.append(probeIDsMinus[i])
                chipsMinus=[]                         
                mediansMinus=mediansMinus[sortIndex]
                strandsMinus=strandsMinus[sortIndex]
                mismatchesMinus=mismatchesMinus[sortIndex]                                                           
                
                #merge Plus and Minus (Minus first)                    
                print 'plus pos: %u\tminus pos: %u'%(len(chipsPlusList),len(chipsMinusList))
                chipsMinusList.extend(chipsPlusList)
                probeIDMinusList.extend(probeIDPlusList)        
                mediansMinus=append(mediansMinus,mediansPlus)
                strandsMinus=append(strandsMinus,strandsPlus)
                mismatchesMinus=append(mismatchesMinus,mismatchesPlus)
                
                path=os.environ['PYDATA']+"/"+species+"/ensembl"    
                # Get an instance of BerkeleyDB environment
                dbEnv = db.DBEnv()
                # Set the automatic deadlock detection mode.
                # The detect parameter configures the deadlock detector. The deadlock detector
                # will reject the lock request with the lowest priority. If multiple lock 
                # requests have the lowest priority, then the detect parameter is used 
                # to select which of those lock requests to reject
                # DB_LOCK_YOUNGEST : Reject the lock request for the locker ID with the youngest lock
                dbEnv.set_lk_detect(db.DB_LOCK_YOUNGEST)               
                dbEnv.open(path,db.DB_INIT_LOCK|db.DB_CREATE|db.DB_INIT_MPOOL)    
                #Aquires a locker id, guaranteed to be unique across all threads and processes that have the DBEnv open
                lockID = dbEnv.lock_id()    
                dbName='%s_positions_by_region.bkdb'%species
                mergedPosDB=db.DB(dbEnv)
                mergedPosDB.open(dbName,None,db.DB_BTREE,db.DB_CREATE)
                    
                #record region
                print '\nASK FOR A LOCK ALLOWING DUMPING REGION %s\n'%region
                while(1):
                    try:
                        #Aquires a lock and returns a handle to it as a DBLock object.
                        #The  first parametr is the locker parameter is an integer representing the entity doing the locking (here lockID)
                        #The second parameter is a string representing the item to be locked
                        lock = dbEnv.lock_get(lockID, "anytid", db.DB_LOCK_WRITE,db.DB_LOCK_NOWAIT)
                        #flags:                 
                        #   DB_LOCK_NOWAIT:If a lock cannot be granted because the requested lock conflicts 
                        #                   with an existing lock, return immediately instead of waiting for
                        #                   the lock to become available.
                        #mode:
                        #   DB_LOCK_WRITE: write (exclusive)
                        break
                    except Exception:                
                        pass
            
                
                mergedPosDB.put(region,cPickle.dumps(ps_class.PositionList(probeIDList=probeIDMinusList,chipsList=chipsMinusList,medianArray=mediansMinus,strandArray=strandsMinus,mismatchNbArray=mismatchesMinus),protocol=-1))
                print '\nDUMPED REGION %s\n'%region                
                #mergedPosDB[region]=cPickle.dumps(ps_class.PositionList(probeIDList=probeIDMinusList,chipsList=chipsMinusList,medianArray=mediansMinus,strandArray=strandsMinus,mismatchNbArray=mismatchesMinus),protocol=-1)
            #    try:
            #        a=cPickle.loads(mergedPosDB[region])
            #        print '\nLOADED REGION %s\n'%region
            #    except:
            #        print '\nFAILED TO LOAD REGION %s\n'%region
            
                #Release the lock
                dbEnv.lock_put(lock)        
                mergedPosDB.close()
                dbEnv.close()
                if log!=0:
                    t2=time.time()
                    for chipName in posByChip.keys():
                        log.write('%s\t%s\t%s\t19\t%s_positions_by_region(%s).bkdb\t%u\t%u\t%.2f\n'%(date.today(),species,chipName,species,region,posByChip[chipName],keptByChip[chipName],t2-t1))                
                
def gop_by_region(species,aceFlag,log=0):
    """
    **Construct groups of probes (GOPs)**
    ::
    
        Probes that are less than 2 kb apart and that are outside gene limits +/- 2 kb 
        are grouped and given an Ensmelb  gene-like name (GOPMUSG00000000001 ...)
    
        Arguments:
            species -- species name
            log -------- handle of a log file for recording messages
    
        Input:
            %species_genes_by_region.bkdb (Ensembl genes positions)
            %species_genes_by_ensembl_region.bkdb (AceView genes positions)
            %species_positions_by_region.bkdb
            %species_chromosomes_by_region.bkdb
    
        Output:
            create %chip_gops_by_region.bkdb (chips used are those used to search positions)
                       
    """    
    
    if aceFlag:    
        path=os.environ['PYDATA']+'/'+species+'/aceview/'+species+'_genes_by_ensembl_region.bkdb'   
        aceGeneDB=bsddb.btopen(path)
        
    path=os.environ['PYDATA']+"/"+species+'/ensembl/'+species+'_genes_by_region.bkdb'    
    ensGeneDB=bsddb.btopen(path)
                    
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_positions_by_region.bkdb"    
    positionDB=bsddb.btopen(path)
        
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_chromosome_by_region.bkdb"                        
    regionDB=bsddb.btopen(path)
    regions=regionDB.keys()
    
        
    gopRank=0
    
    for region in regions:  
        print 'region %s'%region
        t1=time.time()  
        dict={}              
        chips=[]
        gopIndex=-1
        
        
        if region not in ensGeneDB.keys():
            print 'region %s (%s)  not in ensGeneDB'%(region,regionDB[region])
        else:
            if region not in positionDB.keys():
                print 'region %s (%s) not in positionDB'%(region,regionDB[region])
            else:    
                print 'PROCESSING REGION %s (CHROMOSOME %s)'%(region,regionDB[region])          
                  
                genes=cPickle.loads(ensGeneDB[region])
                if aceFlag:
                    #some regions do not exist in acegenes _by_region
                    try:
                        aceGenes=cPickle.loads(aceGeneDB[region])
                        #concatenate arrays
                        genes.IDs=append(genes.IDs,aceGenes.IDs)
                        genes.starts=append(genes.starts,aceGenes.starts)
                        genes.ends=append(genes.ends,aceGenes.ends)
                        genes.strands=append(genes.strands,aceGenes.strands)        
                    except:
                        pass
                  
                               
                #order genes on start position
                sortIndex=genes.ends.argsort()      
                genes.starts=genes.starts[sortIndex]  
                genes.ends=genes.ends[sortIndex]
                genes.strands=genes.strands[sortIndex]
                genes.IDs=[genes.IDs[i] for i in sortIndex]        
                        
                sortIndex=genes.starts.argsort()        
                genes.starts=genes.starts[sortIndex]
                genes.ends=genes.ends[sortIndex]
                genes.strands=genes.strands[sortIndex]        
                genes.IDs=[genes.IDs[i] for i in sortIndex]
                
                #add 2kb in 5' and 3' end
                genes.starts=genes.starts-2000
                genes.ends=genes.ends+2000
                
                #separate plus and minus strand genes
                plusIndex=where(genes.strands==1)
                genesPlus=ps_class.GeneList()        
                genesPlus.IDs=[genes.IDs[i] for i in plusIndex[0]]            
                genesPlus.starts=genes.starts[plusIndex]
                genesPlus.ends=genes.ends[plusIndex]
                genesPlus.strands=genes.strands[plusIndex]
                
                minusIndex=where(genes.strands==-1)
                genesMinus=ps_class.GeneList()
                genesMinus.IDs=[genes.IDs[i] for i in minusIndex[0]]
                genesMinus.starts=genes.starts[minusIndex]
                genesMinus.ends=genes.ends[minusIndex]
                genesMinus.strands=genes.strands[minusIndex]
                    
                '''recover the list of probe positions for the current region'''                                      
                posArray=cPickle.loads(positionDB[region])                
                          
                plusNb=0
                minusNb=0
                
                
                
                gopIDArray=array(())
                startArray=array(())
                endArray=array(())
                strandArray=array(())
                upGeneDistanceArray=array(())                
                downGeneDistanceArray=array(())        
                upGeneIDArray=array(())      
                downGeneIDArray=array(())   
                positionNbArray=array(())   
                gopsByChip={}     
                
                
                #localization of the first position of plus strand        
                plusIndex=posArray.strands.searchsorted(1)
                insideNb=0    
                for  strand in [-1,1]:
                    if strand==-1:
                        scanFlag=1
                        if plusIndex>1:                                                            
                            firstIndex=0
                            lastIndex=plusIndex
                            genes=genesMinus                                    
                        else:
                            scanFlag=0
                        
                    else:
                        scanFlag=1                        
                        if plusIndex<len(posArray.medians):                                                                                                
                            firstIndex=plusIndex
                            lastIndex=len(posArray.medians)
                            genes=genesPlus                                                        
                        else:
                            scanFlag=0                    
                    chips=set()       
                                                            
                    if scanFlag:     
                        positionNb=0                        
                        geneNb=len(genes.IDs)
                        #rank of left gene                    
                        leftGeneRank=-1
                        #test if the first position is inside or outside a gene
                        if posArray.medians[firstIndex]<genes.starts[0]:
                            #first position is on the left side of the first gene
                            # => start a gop
                            '''start a new GOP'''
                            gopFlag=1                
                            chips=set(posArray.chips[firstIndex])
                            start=posArray.medians[firstIndex]
                            end=posArray.medians[firstIndex]                                
                            positionNb=1
                            leftGeneRank=-1    
                            insideNb+=1                       
                        elif posArray.medians[firstIndex]<=genes.ends[0]:
                            #first position is inside first gene
                            gopFlag=0
                            leftGeneRank=0                
                        else:                                                              
                            #first position is on the right side of the first gene
                            leftGeneRank=0
                            while genes.ends[leftGeneRank]<posArray.medians[firstIndex] and leftGeneRank<geneNb-1:
                                leftGeneRank+=1
                            if leftGeneRank<geneNb-1:    
                                if genes.starts[leftGeneRank]<=posArray.medians[firstIndex]:
                                    #first position is inside leftGeneRank gene
                                    gopFlag=0                
                                else:
                                    #first position is on the left side of the leftGeneRank gene
                                    # => start a gop
                                    '''start a new GOP'''
                                    gopFlag=1                
                                    chips=set(posArray.chips[firstIndex])
                                    start=posArray.medians[firstIndex]
                                    end=posArray.medians[firstIndex]                                
                                    positionNb=1
                                    leftGeneRank=leftGeneRank-1  
                                    insideNb+=1                      
                            else:
                                if posArray.medians[firstIndex]<genes.starts[leftGeneRank]:
                                    #first position is on the left side of the last gene
                                    # => start a gop
                                    '''start a new GOP'''
                                    gopFlag=1                
                                    chips=set(posArray.chips[firstIndex])
                                    start=posArray.medians[firstIndex]
                                    end=posArray.medians[firstIndex]                                
                                    positionNb=1
                                    leftGeneRank=geneNb-2
                                    insideNb+=1                        
                                elif posArray.medians[firstIndex]<=genes.ends[0]:
                                    #first position is inside the last gene
                                    gopFlag=0
                                    leftGeneRank=geneNb-1
                                else:
                                    #first position is on right side of the last gene
                                    # => start a gop
                                    '''start a new GOP'''
                                    gopFlag=1                
                                    chips=set(posArray.chips[firstIndex])
                                    start=posArray.medians[firstIndex]
                                    end=posArray.medians[firstIndex]                                
                                    positionNb=1
                                    leftGeneRank=geneNb-1
                                    insideNb+=1
                        if leftGeneRank==-1:
                            leftGeneRank=0            
                            leftGeneID=''
                            leftGeneEnd=-1  
                            rightGeneRank=1              
                            rightGeneID=genes.IDs[rightGeneRank]
                            rightGeneStart=genes.starts[rightGeneRank]                
                        elif leftGeneRank==geneNb-1:                
                            leftGeneID=genes.IDs[leftGeneRank]
                            leftGeneEnd=genes.ends[leftGeneRank]
                            rightGeneID=''
                            rightGeneStart=posArray.medians[lastIndex-1]+1                    
                            rightGeneRank=geneNb-1            
                        else:
                            leftGeneID=genes.IDs[leftGeneRank]
                            leftGeneEnd=genes.ends[leftGeneRank]
                            #find the next non overlaping gene
                            rightGeneRank=leftGeneRank+1
                            while genes.starts[rightGeneRank]<=leftGeneEnd and rightGeneRank<geneNb-1:
                                leftGeneEnd=max(leftGeneEnd,genes.ends[rightGeneRank])
                                rightGeneRank+=1
                            if rightGeneRank<geneNb-1:
                                rightGeneID=genes.IDs[rightGeneRank]
                                rightGeneStart=genes.starts[rightGeneRank]                    
                            else:
                                rightGeneRank= geneNb-1
                                if genes.starts[rightGeneRank]<=leftGeneEnd:
                                    leftGeneEnd=max(leftGeneEnd,genes.ends[rightGeneRank])
                                    rightGeneID=''
                                    rightGeneStart=posArray.medians[lastIndex-1]+1                                         
                                else:
                                    rightGeneID=genes.IDs[rightGeneRank]
                                    rightGeneStart=genes.starts[rightGeneRank]                                                
                                
                                                                                                                                                                            
                        if firstIndex+1<lastIndex:                                         
                            for index in range(firstIndex+1,lastIndex):                               
                                #test if the position is inside or outside a gene
                                #print 'left %u right %u index %u insideNb %u gop %u gopRank %u inInterval %s positonNb %u end %u pos %u start %u'%(leftGeneRank,rightGeneRank,index,insideNb,gopFlag,gopRank,posArray.medians[index]>leftGeneEnd and posArray.medians[index]<rightGeneStart,positionNb,leftGeneEnd,posArray.medians[index],rightGeneStart)
                                if posArray.medians[index]>leftGeneEnd and posArray.medians[index]<rightGeneStart:                         
                                    if gopFlag:
                                        if posArray.medians[index]>end+2000 or index==lastIndex-1:
                                            '''record the current GOP'''
                                            gopRank+=1
                                            gopIndex+=1
                                            if strand==-1:
                                                minusNb+=1
                                            else:         
                                                plusNb+=1                                       
                                            gopIDArray=append(gopIDArray,'GOPMUSG%011u'%gopRank)                        
                                            startArray=append(startArray,start)
                                            if index==lastIndex:
                                                endArray=append(endArray,posArray.medians[index])
                                            else:    
                                                endArray=append(endArray,end)
                                            strandArray=append(strandArray,strand)                        
                                            positionNbArray=append(positionNbArray,positionNb)
                                            for chip in chips:
                                                try:
                                                    gopsByChip['m%u'%chip]=append(gopsByChip['m%u'%chip],gopIndex)
                                                except:
                                                    gopsByChip['m%u'%chip]=array([gopIndex])                                    
                                                
                                            #record the leftGeneID as left genes of the current GOP
                                            if strand==-1:
                                                downGeneIDArray=append(downGeneIDArray,leftGeneID)
                                                downGeneDistanceArray=append(downGeneDistanceArray,start-leftGeneEnd)
                                            else:
                                                upGeneIDArray=append(upGeneIDArray,leftGeneID)
                                                upGeneDistanceArray=append(upGeneDistanceArray,start-leftGeneEnd)
                                                
                                            #record the next contiguous GOP as the right gene of the current GOP                                                                                                                
                                            if strand==-1:
                                                upGeneIDArray=append(upGeneIDArray,'GOPMUSG%011u'%(gopRank+1))
                                                upGeneDistanceArray=append(upGeneDistanceArray,posArray.medians[index]-end)
                                            else:
                                                downGeneIDArray=append(downGeneIDArray,'GOPMUSG%011u'%(gopRank+1))
                                                downGeneDistanceArray=append(downGeneDistanceArray,posArray.medians[index]-end)        
#                                            if len(upGeneIDArray)<len(startArray):
#                                                'stop'                                                                        
                                            #keep in memory the current GOP as the left gene for the next GOP which is contiguous                                        
                                            leftGeneID='GOPMUSG%011u'%gopRank
                                            leftGeneEnd=end
                                            #the right gene is left unchanged
                                            
                                            if index<lastIndex-1:
                                                '''start a new GOP'''
                                                gopFlag=1                
                                                chips=set(posArray.chips[index])
                                                start=posArray.medians[index]
                                                end=posArray.medians[index]                                    
                                                positionNb=1
                                                insideNb+=1                            
                                                                                                                
                                        else:
                                            #continue to scan position for the current GOP
                                            end=posArray.medians[index] 
                                            positionNb+=1
                                            chips=chips.union(set(posArray.chips[index]))
                                            insideNb+=1                                                                                                                                                               
                                    else:                                
                                        '''start a new GOP'''
                                        gopFlag=1                
                                        chips=set(posArray.chips[index])
                                        start=posArray.medians[index]
                                        end=posArray.medians[index]                                
                                        positionNb=1
                                        insideNb+=1
                                        #the right and left genes are left unchanged 
                                        if index==lastIndex-1:
                                            '''record the current GOP'''
                                            gopRank+=1
                                            gopIndex+=1
                                            if strand==-1:
                                                minusNb+=1
                                            else:         
                                                plusNb+=1                                       
                                            gopIDArray=append(gopIDArray,'GOPMUSG%011u'%gopRank)                        
                                            startArray=append(startArray,start)
                                            endArray=append(endArray,end)
                                            strandArray=append(strandArray,strand)                        
                                            positionNbArray=append(positionNbArray,positionNb)
                                            for chip in chips:
                                                try:
                                                    gopsByChip['m%u'%chip]=append(gopsByChip['m%u'%chip],gopIndex)
                                                except:
                                                    gopsByChip['m%u'%chip]=array([gopIndex])                                    
                                                
                                            #record the leftGeneID as left genes of the current GOP
                                            if strand==-1:
                                                downGeneIDArray=append(downGeneIDArray,leftGeneID)
                                                downGeneDistanceArray=append(downGeneDistanceArray,start-leftGeneEnd)
                                            else:
                                                upGeneIDArray=append(upGeneIDArray,leftGeneID)
                                                upGeneDistanceArray=append(upGeneDistanceArray,start-leftGeneEnd)
                                            #record the rightgeneID as right gene of the current GOP                                                                                                                
                                            if strand==-1:
                                                upGeneIDArray=append(upGeneIDArray,rightGeneID)
                                                upGeneDistanceArray=append(upGeneDistanceArray,rightGeneStart-end)
                                            else:
                                                downGeneIDArray=append(downGeneIDArray,rightGeneID)
                                                downGeneDistanceArray=append(downGeneDistanceArray,rightGeneStart-end)
#                                            if len(upGeneIDArray)<len(startArray):
#                                                'stop'    
                                                                                                                                                   
                                else:
                                    
                                    if gopFlag: 
                                        gopFlag=0                       
                                        #quit the current intergenic interval                                                            
                                        #=> record the current GOP
                                        gopRank+=1
                                        gopIndex+=1
                                        if strand==-1:
                                            minusNb+=1
                                        else:         
                                            plusNb+=1                                       
                                        gopIDArray=append(gopIDArray,'GOPMUSG%011u'%gopRank)                        
                                        startArray=append(startArray,start)
                                        endArray=append(endArray,end)
                                        strandArray=append(strandArray,strand)                        
                                        positionNbArray=append(positionNbArray,positionNb)
                                        for chip in chips:
                                            try:
                                                gopsByChip['m%u'%chip]=append(gopsByChip['m%u'%chip],gopIndex)
                                            except:                                     
                                                gopsByChip['m%u'%chip]=array([gopIndex])                                                                                                
                                        #record the leftGeneID as left gene of the current GOP
                                        if strand==-1:
                                            downGeneIDArray=append(downGeneIDArray,leftGeneID)
                                            downGeneDistanceArray=append(downGeneDistanceArray,start-leftGeneEnd)
                                        else:
                                            upGeneIDArray=append(upGeneIDArray,leftGeneID)
                                            upGeneDistanceArray=append(upGeneDistanceArray,start-leftGeneEnd)                                
                                        #record the rightgeneID as right gene of the current GOP                                                                                                                
                                        if strand==-1:
                                            upGeneIDArray=append(upGeneIDArray,rightGeneID)
                                            upGeneDistanceArray=append(upGeneDistanceArray,rightGeneStart-end)
                                        else:
                                            downGeneIDArray=append(downGeneIDArray,rightGeneID)
                                            downGeneDistanceArray=append(downGeneDistanceArray,rightGeneStart-end)                                 
                                                                                                                                                            
                 
                                    #find where is the position (inside or outside a gene ?)
                                    #find the left gene
                                    while genes.ends[leftGeneRank]<posArray.medians[index] and leftGeneRank<geneNb-1:
                                        leftGeneRank+=1
                                    if leftGeneRank<geneNb-1:    
                                        if genes.starts[leftGeneRank]<=posArray.medians[index]:
                                            #current position is inside leftGeneRank gene
                                            gopFlag=0                
                                        else:
                                            #current position is on the left side of the leftGeneRank gene                             
                                            # => start a gop
                                            '''start a new GOP'''
                                            gopFlag=1                
                                            chips=set(posArray.chips[index])
                                            start=posArray.medians[index]
                                            end=posArray.medians[index]                                
                                            positionNb=1
                                            rightGeneRank=leftGeneRank
                                            #=>change leftGeneRank
                                            leftGeneRank=leftGeneRank-1                                
                                            insideNb+=1                                                                                        
                                    else:
                                        rightGeneRank=leftGeneRank
                                        if posArray.medians[index]<genes.starts[leftGeneRank]:
                                            #current position is on the left side of the last gene
                                            #=>change leftGeneRank
                                            leftGeneRank=leftGeneRank-1
                                            # => start a gop
                                            '''start a new GOP'''
                                            gopFlag=1                
                                            chips=set(posArray.chips[index])
                                            start=posArray.medians[index]
                                            end=posArray.medians[index]                                
                                            positionNb=1                                    
                                            insideNb+=1                        
                                        elif posArray.medians[index]<=genes.ends[leftGeneRank]:
                                            #current position is inside the last gene
                                            gopFlag=0                                    
                                        else:
                                            #current position is on right side of the last gene
                                            # => start a gop
                                            '''start a new GOP'''
                                            gopFlag=1                
                                            chips=set(posArray.chips[index])
                                            start=posArray.medians[index]
                                            end=posArray.medians[index]                                
                                            positionNb=1                                    
                                            insideNb+=1                                        
                                                                                                                                        
                                    if leftGeneRank==geneNb-1:                
                                        leftGeneID=genes.IDs[leftGeneRank]
                                        leftGeneEnd=genes.ends[leftGeneRank]
                                        rightGeneID=''
                                        rightGeneStart=posArray.medians[lastIndex-1]+1
                                        rightGeneRank=geneNb-1
                                    else:
                                        leftGeneID=genes.IDs[leftGeneRank]
                                        leftGeneEnd=genes.ends[leftGeneRank]
                                        #find the next non overlaping gene
                                        rightGeneRank=leftGeneRank+1
                                        while genes.starts[rightGeneRank]<=leftGeneEnd and rightGeneRank<geneNb-1:
                                            leftGeneEnd=max(leftGeneEnd,genes.ends[rightGeneRank])
                                            rightGeneRank+=1
                                        if rightGeneRank<geneNb-1:
                                            rightGeneID=genes.IDs[rightGeneRank]
                                            rightGeneStart=genes.starts[rightGeneRank]                                
                                        else:
                                            rightGeneRank=geneNb-1
                                            if leftGeneRank<geneNb-1:
                                                rightGeneRank=  geneNb-1
                                                if genes.starts[rightGeneRank]<=leftGeneEnd:
                                                    leftGeneEnd=max(leftGeneEnd,genes.ends[rightGeneRank])
                                                    rightGeneID=''
                                                    rightGeneStart=posArray.medians[lastIndex-1]+1                                                                                    
                                                else:
                                                    rightGeneID=genes.IDs[rightGeneRank]
                                                    rightGeneStart=genes.starts[rightGeneRank]                                                                                    
                                            else:
                                                rightGeneID=''
                                                rightGeneStart=posArray.medians[lastIndex-1]+1  
                t2=time.time()                                                                                                                                                                                                                                  
                if log!=0:
                    log.write('%s\t%s\tr%s\t20\tgops_by_region.bkdb (%u GOPs)\t%u\t%u\t%.2f\n'%(date.today(),species,region,gopRank,plusNb,minusNb,t2-t1))                                           
                for chip in gopsByChip.keys():
                    dict[chip]={}           
                    dict[chip][region]=ps_class.GopList(idArray=gopIDArray[gopsByChip[chip]],regionArray=array([region for i in range(len(gopsByChip[chip]))]),startArray=uint32(startArray[gopsByChip[chip]]),endArray=uint32(endArray[gopsByChip[chip]]),strandArray=int8(strandArray[gopsByChip[chip]]),upGeneDistanceArray=uint32(upGeneDistanceArray[gopsByChip[chip]]),downGeneDistanceArray=uint32(downGeneDistanceArray[gopsByChip[chip]]),upGeneIDArray=upGeneIDArray[gopsByChip[chip]],downGeneIDArray=downGeneIDArray[gopsByChip[chip]],positionNbArray=uint8(positionNbArray[gopsByChip[chip]]))
                    if log!=0:
                        log.write('%s\t%s\t%s\t20\tgops_by_region.bkdb\tr%s\t%u\t%.2f\n'%(date.today(),species,chip,region,len(gopsByChip[chip]),t2-t1))                                           
            
                    path=os.environ['PYDATA']+"/"+species+'/ensembl/'+chip+'_gops_by_region.bkdb'            
                    makebt(dict=dict[chip],path=path,newFlag=0,objFlag=1)                
    positionDB.close()    
    ensGeneDB.close()
    regionDB.close()
    if aceFlag:
        aceGeneDB.close()

if __name__=='__main__':
    import setenviron    
    setenviron.setenviron()   
    #positions_by_region('mouse',['m4', 'm5', 'm8', 'm27', 'm67', 'm70'],'62_37o','mus_musculus_funcgen_62_37o',['Mu11ksubA', 'MG-U74Av2', 'Mouse430_2', 'MOE430A', 'MoGene-1_0-st-v1', 'Mouse430A_2'],port=5306) 
    gop_by_region('mouse',1)    
    #merge_positions('mouse',['m4','m5','m8','m26','m27','m48','m49','m50','m53','m54','m62','m65','m67','m70'],'20',log=0)