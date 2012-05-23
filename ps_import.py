#!/usr/bin/python
"""**Create berkeley databases for chips, species, probes and probe sets.**"""

import os, sys,  shutil, cPickle, time
from datetime import date
import bsddb
try:
    import MySQLdb
except:
    pass
from dbtools import write_obj_db, makebt
import ps_class
import setenviron
setenviron.setenviron()
                                            
def make_chipdb(newFlag,log=0):
    """
    **Create chips database.**
    ::
    
        Arguments:
            newFlag  -- indicates if the data base is made from scratch (newFlag=1, existing
                        data base is erased) or is to be created or updated (newFlag=0)
            log ------- handle of a log file for recording messages
                
        Input:
            chip.txt
        
        Output:
            Create chip.bkdb
                key ==== chip name
                value == Chip object            
        
    """
           
    chipFile = open(os.environ['RAWDATA']+"/chip.txt")    
    lines=chipFile.readlines()    
    for i in range(len(lines)):        
        lines[i]=lines[i].split("\n")[0]        
        lines[i]=lines[i].split("\t")
    objClass='Chip'    
    objItems='myName=item[0],name=item[1],shortName=item[2],mySpeciesName=item[3],probesetNb=item[4],probeNb=item[5],compName=item[6],ens47Name=item[7],ens48Name=item[8],geoName=item[9]'    
    path=os.environ['COMMONDB']+"/chip.bkdb"    
    write_obj_db(lines, objClass, objItems,'string',path, newFlag,log=log)    
    
def make_speciesdb(newFlag,log=0):
    """
    **Create species database.**
    ::
    
        Arguments:
            newFlag  -- indicates if the data base is made from scratch (newFlag=1, existing
                        data base is erased) or is to be created or updated (newFlag=0)
            log ------- handle of a log file for recording messages
            
        Input:
            species.txt
        
        Output:
            species.bkdb
                key ==== species name
                value == Species object
    
    """
    
    speciesFile = open(os.environ['RAWDATA']+"/species.txt")    
    lines=speciesFile.readlines()    
    for i in range(len(lines)):        
        lines[i]=lines[i].split("\n")[0]        
        lines[i]=lines[i].split("\t")
    objClass='Species'    
    objItems='myName=item[0],officialName=item[1],ensName=item[2]'    
    path=os.environ['COMMONDB']+"/species.bkdb"    
    write_obj_db(lines, objClass, objItems,'string', path, newFlag,log=log)    
    
def make_probesetdb(chipList,log=0):
    """
    **Create probeset database.**
    ::
    
        Arguments:
            chipList  -- list the chips to be used for creating probe set databases
            log ------- handle of a log file for recording messages
                
        Input:
            %chip.bkdb
            %chip_probeset.txt
        
        Output:
            Create %chip_probeset.bkdb
                key ==== probe set id
                value == Probeset object
            Create chip_probeset.bkdb
                key ==== chip name
                value == list of probe set id
    
    """
    
        #open chip to recover species name
    path=os.environ['COMMONDB']+"/chip.bkdb"
    chipDB=bsddb.btopen(path,'r')
    #chip dictionnary                                
    chipDict={}                         
    for chipName in chipList:
        try:
            psFile = open(os.environ['RAWDATA']+"/%s_probeset.txt"%chipName)        
            currChip=cPickle.loads(chipDB[chipName])
            species=currChip.mySpeciesName
            probesetList=[]
            probesetDict={}
            index=-1                            
            for line in psFile.readlines():                    
                index+=1
                currProbeset=line.split("\n")[0]                               
                probesetDict[currProbeset]=ps_class.Probeset(probesetID=currProbeset,probesetIndex=index)                                                
                probesetList.append(currProbeset)
            path=os.environ['PYDATA']+"/" + species + "/" + chipName + "_probeset.bkdb"
            if log!=0:
                log.write('%s\t%s\t%s\t3\t%s_probeset.bkdb'%(date.today(),species,chipName,chipName))                                                                                
            makebt(dict=probesetDict,path=path,newFlag=1,objFlag=1,log=log)        
            print "*"*5                                                
            chipDict[chipName]=probesetList
        except:
            pass    
    path=os.environ['COMMONDB']+"/chip_probeset.bkdb"    
    if log!=0:
        log.write('%s\t\t\t3\tchip_probeset.bkdb'%date.today())
    makebt(dict=chipDict,path=path,newFlag=0,objFlag=1,log=log)    
                                                                                                                                                   
def make_probedb(ensTable,ensVersion,species,chipName,ensChipName,aceFlag,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):   
    """
    **Create probe database.**    
    ::
    
        Arguments:
            ensTable ---- Ensembl table to be interrogated by MySQL
            ensVersion -- Ensembl version
            species ---- species
            host ------- Ensembl database address    
            port ------- Ensembl port
            user ------- user name
            pswd ------- user pasword
            log -------- handle of a log file for recording messages
                
        Output:    
            Create %chip_probe.bkdb (AceView and Ensembl versions)
                key ==== probe id
                value == Probe object
                                
    """
     
    #print ensTable 
    #print host      
    try:
        connection = MySQLdb.connect(host = host,port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 4 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 4 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)
      
    # Create a mySQL cursor         
    cursor = connection.cursor()
    if int(ensVersion.split("_")[0])<48:
        query = """
        SELECT DISTINCT
        oligo_probe.name,
        oligo_probe.oligo_probe_id,    
        oligo_probe.probeset          
        FROM 
        oligo_probe, oligo_array 
        WHERE
        oligo_array.name='%s' 
        AND oligo_probe.oligo_array_id=oligo_array.oligo_array_id                
        ORDER BY oligo_probe.name
        """ %ensChipName        
    else:
        query = '''
        SELECT DISTINCT
        probe.name,
        probe.probe_id,    
        probe.probe_set_id          
        FROM 
        probe, array_chip 
        WHERE
        array_chip.name="%s" 
        AND probe.array_chip_id=array_chip.array_chip_id                    
        ORDER BY probe.name
        ''' %ensChipName
        
    t1=time.time()
    tryQuery=1
    done=0
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '4: lost MySQL connection'
            tryQuery=raw_input('would you like to continue (enter 1) or to cancel (enter 0)')
            if tryQuery=='0':
                tryQuery=0
                if log !=0:
                    log.write('\tCANCELED\t\n')  
                    sys.exit()
                sys.exit()                      
            elif tryQuery!="1":
                print 'enter 1 or 0'                    
    if done:
        t2=time.time()
        print "mySQL database %s delivered results in %.2f s\n"%("ensembl_"+ensTable,t2-t1)        
        # Get query's result
        #lines = cursor.fetchall()    
        # Initialize the object dictionary used to create the berkeley database
        dict={}
        # Put each resulting line in dictionary
        while (1):
            line=cursor.fetchone()
            #for n in lines:
            if line==None:
                break            
            probeID=line[0]            
            #some Affy chips do not have a probeID => probeID is made
            #with x and y position on chip and is terminated with a ;            
            if probeID[-1]==';':
                probeID=probeID[:-1]                 
            if ':' in probeID:          
                positions=probeID.split(":")  
                probeID="%04u:%04u"%(int(positions[0]),int(positions[1]))
            dict[probeID]=ps_class.Probe(probeID=probeID,ensemblID=line[1],ensProbesetID=line[2])                        
        cursor.close ()
        connection.close ()
        
        # create berkeley database
        path=os.environ['PYDATA'] + "/" + species +  "/ensembl/" + chipName+ "_probe.bkdb"            
        makebt(dict=dict,path=path,newFlag=1,objFlag=1,log=log)        
        if aceFlag:  
            #copy probeDB in AceView directoy     
            acePath=os.environ['PYDATA']+"/"+species+"/aceview/"+ chipName+ "_probe.bkdb"
            shutil.copyfile(path,acePath)

def fill_probesetdb(ensTable,ensVersion,species,ensSpecies,chipName,ensChipName,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Update probeset database.**    
    ::
    
        Arguments:
            ensTable ---- Ensembl table to be interrogated by MySQL
            ensVersion -- Ensembl version
            species ----- species name
            host -------- Ensembl database address    
            port -------- Ensembl port
            user -------- user name
            pswd -------- user password
            log --------- handle of a log file for recording messages
                
        Output:    
            Update %chip_probeset.bkdb
                key ==== probe set id
                value == Probeset object
                                    
    """ 
          
    try:    
        connection = MySQLdb.connect(host=host,port=port,user =user,db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 5 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 5 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)                 
    cursor = connection.cursor()
    if int(ensVersion)<50:        
        query = """
        SELECT DISTINCT 
        dbprimary_id,
        gene_stable_id 
        FROM 
        %s_gene_ensembl__xref_affy_%s__dm 
        where 
        dbprimary_id != 'null' 
        order by dbprimary_id;
        """ %(ensSpecies,ensChipName.lower())                            
    elif int(ensVersion) in (50,51,52,53,54,55):             
        query="""
        SELECT DISTINCT %s_gene_ensembl__ox_AFFY_%s__dm.display_label_1074,%s_gene_ensembl__gene__main.stable_id_1023
        FROM %s_gene_ensembl__ox_AFFY_%s__dm,%s_gene_ensembl__transcript__main,%s_gene_ensembl__gene__main
        WHERE %s_gene_ensembl__transcript__main.ox_AFFY_%s_bool=1
        AND %s_gene_ensembl__ox_AFFY_%s__dm.transcript_id_1064_key=%s_gene_ensembl__transcript__main.transcript_id_1064_key
        AND %s_gene_ensembl__transcript__main.gene_id_1020_key=%s_gene_ensembl__gene__main.gene_id_1020_key
        ORDER BY %s_gene_ensembl__ox_AFFY_%s__dm.display_label_1074,%s_gene_ensembl__gene__main.stable_id_1023    
        """%(ensSpecies,ensChipName,ensSpecies,
             ensSpecies,ensChipName,ensSpecies,ensSpecies,
             ensSpecies,ensChipName,
             ensSpecies,ensChipName,ensSpecies,
             ensSpecies,ensSpecies,
             ensSpecies,ensChipName,ensSpecies)            
    elif int(ensVersion)>=56:    
        if int(ensVersion) in (56,57,58,59):
            labelRank=11049
        else:
            labelRank=11056                              
        query="""
        SELECT DISTINCT %s_gene_ensembl__eFG_AFFY_%s__dm.display_label_%u,%s_gene_ensembl__gene__main.stable_id_1023
        FROM %s_gene_ensembl__eFG_AFFY_%s__dm,%s_gene_ensembl__transcript__main,%s_gene_ensembl__gene__main
        WHERE %s_gene_ensembl__transcript__main.eFG_AFFY_%s_bool=1
        AND %s_gene_ensembl__eFG_AFFY_%s__dm.transcript_id_1064_key=%s_gene_ensembl__transcript__main.transcript_id_1064_key
        AND %s_gene_ensembl__transcript__main.gene_id_1020_key=%s_gene_ensembl__gene__main.gene_id_1020_key
        ORDER BY %s_gene_ensembl__eFG_AFFY_%s__dm.display_label_%u,%s_gene_ensembl__gene__main.stable_id_1023    
        """%(ensSpecies,ensChipName,labelRank,ensSpecies,
             ensSpecies,ensChipName,ensSpecies,ensSpecies,
             ensSpecies,ensChipName,
             ensSpecies,ensChipName,ensSpecies,
             ensSpecies,ensSpecies,
             ensSpecies,ensChipName,labelRank,ensSpecies)                        
    else:
        print "mart version %u not supported"%ensVersion
        sys.exit(1)
    try:    
        os.mkdir(os.environ['PYDATA']+"/"+species+'/log')
    except:
        pass          
    log1=open(os.environ['PYDATA']+"/"+species+'/log/'+'%s_fill_probeset_log.txt'%chipName,'a')       
    t1=time.time() 
    tryQuery=1
    done=0
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '4: lost MySQL connection'
            tryQuery=raw_input('would you like to continue (enter 1) or to cancel (enter 0)')
            if tryQuery=='0':
                tryQuery=0
                if log !=0:
                    log.write('\tCANCELED\t\n') 
                    sys.exit()                   
            elif tryQuery!="1":
                print 'enter 1 or 0'                
    if done:           
        t2=time.time()
        print "mySQL database %s delivered results in %.2f s\n"%("ensembl_"+ensTable,t2-t1)    
        lines = cursor.fetchall()
        cursor.close ()
        connection.close ()
        
        # open the chip probeset data base to fill ensembl assignation     
        path=os.environ['PYDATA']+"/" + species+"/" + chipName+"_probeset.bkdb"            
        probesetDB=bsddb.btopen(path)
        psNb=len(probesetDB.keys())
        
        #empty ensemblGeneIDs attribute
        for probesetID in probesetDB.keys():
            probeset=cPickle.loads(probesetDB[probesetID])
            probeset.ensGeneIDs= None        
            probesetDB[probesetID]=cPickle.dumps(probeset,protocol=-1)
        lineNb=0
        addedPs=0
        for n in lines:
            lineNb+=1
            #some probe set registered in Ensembl do not exist on Affy
            try:        
                probeset=cPickle.loads(probesetDB[n[0]])
                if probeset.ensGeneIDs is None:
                    probeset.ensGeneIDs=[n[1]]
                else:
                    probeset.ensGeneIDs.append(n[1])
                probesetDB[n[0]]=cPickle.dumps(probeset,protocol=-1)
            except:
                addedPs+=1
                probeset=ps_class.Probeset(probesetID=n[0],ensGeneIDs=[n[1]],probesetIndex=psNb)
                probesetDB[n[0]]=cPickle.dumps(probeset,protocol=-1)
                psNb+=1
                log1.write('%s added\n'%n[0])                                        
        probesetDB.close()
        t3=time.time()
        if log !=0:
            log.write('%s\t%s\t%s\t5\t%s_probeset.bkdb\tmart_%s\t%u added\t%.2f\n'%(date.today(),species,chipName,chipName,ensVersion,addedPs,t3-t1))
                                        
#def probeset(probeFileName,species,chipName,probeLength,aceFlag,addProbeset,log=0):
def probeset(probeFileName,species,chipName,aceFlag,addProbeset,log=0):
    """
    **Update probeset database (fill probe information).**
    ::
    
        Arguments:
            probeFileName -- Affymetrix file describing probes
            species -------- species name
            chipName  ------ chip name
            probeLength ---- probe length
            log ------- handle of a log file for recording messages
                
        Input:
            %chip_probe.bkdb
            %chip_probeset.txt
        
        Output:
            Update %chip_probeset.bkdb
                key ==== probe set id
                value == Probeset object
            Update or Create  %chip_probe.bkdb (Ensembl and AceView version)
                key ==== probe set id
                value == Probe object
        
    """    
    
    # open the chip probe data base
    log1=open(os.environ['PYDATA']+"/"+species+'/log/'+'%s_probeset_log.txt'%chipName,'a')
    t1=time.time()
    probePath=os.environ['PYDATA']+"/"+species+"/ensembl/"+chipName+"_probe.bkdb"    
    probeDB=bsddb.btopen(probePath)    
    # open the chip probeset data base 
    probesetPath=os.environ['PYDATA']+"/"+species+"/"+chipName+"_probeset.bkdb"        
    probesetDB=bsddb.btopen(probesetPath) 
    psNb=len(probesetDB.keys())
    addedPs=0 
    addedP=0  
    # open the file describing the probes (e.g. MG_U74Av2_probe_tab, an Affymetrix tab file)
    # for Affymetrix the probes are ordered inside the file according to their target position
    # which corresponds to the probe rank (here index= rank -1)
    chipFile = open(probeFileName)            
    #initialize not to have error in pydev        
    probesetID=""    
    probeIDs=[]
    ensemblProbeIDs=[]
    probeIndexes=[]
    targetPositions=[]
    startPositions=[]            
    probesetNb=0 
    #analyse first line    
    fields = chipFile.readline().split('\t') 
    pPos=-1
    psPos=-1
    xPos=-1
    yPos=-1
    seqPos=-1
    targPos=-1
    startPos=-1   
    #in gene chips (e.g. MoGene-1_0-st-v1), probe names are composit = probe ID + transcript cluster ID
    #(probeID in Affymetrix probe_tab is not unique)  
    compositFlag=0
    for i in range(len(fields)):
        if fields[i].lower()=='probe set name' or fields[i].lower()=='probe set id' or fields[i].lower()=='transcript cluster id':
            psPos=i
            if fields[i].lower()=='transcript cluster id':
                compositFlag=1            
        if fields[i].lower()=='probe id':
            pPos=i            
        if fields[i].lower()=='probe x':
            xPos=i            
        if fields[i].lower()=='probe y':
            yPos=i            
        if fields[i].lower()=='probe sequence':
            seqPos=i
        if fields[i].lower()=='probe interrogation position':
            targPos=i
        if fields[i].lower()=='start':
            startPos=i
    exitSys=0        
    if seqPos==-1:        
        log1.write('NO PROBE SEQUENCE\n')        
        print 'NO PROBE SEQUENCE'
        print repr(fields)
        exitSys=1            
    if psPos==-1:        
        log1.write('NO PROBESETID\n')        
        print 'NO PROBESET ID'
        print repr(fields)
        exitSys=1
    if exitSys:        
        sys.exit()                                                                                                                         
    while 1:        
        line = chipFile.readline()        
        if not line:
            break            
        items=line.split("\t")        
        currProbesetID=items[psPos]                                        
        if currProbesetID!=probesetID:
            #save previous probeset 
            #Add it if it is not registered and if it is allowed by addProbeset
            if probesetNb>0:                                            
                #update probe set database
                updateFlag=1  
                try:                     
                    probeset=cPickle.loads(probesetDB[probesetID])
                except:
                    if addProbeset:
                        addedPs+=1
                        probeset=ps_class.Probeset(probesetID=probesetID,probesetIndex=psNb)
                        psNb+=1                         
                        log1.write('%s added\n'%probesetID)
                    else:
                        updateFlag=0
                        #correct probesetNb which has been incremented by default
                        probesetNb-=1
                if updateFlag:                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    probeset.probeIDs=probeIDs
                    probeset.ensProbeIDs=ensemblProbeIDs
                    probeset.probeIndexes=probeIndexes
                    if probeIDs!=[]:
                        if targetPositions!=[]:
                            probeset=set_probeset_info(probeset,probeIDs,targetPositions,probeLength,'transcript')
                        elif startPositions!=[]: 
                            probeset=set_probeset_info(probeset,probeIDs,startPositions,probeLength,'genome')
                        else:
                            probeset=set_probeset_info(probeset,probeIDs,0,0,'none')
                    probesetDB[probesetID]=cPickle.dumps(probeset,protocol=-1)
            #initialize the current probe set
            probesetNb+=1    
            probesetID=currProbesetID            
            probeIDs=[]
            ensemblProbeIDs=[]
            probeIndexes=[]
            targetPositions=[]            
            index = -1            
        index+=1
        probeIndexes.append(index)        
        if pPos>=0:
            if compositFlag:
                probeID=items[pPos]+'-'+items[psPos]
            else:
                probeID=items[pPos]
        elif xPos>=0:
            probeID='%04u:%04u'%(int(items[xPos]),int(items[yPos]))
        else:            
            log1.write('NO PROBE ID\n')            
            print 'NO PROBE ID'
            sys.exit()
        probeIDs.append(probeID)
        #update probe database  
        updateFlag=1          
        try:
            probe=cPickle.loads(probeDB[probeID])
            #verify the probesetID
            if probe.probesetID!=None:
                if probe.probesetID!=probesetID:
                    print 'probe %s (%04u:%04u) belongs to %s in probe_tab and to %s in Ensembl'%(probeID,probe.yPosition,probe.xPosition,probesetID,probe.probesetID)                     
                    log1.write('probe %s (%04u:%04u) belongs to %s in probe_tab and to %s in Ensembl\n'%(probeID,probe.yPosition,probe.xPosition,probesetID,probe.probesetID))
            else:
                probe.probesetID=probesetID                                                                                        
        except:
            if addProbeset:
                #probe was not registered in Ensembl => create it 
                addedP+=1                                    
                log1.write('%s added\t%s\n'%(probeID,probesetID)) 
                probe=ps_class.Probe(probesetID=probesetID,probeID=probeID)
                probe.probesetID=probesetID
            else:
                updateFlag=0
        if updateFlag:                
            ensemblProbeIDs.append(probe.ensemblID)
            if targPos>=0:                                                           
                probe.targetPosition=int(items[targPos])
                targetPositions.append(probe.targetPosition)
            elif startPos>=0:
                try:
                    probe.targetPosition=int(items[startPos])
                    startPositions.append(probe.targetPosition)
                except:
                    #unknown position (--- in probe tab)
                    pass
                 
            probe.sequence=items[seqPos]
    #        if probeLength!=len(probe.sequence):
    #            print "WARNING probe sequence length = %u ~= probeLength = %u"%(len(probe.sequence),probeLength)
            probeLength=len(probe.sequence)
            probe.index=index
            if xPos>=0:
                probe.xPosition=items[xPos]
            if yPos>=0:
                probe.xPosition=items[yPos]            
            probeDB[probeID]=cPickle.dumps(probe,protocol=-1)
                                                                              
    #save last probeset
    if probesetNb>0:                
        #update probe set database
        updateFlag=1     
        try:                            
            probeset=cPickle.loads(probesetDB[probesetID])
        except:     
            if addProbeset:           
                addedPs+=1
                probeset=ps_class.Probeset(probesetID=probesetID,probesetIndex=psNb)
                psNb+=1  
                print 'probeset ',probesetID,'added'   
                log1.write('%s added\n'%probesetID)
            else:
                updateFlag=0
        if updateFlag:                                                    
            probeset.probeIDs=probeIDs
            probeset.probeIndexes=probeIndexes
            if probeIDs!=[]:            
                if targetPositions!=[]:
                    probeset=set_probeset_info(probeset,probeIDs,targetPositions,probeLength,'transcript')
                elif startPositions!=[]: 
                    probeset=set_probeset_info(probeset,probeIDs,startPositions,probeLength,'genome')
                else:
                    probeset=set_probeset_info(probeset,probeIDs,0,0,'none')            
            probesetDB[probesetID]=cPickle.dumps(probeset,protocol=-1)             
    probeDB.close()
    probesetDB.close()  
    if log !=0:
        t2=time.time()
        log.write('%s\t%s\t%s\t6\t%s_probeset.bkdb\t\tadded %u\t%.2f\n'%(date.today(),species,chipName,chipName,addedPs,t2-t1))
        if addedP>0:
            log.write('%s\t%s\t%s\t6\t%s_probe.bkdb\t\tadded %u\t%.2f\n'%(date.today(),species,chipName,chipName,addedP,t2-t1))
    if aceFlag:          
        #copy probeDB in AceView     
        acePath=os.environ['PYDATA']+"/"+species+"/aceview/"+chipName+"_probe.bkdb"
        shutil.copyfile(probePath,acePath)

def set_probeset_info(probeset,probeIDs,position,probeLength,type):
    """
    **Update probe set information (calculate target genomic start and end).**
    ::
    
        Arguments:
            probeset ----- Probeset object
            probeIDs ----- list of probe ids
            position ----- list of probe position in transcript or in genome
            probeLength -- probe length
            type --------- 'transcript' or 'genome'
            
    """    
    
    # probe nb
    probeset.probeNb=len(probeIDs)
    if type!='none':                              
        # position = probe position in transcript or in genome
        # probeLength = probe length      
        minPosition=position[0]
        maxPosition=minPosition    
        for i in range(1,len(position)):
            minPosition=min(minPosition,position[i])
            maxPosition=max(maxPosition,position[i])
        probeset.probesetTargetLength=maxPosition-minPosition+probeLength    
        if type=='transcript':                            
            probeset.probesetTargetStart=int(minPosition-round(probeLength/2))
            probeset.probesetTargetEnd=int(maxPosition+round(probeLength/2))
        elif type=='genome':                    
            probeset.probesetTargetStart=minPosition
            probeset.probesetTargetEnd=maxPosition+probeLength
    return probeset        


if __name__=='__main__':
    make_chipdb(1)
    #fill_probesetdb("ensembl_mart_62",62,'mouse','mmusculus','m65','MoEx_1_0_st_v1',host='martdb.ensembl.org',port=5316)
    #fill_probesetdb("ensembl_mart_62",62,'mouse','mmusculus','m67','MoGene_1_0_st_v1',host='martdb.ensembl.org',port=5316)
    #make_chipdb(1)
    #make_probesetdb(['m65','m67'])
    #make_probedb('homo_sapiens_funcgen_63_37','63_37','human','m2','HG_U95A')
