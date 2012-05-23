#!/usr/bin/python
"""
**Control data import.**
::

    Parameters:
    actions       -- a list of action designed by numbers
            1 : make chip.bkdb
            2 : make species.bkdb
            3 : make chip_probe.bkdb
    mySpeciesName -- species name used inside program
    myChipName    -- chip name used inside program
    ensVersion    -- Ensembl version to be interrogated
    
    Actions:
    
    1 -- 
    Example:
    pscontrol_chip [3] mouse m8 62_37o
"""
import os, time, sys, shutil, bsddb, cPickle
try:
    import MySQLdb
except:
    pass
from datetime import date
import getopt
import setenviron
setenviron.setenviron()
import ps_import
import ps_ensembl
import ps_aceview
import ps_position
import ps_probeset
import ps_export

if __name__=='__main__':
    #control argument list
    if len(sys.argv)==1:
        sys.exit("pscontrol_chip needs at least actions parameters")
        
    myChipName=''
    mySpeciesName=''
    currChip=''     
    myChipList=''
    ensVersion=''
    region=''
    log=''
    subAction='a'
    Ok=0
    aceFlag=0
    addProbeset=0;
    
    usage='psawn uses -a -b -c -e -f -l -p -r -s -t -v -h or --action --sub --chip --ensembl --aceflag --chiplist --addprobeset --region --species --tscript --aceview --help'
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:c:e:f:l:p:r:s:t:v:h:", ["action=","sub=","chip=","ensembl=","aceflag","chiplist=","addprobeset=","region=","species=","tscript=","aceview=","help="])
    except getopt.getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print usage
        sys.exit(2)    
    for opt, arg in opts:        
        if opt in ("-h", "--help"):
            print usage            
            sys.exit()
        elif opt in ("-a", "--action"):
            action= int(arg)
        elif opt in ("-b", "--sub"):
            subAction= arg    
        elif opt in ("-c", "--chip"):
            myChipName=arg
            """
            Recover Ensembl chip name which is in general equal to  
            chip.shortName, the bracketed name of the corresponding
            GPL platform in GEO.        
            """  
            if action>0 or (action==0 and subAction>'a'):              
                path=os.environ['COMMONDB']+"/chip.bkdb"           
                chip=bsddb.btopen(path,'r')                
                currChip=cPickle.loads(chip[myChipName])
                mySpeciesName=currChip.mySpeciesName
                shortChipName=currChip.shortName
                probeNb=currChip.probeNb
                #probeLength=currChip.probeLength
        elif opt in ("-f", "--aceflag"):
            aceFlag= int(arg)
        elif opt in ("-p", "--addprobeset"):
            addProbeset= int(arg)            
        elif opt in ("-r", "--region"):
            region=arg                               
            if action==23:     
                psRange=eval(arg)                                            
            if action==25:     
                psRange=eval(arg)    
        elif opt in ("-s", "--species"):
            mySpeciesName=arg            
        elif opt in ("-t", "--tscript"):
            tscriptFile=arg                
        elif opt in ("-e", "--ensembl"):                        
            ensVersion=arg                        
            if int(ensVersion.split("_")[0])<48:
                if currChip!='':
                    ensChipName=currChip.ens47Name
                ensPort=3306
                martPort=3316                                                                    
            else:
                if currChip!='':
                    ensChipName=currChip.ens48Name
                ensPort=5306
                martPort=5316
        elif opt in ("-l", "--chiplist"):
            myChipList=eval(arg)           
            if action>0 or (action==0 and subAction>'a'):
                path=os.environ['COMMONDB']+"/chip.bkdb"
                chip=bsddb.btopen(path,'r')            
                ensChipList=list()
                for myChipName in myChipList:
                    currChip=cPickle.loads(chip[myChipName])
                    if ensVersion!='':                                
                        if int(ensVersion.split("_")[0])<48:                            
                            ensChipList.append(currChip.ens47Name)
                        else:                            
                            ensChipList.append(currChip.ens48Name)
                currChip=''    
        elif opt in ("-v", "--aceview"):
            aceVersion=arg                                     
        else:
            assert False, "unhandled option"    
    if mySpeciesName!='':        
        """
        Recover the species name used in Ensembl, which is
        the lower case official species name with underscores
        """                            
        path=os.environ['COMMONDB']+"/species.bkdb"
        species=bsddb.btopen(path,'r')                
        currSpecies=cPickle.loads(species[mySpeciesName]) 
        officialSpeciesName=currSpecies.officialName
        ensSpeciesName=officialSpeciesName.lower()
        ensSpeciesName=("_").join(ensSpeciesName.split(" "))    
                               
    if action==4 or action==5:        
        if myChipName=='':
            print("myChipName must be specified with -c or --chip if action 4 or 5 are selected")
            sys.exit()          
        if mySpeciesName=='':
            print("mySpeciesName must be specified with -s or --species if action 4 or 5 are selected")
            sys.exit()              
        if ensVersion=='':
            print("ensVersion must be specified with -e or --ensembl if action 4 or 5 are selected")
            sys.exit()
    
            
    os.chdir(os.environ['PYDATA'])
    if myChipName!='':
        log=open('log_%s.txt'%myChipName,'a')
    else:
        log=open('log_%s.txt'%mySpeciesName,'a')
    
    if action==0:
        
        if subAction=='a':
            print 'PROCESS USER DATA'
        elif subAction=='b':
            print 'PROCESS EXTERNAL DATA'
        elif subAction=='c':
            print 'PROCESS ENSEMBL DATA'
        elif subAction=='d':
            print 'PROCESS ACEVIEW DATA'
        elif subAction=='e':
            print 'PROCESS POSITION DATA'
        elif subAction=='f':
            print 'PROCESS ENSEMBL PROBESET DATA'            
        elif subAction=='g':
            print 'PROCESS ACEVIEW PROBESET DATA'            
        elif subAction=='h':
            print 'OUTPUT DATA'            
   
   
    """CREATE COMMON DATA BASES FROM USER DATA"""
                                                
    if action==1 or (action==0 and subAction=='a'):
        
        print "*"*5+' 1: RECOVER CHIP INFORMATION'
        log.write('%s\t\t\t1\tchip.bkdb\t'%date.today())        
        
        Ok=ps_import.make_chipdb(1,log)
                
        if Ok!=None:
            try:
                os.remove(os.environ['COMMONDB']+"/species.bkdb")
            except:
                pass    
            log.write('ACTION 1 FAILED\n')
            print('ACTION 1 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
                    
        print '%s\n'%('*'*5)                   
        
    if action==2 or (action==0 and subAction=='a'):
          
        print "*"*5+' 2: RECOVER SPECIES INFORMATION'        
        log.write('%s\t\t\t2\tspecies.bkdb\t'%date.today())        
        
        Ok=ps_import.make_speciesdb(1,log)
                
        if Ok!=None:
            try:
                os.readlink(os.environ['COMMONDB']+"/species.bkdb")
            except:
                pass
            log.write('ACTION 2 FAILED\n')
            print('ACTION 2 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
            
        print '%s\n'%('*'*5)        
        
    if action==3 or (action==0 and subAction=='a'):
          
        print "*"*5+' 3: RECOVER PROBESETS OF THE CHIP LIST'
                                                                    
        Ok=ps_import.make_probesetdb(myChipList,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['COMMONDB']+"/chip_probeset.bkdb")
            except:
                pass
            log.write('ACTION 3 FAILED\n')
            print('ACTION 3 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
                                  
        for myChipName in myChipList:
            try:
                oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
                newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step3.bkdb"
                os.rename(oldPath, newPath)
            except:
                log.write('CHIPSET %s CREATION FAILED\n'%myChipName)
                print('CHIPSET %s CREATION FAILED\n'%myChipName)
                    
        print '%s\n'%('*'*5)                    
        
   
    """CREATE COMMON DATA BASES FROM ENSEMBL DATA"""
    
    if action==4 or (action==0 and subAction=='b'):
        
        ensTable=ensSpeciesName+"_funcgen_"+str(ensVersion)                                      
        print "*"*5+' 4 RECOVER PROBES FOR ' + myChipName + ' FROM ' + ensTable                                          
        log.write('%s\t%s\t%s\t4\t%s_probe.bkdb\tens_%s'%(date.today(),mySpeciesName,myChipName,myChipName,ensVersion))                                      
        Ok=ps_import.make_probedb(ensTable,ensVersion,mySpeciesName,myChipName,ensChipName,aceFlag,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb")
            except:
                pass
            if aceFlag:
                try:
                    os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb")
                except:
                    pass
            log.write('ACTION 4 FAILED\n')
            print('ACTION 4 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
        
        oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb"
        newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe_step4.bkdb"
        os.rename(oldPath, newPath)
        if aceFlag:
            oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb"
            newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe_step4.bkdb"
            os.rename(oldPath, newPath)
                
        print '%s\n'%('*'*5)           
    
    if action==5 or (action==0 and subAction=='b'):
                
        print "*"*5+' 5: FILL PROBESET WITH ENSEMBL ASSIGNATION  FOR '+myChipName
                        
        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step3.bkdb"
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            log.write('ACTION 5 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 5 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()                 
                         
        ensVersion=ensVersion.split("_")[0]
        ensTable="ensembl_mart_%s"%ensVersion
        ensSpeciesName=ensSpeciesName.split("_")[0][0]+ensSpeciesName.split("_")[1]
        #in Mart DB only underscore is used  
        ensChipName=ensChipName.replace('-','_')               
                                        
        Ok=ps_import.fill_probesetdb(ensTable,ensVersion,mySpeciesName,ensSpeciesName,myChipName,ensChipName,host='martdb.ensembl.org',port=martPort,log=log)
        
        if Ok!=None:
            os.remove(os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb")                                    
            log.write('ACTION 5 FAILED\n')
            print('ACTION 5 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
        
        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step5.bkdb"
        os.rename(oldPath, newPath)                            
            
        print '%s\n'%('*'*5)
                   
    if action==6 or (action==0 and subAction=='b'):
        
        probeFile = os.environ['AFFYDATA']+"/%s.probe_tab" %shortChipName
        print "*"*5+' 6: FILL '+myChipName +' PROBESET  WITH AFFYMETRIX DATA IN '+probeFile
                        
        oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe_step4.bkdb"
        newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            log.write('ACTION 4 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 4 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()          
        if aceFlag:              
            oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe_step4.bkdb"
            newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb"
            try:
                shutil.copyfile(oldPath,newPath)            
            except:
                os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb")
                log.write('ACTION 4 FAILED (%s does not exist)\n'%oldPath)
                print('ACTION 4 FAILED (%s does not exist)\n'%oldPath)
                print '%s\n'%('*'*5)
                sys.exit()
                            
        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step5.bkdb"
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb")
            os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb")
            log.write('ACTION 5 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 5 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()                                                    
            
        #Ok=ps_import.probeset(probeFile,mySpeciesName,myChipName,probeLength,aceFlag,addProbeset,log=log)
        Ok=ps_import.probeset(probeFile,mySpeciesName,myChipName,aceFlag,addProbeset,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb")
            except:
                pass
            try:                
                os.remove(os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb")
            except:
                pass
            if aceFlag:
                try:                        
                    os.remove(path=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb")
                except:
                    pass
            log.write('ACTION 6 FAILED\n')
            print('ACTION 6 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()

        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step6.bkdb"
        os.rename(oldPath, newPath)                
        oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe.bkdb"
        newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/ensembl/" + myChipName+ "_probe_step6.bkdb"
        os.rename(oldPath, newPath)
        if aceFlag:
            oldPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe.bkdb"
            newPath=os.environ['PYDATA'] + "/" + mySpeciesName +  "/aceview/" + myChipName+ "_probe_step6.bkdb"
            os.rename(oldPath, newPath)            
    
        print '%s\n'%('*'*5)
        
    
    """IMPORT ENSEMBL INFORMATION ON GENES AND TRANSCRIPTS"""
        
    if action==7  or (action==0 and subAction=='c'):
        
        print "*"*5+' 7: ENSEMBL TRANSCRIPTS BY EXON for '+mySpeciesName                                        
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t7\t%s_transcripts_by_exon.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))                                 
             
        Ok=ps_ensembl.transcripts_by_exon(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_transcripts_by_exon.bkdb")
            except:
                pass
            log.write('ACTION 7 FAILED\n')
            print('ACTION 7 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
            
        print '%s\n'%('*'*5)
        
    if action==8 or (action==0 and subAction=='c') :
        
        print "*"*5+' 8: ENSEMBL EXONS BY GENE for '+mySpeciesName                        
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t8\t%s_exons_by_gene.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))                        
             
        Ok=ps_ensembl.exons_by_gene(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb")
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb")
            except:
                pass                
            log.write('ACTION 8 FAILED\n')
            print('ACTION 8 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()  
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb"                                      
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene_step8.bkdb"
        os.rename(oldPath, newPath)
    
        print '%s\n'%('*'*5)             
        
    if action==9 or (action==0 and subAction=='c') :
        
        print "*"*5+' 9: ENSEMBL GENES BY REGION for '+mySpeciesName                
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t9\t%s_genes_by_region.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))                 
        
        Ok=ps_ensembl.genes_by_region(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_genes_by_region.bkdb")
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_region_by_gene.bkdb")
            except:
                pass
            log.write('ACTION 9 FAILED\n')
            print('ACTION 9 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
        
        print '%s\n'%('*'*5)                 
        
    if action==10 or (action==0 and subAction=='c') :
        
        print '*'*5+' 10: EXONS BY PREDICTED TRANSCRIPT for '+mySpeciesName            
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t10\t%s_exons_by_predicted_transcript.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))
        
        Ok=ps_ensembl.exons_by_predicted_transcript(ensTable,mySpeciesName,port=ensPort,log=log)

        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_predicted_transcript.bkdb")
            except:
                pass
            log.write('ACTION 10 FAILED\n')
            print('ACTION 10 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
        
        print '%s\n'%('*'*5)             
           
    if action==11 or (action==0 and subAction=='c'):
        
        print '*'*5+' 11: PREDICTED TRANSCRIPTS BY REGION for '+mySpeciesName        
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t11\t%s_predicted_transcript_by_region.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))
                            
        Ok=ps_ensembl.predicted_transcripts_by_region(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_predicted_transcripts_by_region.bkdb")
            except:
                pass
            log.write('ACTION 11 FAILED\n')
            print('ACTION 11 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
        
        print '%s\n'%('*'*5)             
         
    if action==12 or (action==0 and subAction=='c'):
        
        print '*'*5+' 12: TRANSCRIPTS BY GENE for '+mySpeciesName  
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t12\t%s_transcripts_by_gene.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))
                  
        Ok=ps_ensembl.transcripts_by_gene(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_transcripts_by_gene.bkdb")
            except:
                pass
            log.write('ACTION 12 FAILED\n')
            print('ACTION 12 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
        
        print '%s\n'%('*'*5)             
                 
    if action==13 or (action==0 and subAction=='c'):
        
        print '*'*5+' 13: CHROMOSOME BY REGION for '+mySpeciesName
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)
        log.write('%s\t%s\t\t13\t%s_chromosome_by_region.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))
                
        Ok=ps_ensembl.chromosome_by_region(ensTable,mySpeciesName,port=ensPort,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_chromosome_by_region.bkdb")
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_region_by_chromosome.bkdb")
            except:
                pass                    
            log.write('ACTION 13 FAILED\n')
            print('ACTION 13 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     
        
        print '%s\n'%('*'*5)             
                              
    if action==14 or (action==0 and subAction=='c'):
        
        print '*'*5+' 14: TRANSCRIPTS SEQUENCES for '+mySpeciesName
        ensTable="%s_core_%s"%(ensSpeciesName,ensVersion)        
        try:
            connection = MySQLdb.connect(host='ensembldb.ensembl.org', port=ensPort,  user='anonymous', db = ensTable)
        except MySQLdb.Error, e: 
            log.write('ACTION 14 ENSEMBL CONNECTION FAILED\n')            
            print "ACTION 14 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])
            sys.exit()                 
        cursor = connection.cursor()
        
        query= '''
        SELECT DISTINCT
        version
        FROM
        coord_system
        WHERE attrib="default_version"
        AND name="chromosome";        
        '''        
        cursor.execute (query)
        lines = cursor.fetchall()
        ncbi=lines[0][0]    
        cursor.close ()
        connection.close ()        
        tscriptFile=os.environ['ENSDATA']+'/'+officialSpeciesName.replace(' ','_')+'.'+ncbi+'.'+ensVersion.split('_')[0]+'.cdna.all.fa'        
        log.write('%s\t%s\t\t14\t%s_transcript_sequence.bkdb\tens_%s'%(date.today(),mySpeciesName,mySpeciesName,ensVersion))        
                
        Ok=ps_ensembl.transcript(tscriptFile,ncbi,mySpeciesName,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_transcript_sequence.bkdb")
            except:
                pass
            log.write('ACTION 14 FAILED\n')
            print('ACTION 14 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                     

        print '%s\n'%('*'*5)             
   
        
    """IMPORT ACEVIEW INFORMATION ON GENES AND TRANSCRIPTS"""
             
    if action==15 or (action==0 and subAction=='d'):
        
        print '*'*5+' 15: ACEVIEW TRANSCRIPT SEQUENCES for '+mySpeciesName
        log.write('%s\t%s\t\t15\t%s_transcript_sequence.bkdb\tace_%s'%(date.today(),mySpeciesName,mySpeciesName,aceVersion)) 
                                                    
        Ok=ps_aceview.transcript_sequence(mySpeciesName,aceVersion,log=log)

        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_transcript_sequence.bkdb")
            except:
                pass
            log.write('ACTION 15 FAILED\n')
            print('ACTION 15 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit() 

        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_transcript_sequence.bkdb"
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_transcript_sequence_step15.bkdb"
        os.rename(oldPath, newPath)                                        
        
        print '%s\n'%('*'*5)                    
         
    if action==16 or (action==0 and subAction=='d'):
         
        print '*'*5+' 16: ACEVIEW GENES BY ENSEMBL REGION for '+mySpeciesName
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_transcript_sequence_step15.bkdb"
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_transcript_sequence.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            log.write('ACTION 16 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 16 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()      
                                         
        Ok=ps_aceview.genes_exons(mySpeciesName,aceVersion,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+'/'+mySpeciesName+'/aceview/'+mySpeciesName+'_genes_by_ensembl_region.bkdb')
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+'/'+mySpeciesName+'/aceview/'+mySpeciesName+'_exons_by_gene.bkdb')
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+'/aceview/'+mySpeciesName+'_transcripts_by_gene.bkdb')
            except:
                pass
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+'/aceview/'+mySpeciesName+'_transcripts_by_exon.bkdb')
            except:
                pass            
            log.write('ACTION 16 FAILED\n')
            print('ACTION 16 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()   

        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene.bkdb"                                      
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene_step16.bkdb"
        os.rename(oldPath, newPath)
                                                
        print '%s\n'%('*'*5)
         
    if action==17 or (action==0 and subAction=='d'):
          
        print '*'*5+' 17: ENSEMBL GENES BY ACEVIEW GENE for '+mySpeciesName
                         
        Ok=ps_aceview.ensembl_genes(mySpeciesName,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+'/'+mySpeciesName+'/aceview/'+mySpeciesName+'_ensembl_genes_by_gene.bkdb')
            except:
                pass            

            log.write('ACTION 17 FAILED\n')
            print('ACTION 17 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
                                 
        print '%s\n'%('*'*5)
        
    """PROCESS POSITIONS"""
                            
    if action==18 or (action==0 and subAction=='e'):
        
        print '*'*5+' 18: PROBE POSITIONS BY REGION for chip list for '+ mySpeciesName 
        ensTable=ensSpeciesName+"_funcgen_"+str(ensVersion)       
        
        Ok=ps_position.positions_by_region(mySpeciesName,myChipList,ensVersion,ensTable,ensChipList,port=ensPort,log=log)
                                            
        if Ok!=None:
            log.write('ACTION 18 FAILED\n')
            print('ACTION 18 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()   
                              
        for myChipName in myChipList:
            if not os.path.exists(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_positions_by_region.bkdb"):
                log.write('POSITION %s CREATION FAILED\n'%myChipName)
                print('POSITION %s CREATION FAILED\n'%myChipName)                              
                                            
        print '%s\n'%('*'*5)
              
    if action==19 or (action==0 and subAction=='e'):
        
        print '\n*'*5+' 19: MERGE PROBE POSITIONS BY REGION for chip list for '+ mySpeciesName        
        path=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_chromosome_by_region.bkdb"                    
        regionDB=bsddb.btopen(path)
        regions=regionDB.keys()
        regionDB.close()     
                        
        if subAction=='a' or action==0:            
            t1=time.time()
                 
            for region in regions:                                                                                               
                Ok=ps_position.merge_positions(mySpeciesName,myChipList,region,log=log)
                if Ok!=None:
                    try:
                        os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl"+'%s_positions_by_region.bkdb'%mySpeciesName)
                    except:
                        pass
                    log.write('ACTION 19 FAILED FOR REGION %s\n'%region)
                    print('ACTION 19 FAILED FOR REGION %s\n'%region)                    
                    print '%s\n'%('*'*5)
                    sys.exit()          
                                               
            print '%s\n'%('*'*5)
            t2=time.time()
            log.write('%s\t%s\t\t19\t%s_positions_by_region.bkdb\t\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))
                                
        elif subAction=='b':
            #prepare files for MPI
            #bat file
            batFid=open('%s/%s_mp.bat'%(os.environ['SCRATCH'],mySpeciesName),'w')
            for region in regions:
                batFid.write('python %s/psawn.py -a 19 -b c -l "%s" -r %s -s %s\n'%(os.environ['PROG'],myChipList,region,mySpeciesName))
            #pbs file
            jobFid=open('%s/%s_mp.pbs'%(os.environ['SCRATCH'],mySpeciesName),'w')
            jobFid.write('#PBS -S /bin/bash\n')
            jobFid.write('#PBS -N mp\n')
            jobFid.write('#PBS -e mp.err\n')
            jobFid.write('#PBS -o mp.log\n')
            jobFid.write('#PBS -l walltime=02:00:00\n')
            jobFid.write('#PBS -l select=3:ncpus=8:mpiprocs=8\n')
            jobFid.write('module unload intel\n')
            jobFid.write('module load intel/11.1.072\n')
            jobFid.write('module load numpy/1.5.1\n')
            jobFid.write('module load scipy/0.8.0\n')
            jobFid.write('module load pserie/9.3.26\n')                        
            jobFid.write('cd %s\n'%os.environ['SCRATCH'])        
            jobFid.write('mpiexec pserie %s_mp.bat'%(mySpeciesName))
            jobFid.close()
                        
        elif subAction=='c':    
            #use MPI files
            ps_position.merge_positions(mySpeciesName,myChipList,region)                        
        
    if action==20 or (action==0 and subAction=='e'):
                                    
        print '*'*5+' 20 ENSEMBL POSITION GOP MAPPING FOR CHIP LIST'                                     
        t1=time.time()
        
        Ok=ps_position.gop_by_region(mySpeciesName,aceFlag,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+'/ensembl/'+chip+'_gops_by_region.bkdb')
            except:
                pass
            log.write('ACTION 20 FAILED\n')
            print('ACTION 20 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                             
                                              
        t2=time.time()
        log.write('%s\t%s\t\t20\tgops_by_region.bkdb\t\t\t%.2f\n'%(date.today(),mySpeciesName,t2-t1))
        
        print '%s\n'%('*'*5)                                                 
   
                                                     
    """PROBESET ASSIGNATION ON ENSEMBL GENES"""
    
    if action==21 or (action==0 and subAction=='f'):
        
        
        print '*'*5+' 21: ENSEMBL POSITION MAPPING for '+myChipName
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene_step8.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            log.write('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()
            
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe_step6.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb")
            log.write('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()     
            
        path=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_positions_by_region.bkdb"            
        regionDB=bsddb.btopen(path)
        regions=regionDB.keys()
        regionDB.close()                                                                      
        t1=time.time()     
        
        for region in regions:
            Ok=ps_probeset.position_mapping(mySpeciesName,myChipName,'ensembl',region,log=log)
            if Ok!=None:
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+mySpeciesName+"_exons_by_gene.bkdb")
                except:
                    pass
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb")
                except:
                    pass
                log.write('ACTION 21 FAILED FOR REGION %s\n'%region)
                print('ACTION 21 FAILED FOR REGION %s\n'%region)
                print '%s\n'%('*'*5)
                sys.exit()
                                             
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb"
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe_step21.bkdb"
        os.rename(oldPath, newPath)                 
                          
        t2=time.time()
        log.write('%s\t%s\t\t21\t%s_exons_by_gene.bkdb\t\tensembl update\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))
        log.write('%s\t%s\t%s\t21\t%s_probe.bkdb\t\tensembl update\t%.2f\n'%(date.today(),mySpeciesName,myChipName,myChipName,t2-t1))
        
        print '%s\n'%('*'*5)
            
        
    """PROBESET ASSIGNATION ON GOP GENES"""
            
    if action==22 or (action==0 and subAction=='f'):
        
        print '*'*5+' 22: GOP POSITION MAPPING for '+myChipName
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe_step21.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:            
            log.write('ACTION 22 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 22 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()
                                                  
        path=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_positions_by_region.bkdb"                
        regionDB=bsddb.btopen(path)
        regions=regionDB.keys()                                                                          
        t1=time.time()     
        
        for region in regions:
            Ok=ps_probeset.position_gopmapping(mySpeciesName,myChipName,region)
            if Ok!=None:
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb")
                except:
                    pass
                log.write('ACTION 22 FAILED FOR REGION %s\n'%region)
                print('ACTION 22 FAILED FOR REGION %s\n'%region)
                print '%s\n'%('*'*5)
                sys.exit()
                
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb"
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe_step22.bkdb"
        os.rename(oldPath, newPath)
                         
        t2=time.time()
        log.write('%s\t%s\t\t22\t%s_exons_by_gene.bkdb\t\tensembl update\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))
        log.write('%s\t%s\t%s\t22\t%s_probe.bkdb\t\tensembl update\t%.2f\n'%(date.today(),mySpeciesName,myChipName,myChipName,t2-t1))                                      
                                    
        print '%s\n'%('*'*5)  
        
                        
    if action==23 or (action==0 and subAction=='f'):

        print '*'*5+' 23: ENSEMBL PROBESET ASSIGNATION for '+myChipName
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe_step22.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:            
            log.write('ACTION 22 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 22 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()
            
        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step6.bkdb"        
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probe.bkdb")            
            log.write('ACTION 6 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 6 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()                                
                    
        t1=time.time()
        
        Ok=ps_probeset.assign_probeset(mySpeciesName,myChipName,'ensembl',log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb")
            except:
                pass
            log.write('ACTION 23 FAILED\n')
            print('ACTION 23 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()             
        if aceFlag:                    
            oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"        
            newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step23.bkdb"
            os.rename(oldPath, newPath)
        
        t2=time.time()
        log.write('%s\t%s\t\t23\t%s_probeset.bkdb\t\tensembl update\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))
        print 'processed in %u s'%int(round(t2-t1))        

        print '%s\n'%('*'*5)
        
                     
    """PROBESET ASSIGNATION ON ACEVIEW"""
                
    if action==24  or (action==0 and subAction=='g'):
        
        print '*'*5+' 24: ACEVIEW POSITION_MAPPING for '+myChipName 
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene_step16.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            log.write('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()

        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe_step6.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene.bkdb")
            log.write('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 21 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()     
               
        path=os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_positions_by_region.bkdb"                
        regionDB=bsddb.btopen(path)
        regions=regionDB.keys()                                                         
        t1=time.time()     
        
        for region in regions:
            Ok=ps_probeset.position_mapping(mySpeciesName,myChipName,'aceview',region,log=log)
            if Ok!=None:
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+mySpeciesName+"_exons_by_gene.bkdb")
                except:
                    pass
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb")
                except:
                    pass
                log.write('ACTION 24 FAILED FOR REGION %s\n'%region)
                print('ACTION 24 FAILED FOR REGION %s\n'%region)
                print '%s\n'%('*'*5)
                sys.exit()                             

        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb"
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe_step24.bkdb"
        os.rename(oldPath, newPath)                 
                      
        t2=time.time()
        log.write('%s\t%s\t\t24\t%s_exons_by_gene.bkdb\t\taceview update\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))
        log.write('%s\t%s\t%s\t24\t%s_probe.bkdb\t\taceview update\t%.2f\n'%(date.today(),mySpeciesName,myChipName,myChipName,t2-t1))        

        print '%s\n'%('*'*5)
    
    if action==25 or (action==0 and subAction=='g'):
                                        
        print '*'*5+' 25: ACEVIEW PROBESET ASSIGNATION for '+myChipName
        
        oldPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe_step24.bkdb"        
        newPath=os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:            
            log.write('ACTION 25 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 25 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()
        
        oldPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset_step23.bkdb"        
        newPath=os.environ['PYDATA']+"/" + mySpeciesName + "/" + myChipName + "_probeset.bkdb"
        try:
            shutil.copyfile(oldPath,newPath)            
        except:
            os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb")            
            log.write('ACTION 25 FAILED (%s does not exist)\n'%oldPath)
            print('ACTION 25 FAILED (%s does not exist)\n'%oldPath)
            print '%s\n'%('*'*5)
            sys.exit()            

        t1=time.time()
        
        Ok=ps_probeset.assign_probeset(mySpeciesName,myChipName,'aceview',aceVersion=aceVersion,log=log)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probe.bkdb")
            except:
                pass
            log.write('ACTION 25 FAILED\n')
            print('ACTION 25 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                                     
    
        t2=time.time()
        log.write('%s\t%s\t\t25\t%s_probeset.bkdb\t\taceview update\t%.2f\n'%(date.today(),mySpeciesName,mySpeciesName,t2-t1))        
        print 'processed in %u s'%int(round(t2-t1))
                
        print '%s\n'%('*'*5)
        
                                
    """WRITE RESULTS"""
                                
    if action==26 or (action==0 and subAction=='h'):
            
        print '*'*5+' 26: MAKE PROBESET LIST for '+myChipName                
        t1=time.time()
        
        Ok=ps_probeset.make_probeset_list(mySpeciesName,myChipName)
        
        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/" + mySpeciesName+"/" + myChipName+"_probeset.dump")
            except:
                pass
            log.write('ACTION 26 FAILED\n')
            print('ACTION 26 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                                     
        
        t2=time.time()
        log.write('%s\t%s\t%s\t26\t%s_probeset.dump\t\t\t%.2f\n'%(date.today(),mySpeciesName,myChipName,myChipName,t2-t1))        
        print '%s\n'%('*'*5)
        
    
    if action==27 or (action==0 and subAction=='h'):
        
        print '*'*5+' 27: PROBESETS BY GENE for '+myChipName
                                        
        Ok=ps_probeset.probesets_by_gene(mySpeciesName,myChipName,aceFlag,log=log)

        if Ok!=None:
            try:
                os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/ensembl/"+myChipName+"_probesets_by_gene.bkdb")
            except:
                pass
            if aceFlag:
                try:
                    os.remove(os.environ['PYDATA']+"/"+mySpeciesName+"/aceview/"+myChipName+"_probesets_by_gene.bkdb")
                except:
                    pass            
                
            log.write('ACTION 27 FAILED\n')
            print('ACTION 27 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                                             
                
        print '%s\n'%('*'*5)
                
    if action==28 or (action==0 and subAction=='h'):
                             
        print '*'*5+' 28: WRITE PROBESETS BY GENE for '+myChipName
                         
        Ok=ps_export.write_probesets_by_gene(mySpeciesName,myChipName,log=log)
        
        if Ok!=None:
            log.write('ACTION 28 FAILED\n')
            print('ACTION 28 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()                                     
        
        print '%s\n'%('*'*5)
        
    if action==29 or (action==0 and subAction=='h'):
        
        print '*'*5+' 29: WRITE PROBESET LIST for '+myChipName
                    
        Ok=ps_export.write_probeset_list(mySpeciesName,myChipName,log=log)
         
        if Ok!=None:
            log.write('ACTION 29 FAILED\n')
            print('ACTION 29 FAILED\n')
            print '%s\n'%('*'*5)
            sys.exit()
                                                 
        print '%s\n'%('*'*5)
        
    
    if action==30 or (action==0 and subAction=='h'):
        if aceFlag:
        
            print '*'*5+' 30: WRITE ACEVIEW GENES for '+mySpeciesName
                            
            Ok=ps_export.write_aceview_genes(mySpeciesName,log=log)
            
            if Ok!=None:
                log.write('ACTION 30 FAILED\n')
                print('ACTION 30 FAILED\n')
                print '%s\n'%('*'*5)
                sys.exit()                                     
                    
            print '%s\n'%('*'*5)
        else:       
            
            print '*'*5+' 30: WRITE ENSEMBL GENES for '+mySpeciesName
                            
            Ok=ps_export.write_ensembl_genes(mySpeciesName,log=log)
            
            if Ok!=None:
                log.write('ACTION 30 FAILED\n')
                print('ACTION 30 FAILED\n')
                print '%s\n'%('*'*5)
                sys.exit()                                     
                    
            print '%s\n'%('*'*5)
     
    log.close()   