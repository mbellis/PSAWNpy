#!/usr/bin/python
"""**Functions for importing Ensembl data.**"""

import sys
from dbtools import *
import bintools
from datetime import date
from numpy import array
import ps_class
try:
    import MySQLdb
except:
    pass    
import setenviron
setenviron.setenviron()



    
def transcripts_by_exon(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find transcript(s) targeted by exons of Ensembl genes of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address    
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_transcripts_by_exons.bkdb
                key ==== Ensembl exon stable id
                value == TranscriptList object  
            
    """
              
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 7 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 7 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)        
                 
    cursor = connection.cursor()
      
    query = """
    SELECT DISTINCT
    exon_stable_id.stable_id,
    transcript_stable_id.stable_id
    FROM
    exon_stable_id, transcript_stable_id, exon_transcript
    WHERE
    exon_transcript.exon_id=exon_stable_id.exon_id
    AND exon_transcript.transcript_id=transcript_stable_id.transcript_id
    ORDER BY exon_stable_id.stable_id,transcript_stable_id.stable_id
    """
    
    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '7: lost MySQL connection'
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
            
        # create berkeley database                                                        
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_transcripts_by_exon.bkdb"
        if os.path.exists(path):
            os.remove(path)
            
        objClass='TranscriptList'
        objItems='idList=item[1]'    
        write_obj_db(lines, objClass, objItems,'list', path,1,log=log)                     

def exons_by_gene(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find exon(s) of Ensembl genes of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_exons_by_gene.bkdb
                key ==== Ensembl gene stable id
                value == ExonList
            
    """        

    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 8 ENSEMBL CONNECTION FAILED\n')
        print "ACTION 8 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])            
        sys.exit(1)        
             
    cursor = connection.cursor()
    
    query = """
    SELECT DISTINCT
    gene_stable_id.stable_id,
    exon_stable_id.stable_id,
    exon.seq_region_start,
    exon.seq_region_end,
    exon.seq_region_strand
    FROM
    seq_region, gene, gene_stable_id, transcript, exon_transcript, exon,exon_stable_id
    WHERE                    
    seq_region.seq_region_id=gene.seq_region_id
    AND transcript.seq_region_id=seq_region.seq_region_id    
    AND gene.gene_id=gene_stable_id.gene_id 
    AND transcript.gene_id=gene.gene_id     
    AND exon_transcript.transcript_id=transcript.transcript_id
    ANd exon.exon_id=exon_transcript.exon_id
    AND exon_stable_id.exon_id=exon.exon_id    
    ORDER BY seq_region.seq_region_id,gene_stable_id.stable_id,seq_region_start,seq_region_end
    """
    
    tryQuery=1
    t1=time.time()
    done=0
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '8: lost MySQL connection'
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
        
        #open transcripts_by_exon
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_transcripts_by_exon.bkdb"
        exonDB=bsddb.btopen(path)
        
        #create berkeley database        
        objClass='ExonList'    
        objItems='idList=item[1],exonStartArray=item[2], exonEndArray=item[3], strandArray=item[4]'
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_exons_by_gene.bkdb"
        write_obj_db(lines, objClass, objItems,'list', path,1,log=log)
        
        
        
        #remove doublons
        geneDB=bsddb.btopen(path)
        oldExonNb=0
        newExonNb=0        
        for gene in geneDB.keys():        
                currGene=cPickle.loads(geneDB[gene])            
                oldExonNb+=len(currGene.IDs)
                transcriptIDs=[]
                for i in range(len(currGene.IDs)):
                    #recover transcripts of the current exon
                    exons=cPickle.loads(exonDB[currGene.IDs[i]])                                
                    if len(exons.IDs)>1:
                        transcriptIDs.append(exons.IDs[0])
                        for j in range(1,len(exons.IDs)):
                            transcriptIDs.append(exons.IDs[j])
                            #create new start, ends, strand values
                            currGene.exonStarts.append(currGene.exonStarts[i])
                            currGene.exonEnds.append(currGene.exonEnds[i])
                            currGene.strands.append(currGene.strands[i])          
                            currGene.IDs.append(currGene.IDs[i])
                    else:                         
                        transcriptIDs.append(exons.IDs[0])                                                
                [exonIDList,exonStarts, exonEnds, exonStrands, transcriptsByExon, transcriptIDs, transcriptStarts, transcriptEnds, exonIndexes, exonGroups, intronStarts, intronEnds]=\
                        bintools.process_exons(exonIDs=array(currGene.IDs),exonStarts=array(currGene.exonStarts),exonEnds=array(currGene.exonEnds),exonStrands=array(currGene.strands),transcriptIDs=array(transcriptIDs))
                currGene.indexes=exonIndexes        
                currGene.IDs=exonIDList
                currGene.exonStarts=exonStarts
                currGene.exonEnds=exonEnds
                currGene.strands=exonStrands
                currGene.groups=exonGroups
                currGene.intronStarts=intronStarts
                currGene.intronEnds=intronEnds    
                currGene.transcriptsByExon=transcriptsByExon         
                currGene.transcriptIDs=transcriptIDs            
                currGene.transcriptStarts=transcriptStarts
                currGene.transcriptEnds=transcriptEnds
                        
                newExonNb+=len(currGene.IDs)       
                geneDB[gene]=cPickle.dumps(currGene,protocol=-1)
        geneDB.close()            
        if log!=0:
            log.write('%s\t%s\t\t8\t%s_exons_by_gene.bkdb\t\t%u\t%u\n'%(date.today(),species,species,oldExonNb,newExonNb))                      
    
    

def genes_by_region(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find genes(s) in Ensembl regions of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_genes_by_region.bkdb
                key ==== Ensembl stable region id
                value == GeneList
            
    """        


    # Connect to the ensembl database trough mySQL
    
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 9 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 9 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])
        
        sys.exit (1)
    
    # Create a mySQL cursor   
    cursor = connection.cursor()
      
    query = """
    SELECT DISTINCT
    gene.seq_region_id,
    gene_stable_id.stable_id,
    gene.seq_region_start,
    gene.seq_region_end,
    gene.seq_region_strand 
    from
    gene, gene_stable_id, gene_attrib
    where     
    gene.gene_id=gene_stable_id.gene_id
    and gene_attrib.gene_id=gene.gene_id
    and gene_attrib.attrib_type_id= "142"
    order by seq_region_id,seq_region_start
    """
    # attrib_type_id 142 corresponds to GC content of the gene
    #  no other attribute are available (homo_sapiens_core_49_36k)
    #test with : gene.seq_region_id='13' and
                
    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '9: lost MySQL connection'
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
            
        #create berkeley database                                                    
        objClass='GeneList'
        objItems='idList=item[1],startArray=array(item[2]), endArray=array(item[3]), strandArray=array(item[4])'
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_genes_by_region.bkdb"
        if os.path.exists(path):
            os.remove(path)
        write_obj_db(lines, objClass, objItems,'list', path,1,log=log)
        #create region_by_gene                
        regionDB=bsddb.btopen(path)
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_region_by_gene.bkdb"
        if os.path.exists(path):
            os.remove(path)
        geneDB=bsddb.btopen(path)
        for region in regionDB.keys():
            genes=cPickle.loads(regionDB[region])
            for geneID in genes.IDs:
                geneDB[geneID]=region
        geneDB.close()
        print '\n'     
            
    
    
   
def exons_by_predicted_transcript(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find exon(s) in Ensembl predicted transcripts of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_exons_by_predicted_transcripts.bkdb
                key ==== Ensembl predicted transcript id
                value == ExonList
            
    """
            
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 10 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 10 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)
     
    # Create a mySQL cursor      
    cursor = connection.cursor()
    
    query = """
    SELECT DISTINCT
    prediction_transcript_id,
    prediction_exon_id,
    seq_region_start, 
    seq_region_end,
    seq_region_strand 
    FROM
    prediction_exon
    ORDER BY prediction_transcript_id,seq_region_start
    """
    
    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '10: lost MySQL connection'
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
            
        # create berkeley database                                                    
        objClass='ExonList'
        objItems='idList=item[1],exonStartArray=array(item[2]), exonEndArray=array(item[3]), strandArray=array(item[4])'
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_exons_by_predicted_transcript.bkdb"
        if os.path.exists(path):
            os.remove(path)
        write_obj_db(lines, objClass, objItems,'list', path,1,log=log)
    
def predicted_transcripts_by_region(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find predicted transcripts by Ensembl region of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_predicted_transcripts_by_region.bkdb
                key ==== Ensembl region id
                value == TranscriptList
            
    """
    
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 11 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 11 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)
        
    # Create a mySQL cursor       
    cursor = connection.cursor()
      
    query = """
    SELECT DISTINCT
    prediction_transcript.seq_region_id, 
    prediction_transcript.prediction_transcript_id,
    prediction_transcript.seq_region_start,
    prediction_transcript.seq_region_end,
    prediction_transcript.seq_region_strand  
    FROM
    prediction_transcript    
    ORDER by seq_region_id,seq_region_start
    """
    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '11: lost MySQL connection'
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
            
        # create berkeley database                                                    
        objClass='TranscriptList'
        objItems='idList=item[1],startArray=array(item[2]), endArray=array(item[3]), strandArray=array(item[4])'
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_predicted_transcripts_by_region.bkdb"
        if os.path.exists(path):
            os.remove(path)
        write_obj_db(lines, objClass, objItems, 'list',path,1,log=log)
    
    


def transcripts_by_gene(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find transcripts by Ensembl region of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_transcripts_by_region.bkdb
                key ==== Ensembl region id
                value == TranscriptList
            
    """    
    
    try:
        connection = MySQLdb.connect(host = host, port = port, user = user, db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 12 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 12 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)
                  
    cursor = connection.cursor()
      
    query = """
    SELECT DISTINCT
    gene_stable_id.stable_id,
    transcript_stable_id.stable_id
    FROM
    gene_stable_id, transcript_stable_id, transcript
    WHERE
    transcript.gene_id=gene_stable_id.gene_id
    AND
    transcript.transcript_id=transcript_stable_id.transcript_id
    ORDER BY gene_stable_id.stable_id,transcript_stable_id.stable_id
    """
    
    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '12: lost MySQL connection'
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
            
        # create berkeley database                                                    
        objClass='TranscriptList'
        objItems='idList=item[1]'        
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_transcripts_by_gene.bkdb"
        if os.path.exists(path):
            os.remove(path)
        write_obj_db(lines, objClass, objItems,'list', path,1,log=log)

def chromosome_by_region(ensTable,species,host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0):
    """
    **Find chromosomes by Ensembl region of a particular species.**
    ::
    
        Arguments:  
            ensTable -- Ensembl table to be interrogated by MySQL
            species --- species
            host ------ Ensembl database address
            port ------ Ensembl port
            user ------ user name
            pswd ------ user pasword
            log ------- handle of a log file for recording messages
        
        Output:
            Create %species_chromosome_by_region.bkdb
                key ==== Ensembl region id
                value == chomosome (str)
            Create %species_region_by_chromosome.bkdb
                key ==== chomosome (str)
                value == Ensembl region id            
                    
    """    
        
    try:
        connection = MySQLdb.connect(host='ensembldb.ensembl.org', port=port,  user='anonymous', db = ensTable)
    except MySQLdb.Error, e:
        if log!=0:
            log.write('ACTION 13 ENSEMBL CONNECTION FAILED\n')            
        print "ACTION 13 ENSEMBL CONNECTION FAILED: %d: %s" % (e.args[0], e.args[1])        
        sys.exit (1)                 
    cursor = connection.cursor()        
      
    query = '''
    SELECT DISTINCT
    seq_region.seq_region_id,
    seq_region.name 
    FROM 
    seq_region, coord_system
    WHERE seq_region.coord_system_id=coord_system.coord_system_id
    AND coord_system.attrib="default_version"
    AND coord_system.name="chromosome"
    ORDER BY seq_region_id ;    
    ''' 

    tryQuery=1
    done=0
    t1=time.time()
    while tryQuery:
        try:    
            cursor.execute (query)
            tryQuery=0
            done=1
        except:
            print '13: lost MySQL connection'
            tryQuery=raw_input('would you like to continue (enter 1) or to cancel (enter 0)')
            if tryQuery=='0':
                tryQuery=0
                if log !=0:
                    log.write('\tCANCELED\t\n')   
                    sys.exit()                 
            elif tryQuery!="1":
                print 'enter 1 or 0'                
    if done:
        t1=time.time()    
        cursor.execute (query)
        t2=time.time()
        print "mySQL database %s delivered results in %.2f s\n"%("ensembl_"+ensTable,t2-t1)
           
        lines = cursor.fetchall()
        
        cursor.close ()
        connection.close ()
            
        # create berkeley database                                                        
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_chromosome_by_region.bkdb"
        if os.path.exists(path):
            os.remove(path)    
        write_str_db(lines, path,1,log=log)
        regionDB=bsddb.btopen(path)
        regions=regionDB.keys()
        chromosomes=list()
        for region in regions:
            chromosomes.append(regionDB[region])        
        regionDB.close()
        
        path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_region_by_chromosome.bkdb"
        if os.path.exists(path):
            os.remove(path)
        chromosomeDB=bsddb.btopen(path,'w')
        for i in range(len(chromosomes)):        
            chromosomeDB[chromosomes[i]]=regions[i]
        chromosomeDB.close()                        
        

    
def transcript(tscriptFile,ncbi,species,log=0):
    """
    **Extract data from an Ensembl fasta cDNA file.**
    ::
    
        Arguments:  
            tscriptFile -- Ensembl fasta cDNA file
            ncbi --------- ncbi version
            species ------ species        
            log ---------- handle of a log file for recording messages
        
        Input:
            %species_region_by_chromosome.bkdb
            
        Output:
            Create %species_region_by_chromosome_log.txt
            Create %species_transcript_sequence.bkdb
                key ==== Ensembl transcript id
                value == Transcript object
                                            
    """    
        
    path=os.environ['PYDATA']+"/"+species+"/log/"+species+"_region_by_chromosome_log.txt"
    logFile=open(path,'w')
    #open regions_by_chromosome
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_region_by_chromosome.bkdb"
    chromosomeDB=bsddb.btopen(path)
    #open the cDNA file
    transcriptFile = open(tscriptFile,'r')        
    dict={}
    #not necessary
    #hust initailized not to have warning in pydev
    transcriptID=""
    geneID=""
    sequence=""
    chromosome="" 
    region=""
    start=0
    end=0
    strand=0
    transcriptNb=0
    while 1:
        line = transcriptFile.readline().split("\n")[0]
        if not line:
            if transcriptNb>0:                
                dict[transcriptID]=ps_class.Transcript(ID=transcriptID,geneID=geneID,chromosome=chromosome,region=region,start=start,end=end,strand=strand,sequence=sequence)
            break
        # get some informations 
        if line[0]==">":
            #record the previous transcript            
            if transcriptNb>0:                
                dict[transcriptID]=ps_class.Transcript(ID=transcriptID,geneID=geneID,chromosome=chromosome,region=region,start=start,end=end,strand=strand,sequence=sequence)                
            # Get gene stable name               
            tab = line.split("gene:")[1].split(" ")         
            geneID = tab[0]    
            infoTable = line.rstrip().split(">")[1].split(" ")
            transcriptNb+=1
            # Get transcriptID stable name
            transcriptID = infoTable[0]
            # Get chromosome informations
            if infoTable[2].count('chromosome:%s:'%ncbi)>0 or infoTable[2].count('supercontig::')>0:
                tab=infoTable[2].split(":")
                chromosome=tab[2]
                region=chromosomeDB[chromosome]
                start=tab[3]
                end=tab[4]
                strand=tab[5]
            else:
                logFile.write('no chromosome nor supercontig for transcript %s\n'%transcriptID)
                chromosome=None
                region=None
                start=None
                end=None
                strand=None
                                                                    
            sequence=""            
        else:
            # Construct sequence
            sequence+=line
    #record the last transcript        
    if transcriptNb>0:
        dict[transcriptID]=ps_class.Transcript(ID=transcriptID,geneID=geneID,chromosome=chromosome,region=region,start=start,end=end,strand=strand,sequence=sequence)        
    logFile.close()                
    chromosomeDB.close()        
    # create berkeley database
    path=os.environ['PYDATA']+"/"+species+"/ensembl/"+species+"_transcript_sequence.bkdb"
    if os.path.exists(path):
        os.remove(path)    
    makebt(dict=dict,path=path,newFlag=1,objFlag=1,log=log)    

                  
if __name__=='__main__':
    exons_by_gene('mmusculus_core_62.37o',host='ensembldb.ensembl.org',port=5306,user='anonymous',pswd='',log=0)    
    pass