#!/usr/bin/python
"""**Write tutorial in Sphinx format.**"""
import os,time
import setenviron
setenviron.setenviron()

if __name__=='__main__':
    CHIP='m8'    
    CHIP_IDS="['m4','m5','m8','m26','m27','m48','m49','m50','m53','m54','m62','m65','m67','m70']"
    CHIP_RANGE=''
    SPECIES='mouse'    
    SPECIES_IDS="['mouse','human']"
    SPECIES_RANGE=''    
    PS_IDS="['1415670_at']"
    ENS_PS_ID='1415670_at'
    ACE_PS_ID='1415670_at'
    PS_RANGE=''
    PS_RANK=0    
    ENS_VER='62_37o'
    PROBE_IDS=['0683:0877']
    PROBE_RANGE=[1,10,100]
    ENS_EXON_IDS=['ENSMUSE00000738026']
    ENS_EXON_RANGE=''
    ENS_GENE_IDS="['ENSMUSG00000000049']"
    ENS_GENE_RANGE=''
    REGION='1'
    PRED_TSCRIPT_IDS="['18616','18608']"
    PRED_TSCRIPT_RANGE=''
    ENS_TSCRIPT_IDS="['ENSMUST00000000090']"
    END_TSCRIPT_RANGE=''
    ACE_VER='Sep07'
    ACE_GENE_IDS="['2210417D09Rik']"
    ACE_GENE_RANGE=''
    ACE_TSCRIPT_IDS="['0610007C21Rik.hSep07']"
    ACE_TSCRIPT_RANGE=''
    ACE_EXON_IDS="['0610006L08Rik.exon0']"
    ACE_EXON_RANGE=''
    GOP_IDS=['GOPMUSG00000000086']
    GOP_RANGE=''    
    GOP_ID='GOPMUSG00000000046'
    GOP_RANK=0
    
    t1=time.time()
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','w')    

    tutoFile.write('**PSAWNpy**\n')
    tutoFile.write('-----------\n\n')
    tutoFile.write('All the application can be run in command line by using psawn.py. At each step control can be made with ps_test.py.\n\n')
    tutoFile.write('Importation of user data\n')
    tutoFile.write('++++++++++++++++++++++++\n')
    tutoFile.write('\n- construct chip database (=> data/pydata/common/chip.bkdb)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 1\n\t\n')
    tutoFile.write('\tpython ps_test.py -a 1 -i "%s" -r "%s"\n\t\n'%(CHIP_IDS,CHIP_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')
    os.system('python  ps_test.py -a 1 -i "%s" -r "%s" -t 1'%(CHIP_IDS,CHIP_RANGE))
    
        
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- construct species database (=> data/pydata/common/species.bkdb)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 2\n\t\n')
    tutoFile.write('\tpython ps_test.py -a 2 -i "%s" -r "%s"\n\t\n'%(SPECIES_IDS,SPECIES_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')
    os.system('python ps_test.py -a 2 -i "%s" -r "%s" -t 1'%(SPECIES_IDS,SPECIES_RANGE))
    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- construct probeset databases (=> data/pydata/common/chip_probeset.bkdb & => data/pydata/mouse/ensembl/%s_probeset.bkdb)\n\n'%CHIP_IDS)
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 3 -i %s\n\t\n'%CHIP_IDS)
    tutoFile.write('\tpython ps_test.py -a 3 -c %s -i "%s" -r "%s"\n\t\n'%(CHIP,PS_IDS,PS_RANGE))    
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 3 -c %s -i "%s" -r "%s" -t 1'%(CHIP,PS_IDS,PS_RANGE))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\nImportation of species information\n')
    tutoFile.write('++++++++++++++++++++++++++++++++++\n\n')

    tutoFile.write('Importation of Ensembl data\n')
    tutoFile.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    tutoFile.write('\n- Construct probe database (=> data/pydata/%s/ensembl/%s_probes.bkdb)\n\n'%(SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 4 -c %s -e %s\n\t\n'%(CHIP,ENS_VER))
    tutoFile.write('\tpython ps_test.py -a 4 -c %s -i "%s"\n\t\n'%(CHIP,PROBE_IDS))    
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py --a 4 -c %s -i "%s" -t 1'%(CHIP,PROBE_IDS))


    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Fill probeset with Ensembl gene identificator(s) (=> data/pydata/%s/ensembl/%s_probeset.bkdbk)\n\n'%(SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 5 -c %s -e %s\n\t\n'%(CHIP,ENS_VER))
    tutoFile.write('\tpython ps_test.py -a 5 -c %s -i "%s" -r "%s"\n\t\n'%(CHIP,PS_IDS,PS_RANGE))
    tutoFile.close()    
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 5 -c %s -i "%s" -r "%s" -t 1'%(CHIP,PS_IDS,PS_RANGE))    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Affymetrix informations (=> data/pydata/mouse/ensembl/%s_probeset.bkdb & => data/pydata/mouse/ensembl/%s_probe.bkdb)\n\n'%(CHIP,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 6 -c %s\n\t\n'%CHIP)
    tutoFile.write('\tpython ps_test.py -a 6 -c %s -i %s\n'%(CHIP,ENS_PS_ID))
    tutoFile.write('\tor (because only one probe set is displayed) \n')
    tutoFile.write('\tpython ps_test.py -a 6 -c %s -r %s\n\t\n'%(CHIP,PS_RANK))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 6 -c %s -i "%s" -t 1'%(CHIP,ENS_PS_ID))
        
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl list of transcripts for each exon (=> data/pydata/%s/ensembl/%s_transcripts_by_exon.bkdb)\n\n'%(CHIP,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 7 -s %s -e %s\n\t\n'%(SPECIES,ENS_VER))
    tutoFile.write('\tpython ps_test.py -a 7 -s %s -i "%s" -r "%s"\n\t\n'%(SPECIES,ENS_EXON_IDS,ENS_EXON_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 7 -s %s -i "%s" -r "%s" -t 1'%(SPECIES,ENS_EXON_IDS,ENS_EXON_RANGE))

    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl list of exons by gene (=> data/pydata/%s/ensembl/%s_exons_by_gene.bkdb)\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 8 -s %s -e %s\n\t\n'%(SPECIES,ENS_VER))
    tutoFile.write('\tpython ps_test.py -a 8 -i "%s" -s %s -r "%s"\n\t\n'%(ENS_GENE_IDS,SPECIES,ENS_GENE_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 8 -i "%s" -s %s -r "%s" -t 1'%(ENS_GENE_IDS,SPECIES,ENS_GENE_RANGE))
    
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl list of regions (=> data/pydata/mouse/ensembl/genes_by_region.bkdb & data/pydata/mouse/ensembl/region_by_gene.bkdb)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 9 -s %s -e %s\n\t\n'%(SPECIES,ENS_VER))
    tutoFile.write('\tpython ps_test.py -a 9  -i "%s" -s %s -q %s -r "%s"\n\t\n'%(ENS_GENE_IDS,SPECIES,REGION,ENS_GENE_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 9  -i "%s" -s %s -q %s -r "%s" -t 1'%(ENS_GENE_IDS,SPECIES,REGION,ENS_GENE_RANGE))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl list exons by predicted transcrips (=> data/pydata/ensembl/mouse_exons_by_predicted_transcript.bkdb)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 10 -e %s -s %s\n\t\n'%(ENS_VER,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 10  -i "%s" -r "%s" -s %s\n\t\n'%(PRED_TSCRIPT_IDS,PRED_TSCRIPT_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 10  -i "%s" -r "%s" -s %s -t 1'%(PRED_TSCRIPT_IDS,PRED_TSCRIPT_RANGE,SPECIES))
    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl list of predicted transcripts by region (=> data/pydata/ensembl/%s_predicted_transcript_by_region.bkdb)\n\n'%SPECIES)
    tutoFile.write(':: \n\n')
    tutoFile.write('\tpython psawn.py -a 11 -e %s -s %s\n\t\n'%(ENS_VER,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 11 -i "%s" -q %s -r "%s" -s %s\n\t\n'%(PRED_TSCRIPT_IDS.replace("'",""),REGION,PRED_TSCRIPT_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 11 -i "%s" -q %s -r "%s" -s %s -t 1'%(PRED_TSCRIPT_IDS.replace("'",""),REGION,PRED_TSCRIPT_RANGE,SPECIES))
    
        

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl transcripts by gene (=> data/pydata/mouse/ensembl/transcripts_by_gene.bkdb)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 12 -s %s -e 62_37o\n\t\n'%(SPECIES))
    tutoFile.write('\tpython ps_test.py -a -i "%s" -q %s -r "%s" -s %s\n\t\n'%(ENS_GENE_IDS,REGION,ENS_GENE_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')    
    os.system('python ps_test.py -a 12 -i "%s" -q %s -r "%s" -s %s -t 1'%(ENS_GENE_IDS,REGION,ENS_GENE_RANGE,SPECIES))
    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl region by chromosome (=> data/pydata/ensembl/%s_chromosome_by_region.bkdb & data/pydata/ensembl/%s_region_by_chromosome.bkdb)\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 13 -e "%s" -s %s\n\t\n'%(ENS_VER,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 13 -s %s\n\t\n'%(SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 13 -s %s -t 1'%(SPECIES))
    
        

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import Ensembl transcript sequences (=> data/pydata/%s/ensembl/transcript_sequence.bkdb)\n\n'%SPECIES)
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 14 -e %s -s %s\n\t\n'%(ENS_VER,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 14 -i "%s" -s %s -r "%s"\n\t\n'%(ENS_TSCRIPT_IDS,END_TSCRIPT_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 14 -i "%s" -s %s -r "%s" -t 1'%(ENS_TSCRIPT_IDS,SPECIES,END_TSCRIPT_RANGE))
        

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\nImportation of AceView data\n')
    tutoFile.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')

    tutoFile.write('\n- Import AceView transcript sequences (=> data/pydata/%s/aceview/%s_transcript_sequence.bkdb)\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 15 -s %s -v %s\n\t\n'%(SPECIES,ACE_VER))
    tutoFile.write('\tpython ps_test.py -a 15 -i "%s" -r "%s" -s %s\n\t\n'%(ACE_TSCRIPT_IDS,ACE_TSCRIPT_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 15 -i "%s" -r "%s" -s %s -t 1'%(ACE_TSCRIPT_IDS,ACE_TSCRIPT_RANGE,SPECIES))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Import AceView transcipt information (=> data/pydata/%s/aceview/%s_exons_by_gene.bkdd & data/pydata/%s/aceview/%s_genes_by_ensembl_region.bkdd & data/pydata/%s/aceview/mouse_transcrits_by_exon.bkdd & data/pydata/%s/aceview/%s_transcrits_by_gene.bkdd)\n\n'%(SPECIES,SPECIES,SPECIES,SPECIES,SPECIES,SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 16 -s %s -v %s\n'%(SPECIES,ACE_VER))

    tutoFile.write('\ndata/pydata/mouse/aceview/%s_exons_by_gene.bkdd\n\n'%SPECIES)
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython ps_test.py -a 16 -b a -i "%s" -s %s -q %s -r "%s"\n\t\n'%(ACE_GENE_IDS,SPECIES,REGION,ACE_GENE_RANGE))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 16 -b a -i "%s" -s %s -r "%s" -t 1'%(ACE_GENE_IDS,SPECIES,ACE_GENE_RANGE))
        
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\ndata/pydata/%s/aceview/%s_genes_by_ensembl_region.bkdd\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython ps_test.py -a 16 -b b -i "%s" -q %s -r "%s" -s %s\n\t\n'%(ACE_GENE_IDS,REGION,ACE_GENE_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 16 -b b -i "%s" -q %s -r "%s" -s %s -t 1'%(ACE_GENE_IDS,REGION,ACE_GENE_RANGE,SPECIES))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\ndata/pydata/%s/aceview/%s_transcrits_sequence.bkdd\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython ps_test.py -a 16 -b c -i "%s" -r "%s" -s %s\n\t\n'%(ACE_TSCRIPT_IDS,ACE_TSCRIPT_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 16 -b c -i "%s" -r "%s" -s %s -t 1'%(ACE_TSCRIPT_IDS,ACE_TSCRIPT_RANGE,SPECIES))

    
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\ndata/pydata/%s/aceview/%s_transcrits_by_gene.bkdd\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython ps_test.py -a 16 -b d -i "%s" -r "%s" -s %s\n\t\n'%(ACE_GENE_IDS,ACE_GENE_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 16 -b d -i "%s" -r "%s" -s %s -t 1'%(ACE_GENE_IDS,ACE_GENE_RANGE,SPECIES))
    
    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')        
    tutoFile.write('\ndata/pydata/%s/aceview/%s_transcrits_by_exon.bkdd\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython ps_test.py -a 16 -b e -i "%s" -r "%s" -s %s\n\t\n'%(ACE_EXON_IDS,ACE_EXON_RANGE,SPECIES))    
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 16 -b e -i "%s" -r "%s" -s %s -t 1'%(ACE_EXON_IDS,ACE_EXON_RANGE,SPECIES))
    
    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')    
    tutoFile.write('\n- Find correspondance between AceView genes and Ensembl genes (data/pydata/mouse/aceview/ensembl_genes_by_gene.bkdd)\n\n')
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 17 -s %s\n\t\n'%(SPECIES))    
    tutoFile.write('\tpython ps_test.py -a 17 -i "%s" -r "%s" -s %s\n\t\n'%(ACE_GENE_IDS,ACE_GENE_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 17 -i "%s" -r "%s" -s %s -t 1'%(ACE_GENE_IDS,ACE_GENE_RANGE,SPECIES))
    

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\nProcess chip information\n')
    tutoFile.write('+++++++++++++++++++++++++\n')

    tutoFile.write('\n- Find probe positions by region for a chip list (data/pydata/%s/ensembl/%s_positions_by_region.bkdd)\n\n'%(SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 18 -e %s -l "%s" -s %s\n\t\n'%(ENS_VER,CHIP_IDS,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 18 -c %s -i "%s" -q %s -r "%s" -s %s\n\t\n'%(CHIP,PROBE_IDS,REGION,PROBE_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 18 -c %s -i "%s" -q %s -r "%s" -s %s -t 1'%(CHIP,PROBE_IDS,REGION,PROBE_RANGE,SPECIES))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Merge probe positions of different chips (data/pydata/%s/ensembl/%s_positions_by_region.bkdd)\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 19 -i "%s" -s %s\n\t\n'%(CHIP_IDS,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 19 -c %s -q %s -r "%s" -s %s\n\t\n'%(CHIP,REGION,PROBE_RANGE,SPECIES))        
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')       
    os.system('python ps_test.py -a 19 -c %s -q %s -r "%s" -s %s -t 1'%(CHIP,REGION,PROBE_RANGE,SPECIES))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\n- Find group of probes (GOP) (data/pydata/%s/ensembl/%s_gop_by_region.bkdd)\n\n'%(SPECIES,SPECIES))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 20 -s %s -v %s\n\t\n'%(SPECIES,ACE_VER))
    tutoFile.write('\tpython ps_test.py -a 20 -c %s -i "%s" -q %s -r "%s" -s %s\n\t\n'%(CHIP,GOP_IDS,REGION,GOP_RANGE,SPECIES))
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 20 -c %s -i "%s" -q %s -r "%s" -s %s -t 1'%(CHIP,GOP_IDS,REGION,GOP_RANGE,SPECIES))

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\nprocess Ensembl information\n')
    tutoFile.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    tutoFile.write('\n- Update probe information for probes targeting genes (update data/pydata/%s/ensembl/%s_exons_by_gene.bkdb & data/pydata/%s/ensembl/%s_probe.bkdb)\n'%(SPECIES,SPECIES,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 21 -b a -c %s -s %s\n\t\n'%(CHIP,SPECIES))    
            
    tutoFile.write('\n- Update probe information for probes targeting GOPs (update data/pydata/%s/ensembl/%s_exons_by_gene.bkdb & data/pydata/%s/ensembl/%s_probe.bkdb)\n'%(SPECIES,SPECIES,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 22 -b a -c %s -s %s\n\t\n'%(CHIP,SPECIES))        
                
    tutoFile.write('\n- Assign probe sets to ensembl genes (update data/pydata/%s/%s_probeset.bkdb & data/pydata/%s/ensembl/%s_probe.bkdb\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython psawn.py -a 23 -c %s -s %s\n\t\n'%(CHIP,SPECIES))
                    
    tutoFile.write('\nprocess AceView information\n')
    tutoFile.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    
    tutoFile.write('\n- Update probe information for probes targeting genes (update data/pydata/%s/aceview/%s_exons_by_gene.bkdb & data/pydata/%s/aceview/%s_probe.bkdb)\n'%(SPECIES,SPECIES,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 24-b a -c %s -s %s\n\t\n'%(CHIP,SPECIES))    

    tutoFile.write('\n- Assign probe sets to aceview genes (update data/pydata/%s/%s_probeset.bkdb & data/pydata/%s/aceview/%s_probe.bkdb\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython psawn.py -a 25 -c %s -s %s\n\t\n'%(CHIP,SPECIES))       
    
    tutoFile.write('\n- Updated probe and probeset files (data/pydata/%s/%s_probeset.bkdb & data/pydata/%s/ensembl/%s_probe.bkdb\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython psawn.py -a 23 -c %s -s %s\n\t\n'%(CHIP,SPECIES))
    tutoFile.write('\tpython ps_test.py -a 23 -b b -c %s -i %s\n\t\n'%(CHIP,ENS_PS_ID))   
    tutoFile.write('\tor (because only one probe set is displayed)\n')
    tutoFile.write('\tpython ps_test.py -a 23 -b b -c %s -r %s\n\t\n'%(CHIP,PS_RANK))     
    tutoFile.close()
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')   
    os.system('python ps_test.py -a 23 -b b -c %s -i %s -t 1'%(CHIP,ENS_PS_ID))                  

    os.chdir(os.environ['PYDATA'])
    tutoFile=open('tuto.txt','a')
    tutoFile.write('\nExport data\n')
    tutoFile.write('++++++++++++\n')
    tutoFile.write('\n- Make a dump of %s_probeset.bkdb (pydata/%s/%s_probeset.dump)\n'%(CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 26 -c %s -s %s \n\t\n'%(CHIP,SPECIES))

    tutoFile.write('\n- Update probesets by gene files (pydata/%s/ensembl/%s_probesets_by_gene.bkdb & pydata/%s/aceview/%s_probesets_by_gene.bkdb)\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 27-c %s -s %s \n\t\n'%(CHIP,SPECIES))
    
    tutoFile.write('\n- Write probesets by gene files (pydata/%s/txt/ensembl_%s_probesets_by_gene_xx.txt & pydata/%s/txt/aceview_%s_probesets_by_gene_xx.txt)\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')    
    tutoFile.write('\tpython psawn.py -a 28-c %s -s %s \n\t\n'%(CHIP,SPECIES))
      
    tutoFile.write('\n- Write probeset lists (pydata/%s/txt/%s_probesets_ensembl.txt & pydata/%s/txt/%s_probesets_aceview.txt)\n'%(SPECIES,CHIP,SPECIES,CHIP))
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 29-c %s -s %s \n\t\n'%(CHIP,SPECIES))
    
    tutoFile.write('\n- Write AceView genes (pydata/%s/txt/%s_aceview_genes.txt\n\n'%(SPECIES,SPECIES)) 
    tutoFile.write('::\n\n')
    tutoFile.write('\tpython psawn.py -a 30 -s %s \n\t\n'%SPECIES)
    tutoFile.close()


    t2=time.time()
    print 'tuto made in %.0f s'%(t2-t1)
    