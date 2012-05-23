#!/usr/bin/python
import os
import setenviron
setenviron.setenviron()

if __name__=='__main__':
    CHIP='m5'            
    GOP_CHIP_IDS="['m4','m5','m8','m26','m27','m48','m49','m50','m53','m54','m62','m65','m67','m70']"    
    SPECIES='mouse'   
    ENS_VER='62_37o'
    ACE_VER='Sep07' 
    ACE_FLAG=1 
    
#    CHIP='m3'            
#    GOP_CHIP_IDS="['m2','m3','m21','m22','m24','m43','m44','m45','m46','m47','m63','m64','m66','m69']"        
#    SPECIES='human'   
#    ENS_VER='63_37'
#    ACE_VER='aug10'
#    ACE_FLAG=1
#    ADD_PROBESET=0
    
#    CHIP='m6'            
#    GOP_CHIP_IDS="['m6','m31','m32','m33','m51','m52','m56','m68','m74']"
#    #GOP_CHIP_IDS="['m68']"        
#    SPECIES='rat'   
#    ENS_VER='63_34'
#    ACE_VER=''
#    ACE_FLAG=0    
#    
    os.chdir('/home/mbellis/sosma/projet/sosma/src/psawn')        
    make=open('make_%s'%CHIP,'w')
#    make.write('python  psawn.py -a 0 -b a -l "%s" -s %s\n'%(GOP_CHIP_IDS,SPECIES))
#    make.write('python  psawn.py -a 0 -b b -c %s -e %s -f %u -p %u\n'%(CHIP,ENS_VER,ACE_FLAG,ADD_PROBESET))
#    make.write('python  psawn.py -a 0 -b c -e %s -s %s\n'%(ENS_VER,SPECIES))
#    make.write('python  psawn.py -a 0 -b d -s %s -v %s\n'%(SPECIES,ACE_VER))    
    make.write('python  psawn.py -a 0 -b f -c %s -f %u\n'%(CHIP,ACE_FLAG))
    make.write('python  psawn.py -a 0 -b g -c %s -v %s\n'%(CHIP,ACE_VER))
#    make.write('python  psawn.py -a 0 -b e -e %s -f %u -l "%s" -s %s\n'%(ENS_VER,ACE_FLAG,GOP_CHIP_IDS,SPECIES))    
    make.write('python  psawn.py -a 0 -b h -c %s -f %u\n'%(CHIP,ACE_FLAG))
    make.close()