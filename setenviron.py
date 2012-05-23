#!/usr/bin/python
"""**Set environment variables containing data paths according to the host used.**"""

import os
from socket import gethostname

def setenviron():
    hostName=gethostname()
    #print 'host name: %s\n'%hostName
    jade=1    
    if hostName=='Cybele':
        rootName='E:/sosma/data/psawn'
        jade=0        
    elif hostName=='gauss':
        rootName='/home/mbellis/sosma/data/psawn'
        jade=0        
    elif hostName=='TRANSPORTABLE':
        rootName='E:/sosma/data/psawn'
        jade=0  
    if jade:
        rootName='/work/cinbell/data/psawn'
        os.environ['PROG']='/work/cinbell/projet/python/sosma/psawn'
        os.environ['SCRATCH']='/scratch/cinbell/psawn'        
                   
                    
    os.environ['AFFYDATA']=rootName+'/affydata'
    os.environ['ACEDATA']=rootName+'/acedata'
    os.environ['ENSDATA']=rootName+'/ensdata'
    os.environ['COMMONDB']=rootName+'/pydata/common'
    os.environ['PYDATA']=rootName+'/pydata'        
    os.environ['RAWDATA']=rootName+'/rawdata'
                        
if '__name__'=='__main__':
    setenviron()    