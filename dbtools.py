#!/usr/bin/python
"""**Functions for creating Berkeley data bases.**"""
import os, time, cPickle
import bsddb
import ps_class
from numpy import array
    
def write_obj_db(lines, objClass, objItems,objType, path, newFlag,log=0):
    """
    **Prepare data for creating data base with objects as values.**
    ::
    
        Arguments:
        lines ----- is a list of tuples (line) which is ordered on line[0] used as key
                    to construct a dictionary. Several lines may have the same key. In this
                    case, class properties are lists.             
        objClass -- class of the object stored as a  value in the constructed dictionary             
        objItems -- a string which is evaluated to assign values stored in a dictionary
                    (item[1],item[2],...) to the right object property
        path ------ name of the btree berkeley data base constructed with the keys and values lists
        newFlag  -- indicates if the data base is new (newFlag=1 and in this case if it exists,
                    it is erased) or is to be updated (newFlag=0)
        log ------- handle of a log file for recording messages             
    
    """        
    
    keys=[]
    values=[]
    item={}        
    
    currKey=''
    keyNb=0
    for line in lines: #                             
        if str(line[0])!=currKey:
            #write previous key in keys list
            # and previous object in value list if they exist
            if keyNb>0:               
                keys.append(currKey)                                                           
                eval('values.append(ps_class.%s(%s))'%(objClass,objItems))                                                            
            #initialize item lists
            keyNb+=1
            currKey=str(line[0])                            
            for i in range(0,len(line)):                    
                if objType=='list':
                    item[i]=list()
                    item[i].append(line[i])
                elif objType=='string':                            
                    item[i]=line[i]                    
        else:
            for i in range(1,len(line)):
                if objType=='list':                
                    item[i].append(line[i])
                elif objType=='string':
                    print 'error string object %s cannot be appended'%item[i]                                                                  
                
    #write last key in keys list
    # and last object in value list if they exist    
    if keyNb>0:
        keys.append(currKey)                       
        eval('values.append(ps_class.%s(%s))'%(objClass,objItems))        
               
    makebt(keys=keys,values=values,path=path,newFlag=newFlag,objFlag=1,log=log)
    

def write_str_db(lines, path, newFlag,log=0):
    """
    **Prepare data for creating data base with strings as values.**
    ::
    
        Arguments:
        lines ----- is a list of tuples (line) which is ordered on line[0] used as key
                    to construct a dictionary. All lines have different keys.
        path ------ name of the btree berkeley data base constructed with the keys and values lists
        newFlag  -- indicates if the data base is new (newFlag=1 and in this case if it exists,
                    it is erased) or is to be updated (newFlag=0)
        log ------- handle of a log file for recording messages            
    
    """        
    
    keys=[]
    values=[]    
    for line in lines:                    
        keys.append(str(line[0]))       
        values.append(line[1])                     
    makebt(keys=keys,values=values,path=path,newFlag=newFlag,objFlag=0,log=log)

def makebt(dict=None,keys=None,values=None,path=None,newFlag=0,objFlag=0,log=0):
    """
    **Create a btree Berkeley database.**
    ::
    
        Arguments:
            dict ----- dictionary to be saved
            keys ----- list of keys
            values --- list of values
            path ----- path of the berkeley database file to be created            
            newFlag -- indicates if the data base is new (newFlag=1 and in this case if it exists,
                       it is erased) or is to be updated (newFlag=0)
            objFlag -- indicate if the value must be stored as a string (0) or
                       as an object by using cPickle
            log ------- handle of a log file for recording messages
                
    """
    
    #either dict or keys and values must be different of None
    if dict is None:   
        assert keys is not None, "keys must be passed as parameter in makebt"
        assert values is not None, "values must be passed as parameter in makebt"
    if keys is None:
        assert dict is not None, "dict must be passed as parameter in makebt"    
    assert path is not None, "path must be passed as parameter in makebt"
    assert objFlag in [0,1], "objFlag must be equal to 0 or 1"    
    
    t1=time.time()
    if os.path.exists(path):        
        if newFlag==0:
            print 'open DB '+ os.path.basename(path)
            #berkDB.open(path, None)
            berkDB=bsddb.btopen(path,'w')
        else:            
            print 'recreate DB '+ os.path.basename(path)
            #berkDB.open(path, None, db.DB_BTREE, db.DB_CREATE)
            os.remove(path)            
            berkDB=bsddb.btopen(path,'w')            
    else: 
        dirName=os.path.dirname(path)
        if not os.path.exists(dirName):
            os.makedirs(dirName)                  
        print 'create DB '+ os.path.basename(path)
        #berkDB.open(path, None, db.DB_BTREE, db.DB_CREATE)
        berkDB=bsddb.btopen(path,'w')

# Insert new elements in database    
    if dict is None:
        if objFlag==0:
            for i in range(len(keys)):                       
                berkDB[keys[i]]=str(values[i])
        else:            
            for i in range(len(keys)):                           
                berkDB[keys[i]]=cPickle.dumps(values[i],protocol=-1)
    else:
        if objFlag==0:                
            for key in dict.keys():                                        
                berkDB[key]=str(dict[key])
        else:                            
            for key in dict.keys():                                                
                berkDB[key]=cPickle.dumps(dict[key],protocol=-1)
# Close database 
    berkDB.close()
    t2=time.time()
    if keys is not None:
        print '%u records created in %.2f s'%(len(keys),t2-t1)
        if log!=0:
            log.write('\t%u\t%.2f\n'%(len(keys),t2-t1))
    else:
        print '%u records created in %.2f s'%(len(dict.keys()),t2-t1)
        if log!=0:
            log.write('\t%u\t%.2f\n'%(len(dict.keys()),t2-t1))        