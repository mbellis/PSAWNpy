'''
Created on 4 juin 2011

@author: mbellis
'''

from bsddb import db                   # the Berkeley db data base
import sys, os, time


# fork into 2 instances and collide
pid = os.fork()

# Part 1: Create database and insert 4 elements
#
filename = 'fruit'

# Get an instance of BerkeleyDB environment
dbEnv = db.DBEnv()
# Set the automatic deadlock detection mode.
# The detect parameter configures the deadlock detector. The deadlock detector
# will reject the lock request with the lowest priority. If multiple lock 
# requests have the lowest priority, then the detect parameter is used 
# to select which of those lock requests to reject
# DB_LOCK_YOUNGEST : Reject the lock request for the locker ID with the youngest lock
dbEnv.set_lk_detect(db.DB_LOCK_YOUNGEST)
dbEnv.open("/home/mbellis/sosma/data",db.DB_INIT_LOCK|db.DB_CREATE| db.DB_INIT_MPOOL)

# force a collision
time.sleep(2)

#Aquires a locker id, guaranteed to be unique across all threads and processes that have the DBEnv open
lockID = dbEnv.lock_id()
lockID1 = dbEnv.lock_id()
print lockID
print lockID1

print 'pid:%u'%pid
print 'lockID:%u'%lockID

fruitDB = db.DB(dbEnv)
# Create a database in file "fruit" with a Hash access method
#     There are also, B+tree and Recno access methods
fruitDB.open(filename, None, db.DB_HASH, db.DB_CREATE)

# Print version information
print '\t', pid, db.DB_VERSION_STRING

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
    except Exception,error:
        print error
        pass

print pid, "here"

# force another one
time.sleep(2)

# Insert new elements in database
fruitDB.put("apple","red")
fruitDB.put("orange","orange")
fruitDB.put("banana","yellow")
fruitDB.put("tomato","red")

#Release the lock
dbEnv.lock_put(lock)
print pid, "there"

# Close database
fruitDB.close()

# Part 2: Open database and write its contents out
#
fruitDB = db.DB(dbEnv)
# Open database
#    Access method: Hash
#    set isolation level to "dirty read (read uncommited)"
fruitDB.open(filename, None, db.DB_HASH, db.DB_DIRTY_READ)

# get database cursor and print out database content
cursor = fruitDB.cursor()
rec = cursor.first()
while rec:
    print rec
    rec = cursor.next()
fruitDB.close()

#FROM http://web.archiveorange.com/archive/v/5H3d1EfqQF52D0YGBAKl
#----------
#
#still unanswered, under what conditions do you change the second 
#argument of lock_get?  Is that a simple lock identifier within the 
#context of the lock ID so that you can lock different records?  Anyway, 
#this locking is close enough for what I need to do.



for chip in chips:
    chipdb=bsddb.btopen('/home/mbellis/sosma/data/psawn/pydata/mouse/ensembl/%s_gops_by_region.bkdb'%chip)
    print '\nchip %s\n'%chip
    for region in chipdb.keys():
        print '%s\t%u'%(region,len(cPickle.loads(chipdb[region]).IDs))
   
        
    
    
os.chdir('/home/mbellis/download')
fout=open('GPL6193.txt','w')
fid=open('GPL6193-132.txt','r')

rank=0
for line in fid.readlines():    
    if rank>16:        
        fout.write('%s\n'%line.split('\t')[0])
    rank+=1           
print '%s'%line.split('\t')[0]
        
fid=open('GPL6193.txt','r')
rank=-1
for line in fid.readlines():
    rank+=1
    pass
        print line        

    


