#!/usr/bin/python

from string import *
from UserList import *
from ConfigParser import *
from vasp.parsing.Structure import *

    
class DBI:
  def __init__(self,db=None,name=None,connect=1):
    self.db=db
    self.blocksize=10
    self.subst={}
    self.cfg={}
    self.name=name
    if connect:
      self.connect()

    if hasattr(self,"init"):
      self.init()
    
    
  def assure(self,key):
    if self.subst.has_key(key):
      return 1
    raise "The %s is missing."%(key[1:])

  def config(self,c,name=None):
    if name is not None:
      self.name=name
      
    self.cfg.update(c)
    for k,v in self.cfg.items():
      if lower(k) in ["user_uid","username"]:
        self.subst["#"+upper(k)]='"%s"'%v
	continue
      key=None
      if lower(k[:7])=="define_":
        key=upper(k[7:])
      elif lower(k[:6])=="table_":
        key=upper(k[6:])
      if key is not None:
	if "#" not in key:	
          self.subst["#"+key]=v
                
	
  def connect(self):
    pass
  
  def disconnect(self):
    self.db=None
    self.user_id=None
    
  def isConnected(self):
    return self.db is not None
    
  def subs(self,s):
    keys=self.subst.keys()
    keys.sort(key=len,reverse=1)
    for k in keys:
      s=replace(s,k,self.subst[k])
    print s
    return s

  def exe(self,s,ignore_errors=0):
    try:
      c=self.db.cursor()
      c.execute(self.subs(s))
    except:
      if not ignore_errors:
        raise
      

  def fetchall(self,s):
    c=self.db.cursor()
    c.execute(self.subs(s))
    return c.fetchall()

  def fetchone(self,s):
    c=self.db.cursor()
    c.execute(self.subs(s))
    return c.fetchone()

  def fetchiter(self,s):
    c=self.db.cursor()
    c.execute(self.subs(s))
    while 1:
      d=c.fetchone()
      if d is not None:
        yield d
      else:
        break

  def fetchvalue(self,s):
    return self.fetchone(s)[0]
  def fetchvalues(self,s):
    return map(lambda x:x[0],self.fetchall(s))

  def run(self,s,ignore_errors=0):
    for ss in split(s,";"):
      self.exe(ss,ignore_errors=ignore_errors)

  def runfile(self,f):
    if type(f)==type(""):
      f=open(f,"r")
    s=f.read()
    f.close()
    self.run(s)

  def insertRecord(self,table,indexname,keys,values):
    self.exe("INSERT INTO %s (%s) VALUES (%s)"%(table,keys,values))
    return self.fetchvalue("SELECT MAX(%s) FROM %s"%(indexname,table))

  def commit(self):
    self.db.commit()
    
      

class Query:
  def __init__(self,container,data,condition,cachesize=200):
    self.container  = container
    self.data       = data
    self.condition  = condition
    self.cachesize  = cachesize
    self.refresh()

  def __len__(self):
    return self.cumulative[-1]

  def refresh(self):
    self.cache      = []
    self.hash       = {}
    self.count()
    
  def finddb(self,i):
    for j in xrange(len(self.cumulative)-1):
      if self.cumulative[j]<=i and self.cumulative[j+1]>i:
        return j	        
  
  def shrinkCache(self,n):
    "Remove (up to) n elements from the cache."
    n=min(n,len(self.cache))
    if n>0:
      del self.cache[:n]

      for k,v in self.hash.items():      
	if v<n:
	  del self.hash[k]
	else:
	  self.hash[k]=v-n

  def addToCache(self,i,r):
    v=len(self.cache)
    self.cache.append(r)
    self.hash[i]=v
    
  def fetchblock(self,i):
    di=self.finddb(i)
    if di is not None:
      d=self.container[di]
      blocksize=d.blocksize
      offset=i-self.cumulative[di]
      self.shrinkCache(len(self.cache)+blocksize-self.cachesize)
      for x in d.fetchiter("SELECT %s %s LIMIT %d, %d"%(self.data,self.condition,offset,blocksize)):
        self.addToCache(i,(di,x))
	i+=1

  def get(self,i):
    if i>=0 and i<len(self):
      if self.hash.has_key(i):
#        print "cache hit",i,self.hash[i],self.cache[self.hash[i]]
	return self.cache[self.hash[i]]
      else:
#        print "not in cache -> fetch",i
	self.fetchblock(i)
	if self.hash.has_key(i):      
#          print "      got",i,self.hash[i],self.cache[self.hash[i]]
	  return self.cache[self.hash[i]]
	else:
#	  print "error -> refresh"
          self.refresh()
	  self.fetchblock(i)
	  if self.hash.has_key(i):      
#            print "      got",i,self.hash[i],self.cache[self.hash[i]]
	    return self.cache[self.hash[i]]
	  else:
#	    print "FAILED"
            return None	

  def getiter(self):
    return self.container.fetchiter("SELECT %s %s"%(self.data,self.condition))
        
  def __getitem__(self,i):
#    print "getitem",i,"/",len(self),":",self.get(i)
    return self.get(i)[1]  

  def __iter__(self,i):
    for x in self.getiter():
      yield x[1]
        
  def count(self):
    self.counts     = []
    self.cumulative = [0]
    s=0
    
    for d in self.container:
      l=d.fetchvalue("SELECT count(*) %s"%self.condition)
      s+=l
      print d.name,l,s
      self.counts.append(l)
      self.cumulative.append(s)      
      
  
class DBIContainer(UserList):
  def exe(self,s):
    for x in self:
      x.exe(s)

  def fetchone(self,s):
    for x in self:
      v=x.fetchone(s)
      if v is not None:
        return v
    return None	

  def fetchiter(self,s):
    for d in self:
      for r in d.fetchiter(s):
        yield r

  def query(self, data, condition="", cachesize=200):
    return Query(self,data,condition,cachesize)
    
  def fetchvalue(self,s):
    return self.fetchone(s)[0]

  def run(self,s,ignore_errors=0):
    for x in self:
      x.run(s,ignore_errors=ignore_errors)

  def commit(self):
    for x in self:
      x.commit()
      
  def connect(self):
    for x in self:
      x.connect()

  def disconnect(self):
    for x in self:
      x.disconnect()

  def isConnected(self):
    for x in self:
      if x.isConnected():
        return 1
    return 0

class SQLiteDBI(DBI):
  def __init__(self,db=None,connect=1):
    if type(db)==type(""):
      self.path=db
      DBI.__init__(self,connect=connect)
    else:
      DBI.__init__(self,db,connect=connect)

  def config(self,c,name=None):
    cc={}
    cc.update(c)
    c=cc
    DBI.config(self,c,name)
    if c.has_key("path"):
      if self.isConnected():
	self.disconnect()
	self.path=c["path"]
	self.connect()
      else:
	self.path=c["path"]
        
  def connect(self):
    from pysqlite2 import dbapi2 as sqlite
    if self.isConnected():
      self.disconnect()
    self.db=sqlite.connect(self.path)
  
  def disconnect(self):
    self.db=None


  

def von(x):
  "str(x) or NULL"  
  if x is None:
    return "NULL"
  else:
    return str(x)

def son(x):
  '"str(x)" or NULL'  
  if x is None:
    return "NULL"
  else:
    return '"%s"'%str(x)

class UserManagementMixin:
  def init(self):
    self.user_id=None
    
  def login(self):
    self.user_id=self.dbi.getUserID()
    return self.user_id

  def getUserName(self):
    return self.cfg["username"]

  def getUserUID(self):
    return self.cfg["user_uid"]
  
  def login(self):
    self.assure("#USER_UID")
    try:
      self.user_id=self.fetchvalue('SELECT id FROM #USERINFO WHERE user_uid=#USER_UID')
    except:
      self.user_id=None
    return self.user_id

  def addUser(self):
    self.assure("#USER_UID")
    self.assure("#USERNAME")
    if self.getUserID() is None:
      self.exe("INSERT INTO #USERINFO (user_uid,username) VALUES (#USER_UID, #USERNAME)")
      self.commit()    
      self.login()

  def getUserID(self):
    if self.user_id is None:
      self.login()
    return self.user_id
  
class P4VMixin(UserManagementMixin):
  def init(self,calc_id=None,pm=None):
    UserManagementMixin.init(self)
    self.calc_id=calc_id
    self.pm=pm
    
        
  def addRecord(self,pm=None):
    if pm is None:
      pm=self.pm
    else:
      self.pm=pm
    self.login()
    name    = pm.NAME
    if name is None:
      name=""
    user_id = self.user_id
    if user_id is None:
      raise "Can not login to %s"%self.dbi.name
    self.calc_id=self.insertRecord("#CALC","id","name,user_id",'"%s",%d'%(name,user_id))
    self.commit()    
    return self.calc_id
    
  def storePM(self,pm):
    self.addRecord(pm)
    s=pm.INITIAL_STRUCTURE
    if s is not None:
      self.storeStructure(s,-1)
    s=pm.FINAL_STRUCTURE
    if s is not None:
      self.storeStructure(s,-2)
    seq=pm.STRUCTURE_SEQUENCE_L
    if seq is not None:
      for i in range(len(seq)):
        self.storeStructure(seq[i],i)
    self.commit()    
            

  def removeStructure(self,i):
    if i is not None:
      self.exe("DELETE FROM #STRUCT WHERE id=%d"%(i))
      self.exe("DELETE FROM #STRUCTPOS WHERE structure_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTFORCE WHERE structure_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTVELOCITY WHERE structure_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTCONSTRAINTS WHERE structure_id=%d"%(i))
  def removeCalculation(self,i):
    if i is not None:    
      self.exe("DELETE FROM #CALC WHERE id=%d"%(i))
      self.exe("DELETE FROM #STRUCTPOS WHERE calc_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTFORCE WHERE calc_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTVELOCITY WHERE calc_id=%d"%(i))
      self.exe("DELETE FROM #STRUCTCONSTRAINTS WHERE calc_id=%d"%(i))

  def storeStructure(self,s,step=None):
    calc_id=self.calc_id
    if calc_id is not None and step is not None:
      for l in self.fetchall("SELECT id FROM #STRUCT WHERE calc_id=%d AND step=%d"%(calc_id,step)):      
	removeStructure(db,l[0])

    info=s.info

    if len(s.scaling)!=1:
      s.correctScaling()
      
    cmd ="%s, %s,"%(von(calc_id), von(step))
    cmd+="%14.12f,'%s',"%(    s.scaling[0],s.comment)
    cmd+="""%14.12f,%14.12f,%14.12f,
    %14.12f,%14.12f,%14.12f,
    %14.12f,%14.12f,%14.12f,%d"""%(
    s.basis[0][0],s.basis[0][1],s.basis[0][2],
    s.basis[1][0],s.basis[1][1],s.basis[1][2],
    s.basis[2][0],s.basis[2][1],s.basis[2][2],s.types)

    Id=self.insertRecord("#STRUCT","id","calc_id, step, scale,comment,a11,a12,a13,a21,a22,a23,a31,a32,a33,species",cmd)

    if s.isSelective():
      sel    = s.selective
      for i in range(len(sel)):
	self.exe("""INSERT INTO #STRUCTCONSTRAINTS (structure_id,calc_id,atomnumber,x,y,z)
	VALUES (%d,%s,%d,%d,%d,%d)"""%(Id,von(calc_id),i,sel[i][0],sel[i][1],sel[i][2]))
        

    sc=Structure()
    sc.setStructure(s)
    sc.setCartesian()

          
    for i in range(len(s)):
      si=s.speciesIndex(i)

      self.exe("""INSERT INTO #STRUCTPOS (structure_id,calc_id,atomnumber,specie,element,x,y,z)
      VALUES (%d,%s,%d,%d,'%s',%14.12f,%14.12f,%14.12f)"""%(Id,von(calc_id),i,si,info[si].element,
      sc[i][0],            sc[i][1],            sc[i][2]))
    return Id

  def readStructure(self,Id):
    s=Structure()
    l=self.fetchone("""SELECT calc_id,scale,comment,a11,a12,a13,a21,a22,a23,a31,a32,a33,
      species FROM #STRUCT WHERE id=%d"""%(Id))
    print len(l),l
    (calc_id,scale ,s.comment,
     s.basis[0][0],s.basis[0][1],s.basis[0][2],
     s.basis[1][0],s.basis[1][1],s.basis[1][2],
     s.basis[2][0],s.basis[2][1],s.basis[2][2],spec)=l
     
    s.scaling=[scale]
    s.setCartesian()

    l=self.fetchall("SELECT specie,element,x,y,z FROM #STRUCTPOS WHERE structure_id=%d ORDER BY atomnumber"%(Id))
    for spec,element,x,y,z in l:
      i=s.appendAtom(spec,Vector(x,y,z))
      #s.info[i].element=element

    l=self.fetchall("SELECT x,y,z FROM #STRUCTCONSTRAINTS WHERE structure_id=%d ORDER BY atomnumber"%(Id))
    if len(l):
      s.setSelective()
      for i in range(len(l)):
        s.selective[i]=[int(l[i][0]),int(l[i][1]),int(l[i][2])]

    return s

def createFromConfig(l,prototypes={"sqlite":SQLiteDBI},mixin=None):
  import os.path
  import getpass
  defaults={}
  defaults["home"]=os.path.expanduser("~")
  defaults["user"]=getpass.getuser()
  c=ConfigParser(defaults)
  c.read(l)
  c.write(open("out","w"))
  cont=DBIContainer()
  for s in c.sections():
#    try:
    t=lower(c.get(s,"type"))
#    except:
#      print "Section %s - type missing"%s
#      continue
#    try:
    DBI=prototypes[t]
    if mixin is not None:
      l=[mixin]
      l.extend(DBI.__bases__)
#        l=list(DBI.__bases__)
#	l.append(mixin)
      DBI.__bases__=tuple(l)

    d=DBI(connect=0)
      
#    except:
#      print "Section %s - unsupported type %s"%(s,t)
#      continue
#    try:
#      print "config",c
    d.config(c.items(s),s)
#    except:
#      print "Section %s - configuration failed"%(s)
#      continue
    cont.append(d)
  return cont    

  
if __name__=="__main__":    
  import vasp.parsing.SystemPM as pm
  d=createFromConfig(["tables.cfg"],mixin=P4VMixin)
  d.connect()
  cmd=open("tables1.sql").read()
  d.run(cmd,ignore_errors=1)
  
  d[0].addUser()
  s=pm.XMLSystemPM("vasprun.xml")
  d[0].storePM(s)
