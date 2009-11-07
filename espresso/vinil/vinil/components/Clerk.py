# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                        (C) 2007  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

maindom = 'vinil.dom'

# Contains list of classes
tablenames = ('Configuration',
              'User')
#tablenames = ('Simulation',
#              'Server',
#              'Matter',
#              'Job',
#              'Atom',
#              'Configuration',
#              'VTable',
#              'Session',
#              'User')
#              #'Authorization')

from pyre.components.Component import Component as base

class Clerk( base ):

    class Inventory( base.Inventory):

        import pyre.inventory
        db = pyre.inventory.str('db', default = 'vinil' )
        db.meta['tip'] = "the name of the database"

        dbwrapper = pyre.inventory.str(name='dbwrapper', default='psycopg2')
        dbwrapper.meta['tip'] = "the python package that provides access to the database back end"

    # parameter 'facility' doesn't make any sense here
    def __init__(self, name = 'clerk', facility = 'clerk'):
        base.__init__(self, name, facility)
        return

#    def __init__(self, *args, **kwds):
#        Component.__init__(self, *args, **kwds)

    def _init(self):
        base._init(self)

        # connect to the database
        import pyre.db
        dbkwds = DbAddressResolver().resolve(self.inventory.db)
        self.db = pyre.db.connect(wrapper=self.inventory.dbwrapper, **dbkwds)
        self.db.autocommit(True)    # ?

    def _configure(self):
        base._configure(self)
        self.db = self.inventory.db

    # Table specific methods
    def getJobs(self, id=None, where=None):
        '''retrieve job record specified by id'''
        return self._getEntry('Job', id=id, where=where)

    def getSimulations(self, id=None, where=None):
        '''retrieve simulation record specified by id'''
        return self._getEntry('Simulation', id=id, where=where)

    def getVTables(self, id=None, where=None):
        '''retrieve vtable record specified by id'''
        return self._getEntry('VTable', id=id, where=where)

    def getSessions(self, id=None, where=None):
        '''retrieve session record specified by id'''
        return self._getEntry('Session', id=id, where=where)

    def getServers(self, id=None, where=None):
        '''retrieve user server specified by id'''
        return self._getEntry('Server', id=id, where=where)

    def getConfigurations(self, id=None, where=None):
        '''retrieve user configuration specified by id'''
        return self._getEntry('Configuration', id=id, where=where)

    def getAtoms(self, id=None, where=None):
        '''retrieve atom specified by id'''
        return self._getEntry('Atom', id=id, where=where)

    def getMatters(self, id=None, where=None):
        '''retrieve matter specified by id'''
        return self._getEntry('Matter', id=id, where=where)

    # Not used at this point
    def getUsers(self, id=None, where=None):
        '''retrieve user record specified by id'''
        return self._getEntry('User', id=id, where=where)
    
    # Not used at this point
    def getAuthorizations(self, id=None, where=None):
        '''retrieve user authorization specified by id'''
        return self._getEntry('Authorization', id=id, where=where)


    def updateRecord(self, record):
        """Updates row in the database specified by record.id"""
        id          = record.id
        where       = "id='%s'" % id
        assignments = []

        # get the column names and couple them with the new values
        for column in record.getColumnNames():
            value = getattr( record, column )
            assignments.append( (column, value) )
            continue
            
        # update the row, or in other words, record
        self.db.updateRow(record.__class__, assignments, where)
        return record


    def insertRecord(self, record):
        """insert a new record into db"""
        try:
            self.db.insertRow( record )
        except:
#            columns = record.getColumnNames()
#            values = [ record.getColumnValue( column ) for column in columns ]
#            s = ','.join(
#                [ '%s=%s' % (column, value)
#                  for column, value in zip(columns, values)
#                  ] )
#            self._debug.log( 'failed to insert record: %s' % s)
            raise
        return record


    def deleteRecord(self, record, id=None):
        "Deletes a record"
        
        if id:
            self.db.deleteRow(record.__class__, where="id='%s'" % record.id)

        # else no effect?
        return


    # General methods working on database tables
    def getTable(self, tablename):
        for t in tablenames:
            table   = self._getClass(t)
            if table.name == tablename:
                return self._getAll(table)

        return None

    def _getClass(self, classname):
        module  = _import("%s.%s" % (maindom, classname))
        return getattr(module, classname)


    def _getEntry(self, classname, id=None, where=None):
        vclass = self._getClass(classname)
        if id is not None:
            return self._getRecordByID( vclass, id )

        return self._getAll(vclass, where)

    def setDb(self,dbName):
        self.db = dbName

    def _getAll(self, table, where = None):
        index = {}
        all = self.db.fetchall(table, where=where)
        return all

    def _getRecordByID(self, table, id ):
        all = self.db.fetchall( table, where = "id='%s'" % id )
        if len(all) == 1:
            return all[0]
        raise RuntimeError, "Cannot find record of id=%s in table %s" % (
            id, table.__name__)

    #def importTable(self, tablename):
        

def _import(package):
    return __import__(package, {}, {}, [''])

# Don't understand why I need this complicated class.
class DbAddressResolver:

    def resolve(self, address):
        tmp = address.split('@')
        if len(tmp)==1:
            svr = tmp[0]
            up = ''
        elif len(tmp)==2:
            up,svr = tmp
        else:
            raise ValueError, 'Invalid db address: %r' % address

        host,port,database = self._resolve_svr(svr)
        user, pw = self._resolve_up(up)
        ret = {
            'database': database,
            }
        if host: ret['host'] = host
        if port: ret['port'] = port
        if user: ret['user'] = user
        if pw: ret['password'] = pw
        return ret


    def _resolve_up(self, up):
        separator = ':'
        tmp = up.split(separator)
        if len(tmp) == 1:
            user = tmp[0]
            if user == '': user = None
            pw = None
        elif len(tmp) == 2:
            user, pw = tmp
        else:
            raise ValueError, 'Invalid user, password: %r' % up
        return user, pw


    def _resolve_svr(self, svr):
        separator = ':'

        if svr.find(separator) == -1:
            # unix socket
            return None, None, svr
        splits = svr.split(separator)
        if len(splits)==2:
            host, database = splits
            # default port: 5432
            return host, 5432, database
        elif len(splits)==3:
            host, port, database = splits
            return host, port, database
        raise ValueError, 'Invalid db svr: %r' % (svr,)



    """
    def indexActiveUsers(self):
        from jazzclub.dom.User import User
        index = {}
        users = self.getRecords(User)
        for key, user in users:
            index[user.username] = user
            continue
        return index


    def getUser(self, username):
        users = self.indexActiveUsers()
        return users[username]


    def getUserInfo(self, username):
        from jazzclub.dom.Registrant import Registrant
        registrants = self.getRecords(Registrant)

        candidates = filter(
            lambda r: r[1].username==username, registrants)

        assert len(candidates) < 2

        if not candidates: return
        return candidates[0][1]


    def newRecord(self, record):
        tablename = record.name
        store = self._retrieve_store( tablename )
        store[ record.id ] = record
        self._update_store( tablename, store )
        return


    def updateRecord(self, record):
        tablename = record.name
        store = self._retrieve_store( tablename )
        store[ record.id ] = record
        self._update_store( tablename, store )
        return


    def getRecordByID(self, table, id):
        if not isinstance(table, basestring): table = table.name
        store = self._retrieve_store( table )
        return store[ id ]


    def getRecords(self, table):
        if not isinstance(table, basestring): table = table.name
        store = self._retrieve_store(table)
        return store.iteritems()


    def _retrieve_store(self, tablename):
        path = self._store_path( tablename )
        if not os.path.exists( path ):
            ret = dict()
        else:
            f = open(path)
            ret = pickle.load( f )
            f.close()

        return ret


    def _update_store(self, tablename, newstore):
        store = self._retrieve_store( tablename )
        store.update( newstore )
        path = self._store_path(tablename)
        f = open(path, 'w')
        pickle.dump(store, f)
        f.close()
        return


    def _store_path(self, tablename):
        return os.path.join( self.dbroot, tablename )
    """



# version
__id__ = "$Id$"

# End of file

__date__ = "$Jul 29, 2009 7:18:52 PM$"


