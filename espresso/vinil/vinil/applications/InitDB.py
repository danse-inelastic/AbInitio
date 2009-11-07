#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Initialize vinil database tables. 
This will remove all existing tables, so be careful!
"""

from pyre.applications.Script import Script

class InitDB(Script):
    class Inventory(Script.Inventory):

        import pyre.inventory

        import vinil.components
        clerk = pyre.inventory.facility(name="clerk", factory=vinil.components.clerk)
        clerk.meta['tip'] = "the component that retrieves data from the various database tables"

        wwwuser     = pyre.inventory.str(name='www-data', default='')
        tables      = pyre.inventory.list(name='tables', default=[])


    def main(self, *args, **kwds):
        self.db.autocommit(True)

        tables = self.tables
        if not tables:
            from vinil.dom import tables as alltables
            tables = alltables()    # list of database table objects

        for t in tables[::-1]:
            self.dropTable( t )     # Reversed order (needed to handle references)

        for t in tables:
            self.createTable( t )
            if self.wwwuser: self.enableWWWUser( t )

        for t in tables:
            self.initTable( t )

        return


    def createTable(self, table):
        # create the component table
        print " -- creating table %r" % table.name
        try:
            self.db.createTable(table)
        except self.db.ProgrammingError, msg:
            print "    failed; table exists?"
            print msg
        else:
            print "    success"

        return


    def dropTable(self, table):
        print " -- dropping table %r" % table.name
        try:
            self.db.dropTable(table)
        except self.db.ProgrammingError:
            print "    failed; table doesn't exist?"
        else:
            print "    success"

        return


    def initTable(self, table):
        module = table.__module__
        m = __import__( module, {}, {}, [''] )
        inittable = m.__dict__.get( 'inittable' )
        if inittable is None: return
        print " -- Inialize table %r" % table.name
        try:
            inittable(self.clerk)    #self.db )
        except self.db.IntegrityError:
            print "    failed; records already exist?"
        else:
            print "    success"

        return


    def enableWWWUser(self, table):
        print " -- Enable www user %r for table %r" % (self.wwwuser, table.name)
        sql = 'grant all on table "%s" to "%s"' % (table.name, self.wwwuser)
        c = self.db.cursor()
        c.execute(sql)
        return


    def __init__(self):
        Script.__init__(self, 'initdb')
        self.db = None
        return


    def _configure(self):
        Script._configure(self)
        self.clerk          = self.inventory.clerk
        self.clerk.director = self
        self.wwwuser        = self.inventory.wwwuser
        self.tables         = self.inventory.tables
        return


    def _init(self):
        Script._init(self)

        self.db = self.clerk.db
        return


    def _getPrivateDepositoryLocations(self):
        return ['../config']

    
__date__ = "$Nov 6, 2009 4:06:51 PM$"

# **************** DEAD CODE *****************
#        self.idd = self.inventory.idd


#        import pyre.idd
#        idd = pyre.inventory.facility('idd-session', factory=pyre.idd.session, args=['idd-session'])
#        idd.meta['tip'] = "access to the token server"

#        else:
#            tables = [self.clerk._getTable(t) for t in tables]

        # initialize table registry
#        import vnf.dom
#        vnf.dom.register_alltables()

        # id generator
#        def guid(): return '%s' % self.idd.token().locator
#        import vnf.dom
#        vnf.dom.set_idgenerator( guid )

