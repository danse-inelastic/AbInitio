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


"""
def findClerks(extensions):
    s = 'from vnf.components.%s.Clerk import Clerk'
    def _(ext):
        exec s % ext in locals()
        return Clerk
    return [ _(ext) for ext in extensions ]


def findDeepCopiers(extensions):
    s = 'from vnf.components.%s.Clerk import DeepCopier'
    def _(ext):
        exec s % ext in locals()
        return DeepCopier
    return [ _(ext) for ext in extensions ]

from copy import copy
from vnf import Undef
from vnf.utils.variables import Variable, LazyValue
from vnf.utils.expr import (Expr, Select, Insert, Update, Delete, Column, Count, Max, Min,
    Avg, Sum, Eq, And, Asc, Desc, compile_python, compare_columns, SQLRaw,
    Union, Except, Intersect, Alias, SetExpr)
from vnf.utils.exceptions import (
    WrongStoreError, NotFlushedError, OrderLoopError, UnorderedError,
    NotOneError, FeatureError, CompileError, LostObjectError, ClassInfoError)
from vnf.utils.info import get_cls_info, get_obj_info
"""

from pyre.components.Component import Component

class Clerk(Component):

    class Inventory(Component.Inventory):

        import pyre.inventory

        # properties
        db = pyre.inventory.str(name='db', default='ovini')
        db.meta['tip'] = "the name of the database"

        dbwrapper = pyre.inventory.str(name='dbwrapper', default='psycopg')
        dbwrapper.meta['tip'] = "the python package that provides access to the database back end"

    def __init__(self, *args, **kwds):
        Component.__init__(self, *args, **kwds)

    def _init(self):
        Component._init(self)

        # connect to the database
        import pyre.db
        dbkwds = DbAddressResolver().resolve(self.inventory.db)
        self.db = pyre.db.connect(wrapper=self.inventory.dbwrapper, **dbkwds)

    def _configure(self):
        Component._configure(self)
        self.db = self.inventory.db

    def getJob(self, id):
        '''retrieve job record specified by id'''
        from ovini.dom.Job import Job
        return self._getRecordByID( Job, id )

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

    """
    def getRecordByID(self, tablename, id):
        from pyre.db.Table import Table as TableBase
        if isinstance(tablename, basestring):
            Table = self._getTable(tablename)
        elif issubclass(tablename, TableBase):
            Table = tablename
        else:
            raise ValueError, 'tablename must be a string or a table class: %s' % tablename
        return self._getRecordByID(Table, id)

    def _getTable(self, tablename):
        from vnf.dom.registry import tableRegistry as registry
        try: return registry.get(tablename)
        except KeyError:
            # backward compatibility
            candidate = tablename.lower() + 's'
            return registry.get(candidate)
    """

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

__date__ = "$Jul 29, 2009 7:18:52 PM$"


