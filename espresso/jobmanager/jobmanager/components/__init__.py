# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                      (C) 2006-2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

def dds():
    from DistributedDataStorage import DistributedDataStorage
    return DistributedDataStorage('dds', 'dds')


def clerk():
    from Clerk import Clerk
    return Clerk('clerk', 'clerk')


def usersFromDB():
    from UsersFromDB import UsersFromDB
    return UsersFromDB('usersFromDB')


# version
__id__ = "$Id$"

# End of file 