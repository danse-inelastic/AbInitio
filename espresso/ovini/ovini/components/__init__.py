# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

#from vnf import extensions

def clerk():
    from Clerk import Clerk
    """
    from Clerk import Clerk as base, findClerks
    Clerks = findClerks(extensions)

    from _extend_class import subclassOf
    Clerk = subclassOf([base]+Clerks)

    from Clerk import DeepCopier as base, findDeepCopiers
    DeepCopiers = findDeepCopiers(extensions)
    DeepCopier = subclassOf([base]+DeepCopiers)

    Clerk.DeepCopier = DeepCopier
    """

    return Clerk( 'clerk', 'clerk' )

__date__ = "$Jul 19, 2009 9:52:39 PM$"


