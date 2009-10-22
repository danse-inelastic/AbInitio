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

# Can't find SSHer.csaccessor()
def ssher(name='ssher', facility='csaccessor'):
    from SSHer import SSHer
    return SSHer(name=name, facility=facility)

def buildjob(*args, **kwds):
    from JobBuilder import JobBuilder
    builder = JobBuilder()
    return builder(*args, **kwds)

# ?
def retrieveresults(*args, **kwds):
    from ComputationResultsRetriever import ComputationResultsRetriever
    retriever = ComputationResultsRetriever()
    return retriever(*args, **kwds)

# version
__id__ = "$Id$"

# End of file 
