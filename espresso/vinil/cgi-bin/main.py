#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                 Jiao Lin
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

webapp = 'webmain.py'

import os

# I don't know what this is hack is
#query_string = '&'.join( '%s=%s' % (k, ','.join(v)) for k,v in request.iteritems() )
#os.environ['QUERY_STRING'] = query_string

import tempfile
d = tempfile.mkdtemp()
out = os.path.join(d, 'out.html')
err = os.path.join(d, 'err.html')

cmd = "%s >%s  2>%s" % (webapp, out, err)
if os.system( cmd ):
    print open( err ).read()
else:
    lines = open( out ).readlines()
    print ''.join( lines[1:] )

import shutil
shutil.rmtree(d)

# version
__id__ = "$Id$"

# End of file 
