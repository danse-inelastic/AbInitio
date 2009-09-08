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


# this script is to be used with the SimpleHttpServer.
# Should try to improve the server implementation.

webapp = 'webmain.py'


import os


#The "request" object passed from simple http server
#request = {}
#print request
#convert it to a query string
query_string = '&'.join( '%s=%s' % (k, ','.join(v)) for k,v in request.iteritems() )
#print 'query_string: %s<br><br>' % query_string
os.environ['QUERY_STRING'] = query_string


#headers
#print 'headers: %s<br><br>' % headers
#os.environ[ 'CONTENT_TYPE' ] = headers['content-type']


#posted data
#posted = file_handle_for_posted_data.read()
#print posted


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
    
#from  opaldemo_main import main
#main()


# version
__id__ = "$Id$"

# End of file 
