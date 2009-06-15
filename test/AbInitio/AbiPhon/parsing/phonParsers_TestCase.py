#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                        (C) 2009 All Rights Reserved 
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


import os

import unittest

from unittest import TestCase
class TestCase(TestCase):


    def testDOSmeV(self):
        'AbInitio.AbiPhon.parsing.phonParsers: parseDOS_meV2IDF'
        
        outfile = 'DOS_meV.idf'
        if os.path.exists(outfile):
            os.remove(outfile)
            
        infile = 'DOS.meV'
        self.assert_( os.path.exists(infile) )
        
        from AbInitio.AbiPhon.parsing import phonParsers
        phonParsers.parseDOS_meV2IDF()

        self.assert_(os.path.exists(outfile))
        return
        
    
    pass # end of TestCase


def pysuite():
    suite1 = unittest.makeSuite(TestCase)
    return unittest.TestSuite( (suite1,) )

def main():
    pytests = pysuite()
    alltests = unittest.TestSuite( (pytests, ) )
    unittest.TextTestRunner(verbosity=2).run(alltests)
    return


if __name__ == '__main__': main()
    

# version
__id__ = "$Id$"

# End of file 
