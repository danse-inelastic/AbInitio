#!/usr/bin/env python
'''
module to parse files like INCAR
and to write them out

the INCAR file is a tagged format free-ASCII file: Each line consists
of a tag (i.e. a string) the equation sign'=' and a number of
values. It is possible to give several parameter-value pairs ( tag =
values ) on a single line, if each of these pairs are separated by a
semicolon ';'. If a line ends with a backslash the next line is a
continuation line. Comments are normally preceded by the number sign
#, but in most cases comments can be append to a parameter-value
pair without the #. In this case semicolons should be avoided within
the comment.

This module does not support the following features:
1. you can not use a backslash to continue a line

This class can be used for any input file like the INCAR file

John Kitchin jkitchin@cmu.edu

Modified to derive from an ordered dictionary class to maintain same order
between read and written files. 
Olivier Delaire - June 2007
'''

import os
from OrderedDict import OrderedDict

class INPUT(OrderedDict):
    '''
    we store all the key=vals in a dictionary to prevent
    having multiple entries.
    '''
    def __init__(self,filename='INCAR'):
        OrderedDict.__init__(self)

        self.filename = filename
        
        if not os.path.exists(filename):
            return

        f = open(filename,'r')
    
        lines = f.readlines()
        f.close()

        for line in lines:
    
            line = line.strip()
            if len(line) == 0:
                #empty line, ignore it
                continue
            
            #strip out comments.  on finding one, i chop off
            #everything to the right of the # sign
            # this mainly happens when reading in other peoples files
            # we never put comments in the input file when auto-generating them.
            if '#' in line:
                i = line.index('#')
                line = line[0:i]
                                
            if '=' not in line:
                #not a keyword=value line
                # we can ignore the line
                continue

            elif ';' in line:
                #need to split up the line
                keyvals = line.split(';')
                for keyval in keyvals:
                    if '=' in keyval:
                        #make sure it is key=val
                        fields = keyval.split('=')
                        if len(fields) != 2:
                            raise Exception, 'key=val split did not give 2 fields:',keyval

                        key = fields[0]
                        val = fields[1]
                        self[key]=val
            else:
                fields = line.split('=')
                if len(fields) != 2:
                    raise Exception, 'key=val split did not give 2 fields:',keyval
        
                key = fields[0]
                val = fields[1]
                self[key]=val

    def __setitem__(self,name,value):
        '''modify the builtin setitem behavior so that everytime
        a parameter is added or changed the INCAR file is written out
        so it is always up to date.
        '''
        OrderedDict.__setitem__(self,name,value)

        if __debug__:
            print 'Updating after setting an item!'
        self.Write()

    def Write(self):
        '''
        write each key=val pair back out.

        the order is random due to the use of the dictionary.
        [ This is improved by using an OrderedDict class,
        that keeps track of the order in which items are
        added to the input dictionary.]
        additional code would be required here to provide any structure.
        '''
        
        f = open(self.filename,'w')
        for key in self:
            f.write('%s = %s\n' % (key,self[key]))
        f.close()
                    


if __name__ == '__main__':

    c = INPUT('tests/CO/INCAR')
    print c

    c2 = INPUT('tests/INCARS/incar2')
    print c2

    c3 = INPUT('tests/INCARS/incar3')
    print c3
    c3.Write('testincar3')
