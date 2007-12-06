#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Patrick Hung, Caltech
#                        (C) 1998-2006  All Rights Reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
fortraninterface: creates fortran name translations

processes the file fortrannames and create the file phoninterface.h
"""

def get_func_list(fil):
    return [l.strip() for l in open(fil,'r').readlines()]

def F77EXTERNS_LOWERCASE_TRAILINGBAR(str):
    cfunc = "%s" % str.lower()
    ffunc = "%s_" % str.lower()
    return "#define %s %s" % (cfunc, ffunc)

def F77EXTERNS_NOTRAILINGBAR(str):
    cfunc = "%s" % str.lower()
    ffunc = "%s" % str
    return "#define %s %s" % (cfunc, ffunc)

def F77EXTERNS_EXTRATRAILINGBAR(str):
    cfunc = "%s" % str.lower()
    ffunc = "%s__" % str
    return "#define %s %s" % (cfunc, ffunc)

def F77EXTERNS_UPPERCASE_NOTRAILINGBAR(str):
    cfunc = "%s" % str.lower()
    ffunc = "%s" % str.upper()
    return "#define %s %s" % (cfunc, ffunc)

def F77EXTERNS_COMPAQ_F90(str):
    # symbols that contain underbars get two underbars at the end
    # symbols that do not contain underbars get one underbar at the end
    cfunc = "%s" % str.lower()
    if '_' in str:
        ffunc = "%s__" % str
    else:
        ffunc = "%s_" % str
    return "#define %s %s" % (cfunc, ffunc)


if __name__=='__main__':
    import datetime
    file = open('phoninterface.h','w')

    file.write("#ifndef __PHONINTERFACE_H__\n")
    file.write("#define __PHONINTERFACE_H__\n")

    file.write("\n/* Automatically generated on %s */\n" % datetime.date.today() )
    file.write("/* Do not include this file directly.*/\n");
    file.write("/* It should be included via going through phonexterns.h */\n\n");

    functions = get_func_list('fortrannames')
    translators =  ["%s" % f for f in dir() if f[0:3] == "F77"]
    for i, T in enumerate(translators):

        if i == 0:
            file.write( "#if defined(%s)\n\n" % T )
        else:
            file.write( "#elif defined(%s)\n\n" % T )

        for f in functions:
            file.write( eval('%s("%s")' %(T,f)) + "\n") 

        file.write('\n\n')
 
    file.write("#endif\n")

    file.write("#endif\n")

    file.close()

# version
__id__ = "$Id: dynainterface.py,v 1.2 2006/09/06 23:09:59 cummings Exp $"

# End of file
