from vasp.parsing.eigenvectorParser import parseModes, parseEsModes

if __name__=='__main__':
    for mode in parseModes():
        print mode
    #print parseEsModes()