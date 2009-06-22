from configParser import *
from configPW import *
from MatterBase import MatterBase

qe = qe
def test():
    poly = MatterBase()
    poly.fractional_coordinates=[0,0,0]
    poly.cartesian_lattice=[]
    
    #use svn+ssh://svn@danse.us/diffraction/diffraction/diffpy.Structure/trunk
    # if you feed it fractional coordinates and lattice params, it will generate 
    #ibrav=2,    # 
    #celldm(1) =6.65,
    #nat=  1,
    #ntyp= 1,
    
    addNamelistParam('control', 'calculation', "'scf'")
    addNamelistParam('system', 'ibrav', 2)
    save("testni.in")

if __name__ == "__main__":
    test()


