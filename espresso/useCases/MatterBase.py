# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                      California Institute of Technology
#                        (C) 2007  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from vnf.dom.OwnedObject import OwnedObject as base
class MatterBase(base):
    """ stores data for matter object """

    name = 'matterbase'

    import pyre.db

    cartesian_lattice = pyre.db.doubleArray(
        name = 'cartesian_lattice', default = [1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0])
    cartesian_lattice.meta['tip'] = 'array of cartesian lattice vectors in Angstroms'
    
    fractional_coordinates = pyre.db.doubleArray(
        name = 'fractional_coordinates', default = [])
    fractional_coordinates.meta['tip'] = 'array positions as fractional values of unit cell'
    
    atom_symbols = pyre.db.varcharArray(
        name = 'atom_symbols', length = 2, default = [] )
    atom_symbols.meta['tip'] = 'atom symbols for each position in the unit cell'

    chemical_formula = pyre.db.varchar(name = 'chemical_formula', length = 1024)
    
    #add space group
    
    #add bravais lattice type
    
    #see pyre.db.__init__().py for types of variables

# version
__id__ = "$Id$"

# End of file 
