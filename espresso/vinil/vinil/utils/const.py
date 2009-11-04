#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Available packages
PACKAGES    = ("Quantum Espresso",)  #, "VASP", "GULP"]  # Packages

# Type of configuration files
TYPES       = ("PW", "PH", "PP")  # "BANDS", "CPPP", "D3", "DOS", "DYNMAT", "INITIAL_STATE", "GIPAW", "D1", "MATDYN", "PROJWFC", "PWCOND", "Q2R" 

# Steps of job creation
STEPS       = ("Create Simulation",
               "Create Configuration",
               "Set Simulation Parameters",
               "Add to Queue")

# Types of simulations
SIMULATIONS = ( "Total Energy",
                "Electron DOS",
                "Electron Dispersion",
                "Geometry Optimization",
                "Single-Phonon",
                "Multi-Phonon DOS",
                "Multi-Phonon Dispersion" )

SERVERS     = ("foxtrot.danse.us",
               #"octopod.danse.us",
               #"upgrayedd.danse.us",
               #"teragrid"
              )

__date__ = "$Nov 3, 2009 3:12:34 PM$"


