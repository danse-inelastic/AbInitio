import numpy as np
import pylab as pl
import py4cs
from py4cs.filetable import read

from elph import *

V3Sifile = open('DOSCAR_V3Si_ecut320_mp30_tetra.dat', 'r')

# first line is a comment:
V3Sifile.readline()

# read in the data
V3Sidata = read(V3Sifile)
V3Sifile.close()


############################
### V3Si
############################


eV3Si = V3Sidata[:,0]
dV3Si = V3Sidata[:,1]
#depV3Si_ref = V3Sidata[:,2]

#############

EfV3Si = 7.20274569  # <<<<<<<<<<<<<<<<< adjust accordingly

############
pl.plot(eV3Si-EfV3Si, dV3Si, 'r-')
pl.show()

pl.xlim((-1,1))

#############

lambda_V3Si = 1.0 

kerV3Si25 = np.array([Lorentzian(i*de, 0.0, 2.0*gamma(lambda_V3Si, 25)) for i in range(len(eV3Si))])
kerV3Si300 = np.array([Lorentzian(i*de, 0.0, 2.0*gamma(lambda_V3Si, 300)) for i in range(len(eV3Si))])
kerV3Si1000 = np.array([Lorentzian(i*de, 0.0, 2.0*gamma(lambda_V3Si, 1000)) for i in range(len(eV3Si))])

minT = 100 # cannot be zero
maxT = 1500
deltaT = 100

Ts = range(minT, maxT, deltaT)
numTs = len(Ts)


kerV3Si_allT = []
for n in range(numTs):
    T = Ts[n]
    kerV3Si = np.array([Lorentzian(i*de, 0.0, 2.0*gamma(lambda_V3Si, T)) for i in range(len(eV3Si))])
    kerV3Si_allT.append(kerV3Si)

dV3Siep_allT = []
for n in range(numTs):
    T = Ts[n]
    dV3Siep = convolveSymArr(eV3Si, dV3Si, kerV3Si_allT[n])
    dV3Siep_allT.append(dV3Siep)

intdV3Si_allT = []
for n in range(numTs):
    T = Ts[n]
    intdV3Si_allT.append(antiderivate(eV3Si, dV3Siep_allT[n]))

muV3Si_allT = []
for n in range(numTs):
    T = Ts[n]
    muV3Si = findFillEnergy(eV3Si, intdV3Si_allT[n], NelV3Si)
    muV3Si_allT.append(muV3Si)

nefV3Si_allT = []
for n in range(numTs):
    T = Ts[n]
    nef = findNEF(eV3Si, dV3Siep_allT[n], muV3Si_allT[n])
    nefV3Si_allT.append(nef)




# convolve the e-dos with the electronic broadening:

dV3Siep25 = convolveSymArr(eV3Si, dV3Si, kerV3Si25)
dV3Siep300 = convolveSymArr(eV3Si, dV3Si, kerV3Si300)
dV3Siep1000 = convolveSymArr(eV3Si, dV3Si, kerV3Si1000)

# calculate the chemical potentials:
# integrate the dos:
intdV3Si = antiderivate(eV3Si, dV3Si)
intdV3Siep25 = antiderivate(eV3Si, dV3Siep25)
intdV3Siep300 = antiderivate(eV3Si, dV3Siep300)
intdV3Siep1000 = antiderivate(eV3Si, dV3Siep1000)


# calc mu

NelV3Si = 38 #check this

# mu for the non-smeared DOS should be close to zero
# but calculate it for consistency
muV3Si = findFillEnergy(eV3Si, intdV3Si, NelV3Si)
muV3Si25 = findFillEnergy(eV3Si, intdV3Siep25, NelV3Si)
muV3Si300 = findFillEnergy(eV3Si, intdV3Siep300, NelV3Si)
muV3Si1000 = findFillEnergy(eV3Si, intdV3Siep1000, NelV3Si)

energyV3Si = np.array(eV3Si)

energyV3Si25 = energyV3Si - muV3Si25
energyV3Si300 = energyV3Si - muV3Si300
energyV3Si1000 = energyV3Si - muV3Si1000

####### Plot results:###########

### define plotting preferences:
fig_width_pt = 600.0  # Sit this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 20,
          'legend.fontsize': 16,
          'text.fontsize': 20,
          'font.size': 20,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'figure.figsize': fig_size}
pl.rcParams.update(params)

###

pl.figure(3)
pl.axes([0.125,0.15,0.95-0.125,0.95-0.15])
pl.ylim((0,3.0))
pl.plot(eV3Si - muV3Si, dV3Si/8, 'k-', label='V$_{3}$Si 0K')
pl.plot(eV3Si - muV3Si, dV3Siep25/8, 'b-', label='V$_{3}$Si 25K')
pl.plot(eV3Si - muV3Si300, dV3Siep300/8, 'g-', label='V$_{3}$Si 300K')
pl.plot(eV3Si - muV3Si1000, dV3Siep1000/8, 'r-', label='V$_{3}$Si 1000K')
pl.legend()
#pl.plot(eV3Si, dV3Siep300, 'g--', label=None)
#pl.plot(eV3Si, dV3Siep1000, 'r--', label=None)
pl.xlabel('$E-\mu(T)$ (eV)')
pl.ylabel('(states/eV/at.)')
pl.xlim((-0.7,0.7))
pl.ylim((0,3.0))
pl.grid(True)
