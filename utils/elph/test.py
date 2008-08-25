import numpy as np
import pylab as pl

# this file uses the py4cs 'read' to read in column-wise input data
# replace with custom 'read' function if not available
import py4cs
from py4cs.filetable import read

from elph import *

VVfile = open('VV_edos_vasp_tetra.dat', 'r')
VCofile = open('V15Co_edos_vasp_tetra.dat', 'r')

VVfile.readline()
#VVdata=eval(VVfile.read())
VVdata = read(VVfile)
VVfile.close()

VCofile.readline()
VCodata = read(VCofile)
VCofile.close()

eVV = VVdata[:,0]
dVV = VVdata[:,1]
depVV_ref = VVdata[:,2]

eVCo = VCodata[:,0]
dVCo = VCodata[:,1]
depVCo_ref = VCodata[:,2]

#############
pl.plot(eVV, dVV, 'b-', eVCo, dVCo, 'r-')
pl.show()

pl.xlim((-2,2))

#############

de = eVV[1] - eVV[0]
# de = 0.005

lambda_VV = 0.7
lambda_VCo = 0.5

kerVV300 = np.array([ep.Lorentzian(i*de, 0.0, 2.0*ep.gamma(lambda_VV, 300)) for i in range(len(e))])
kerVV1000 = np.array([ep.Lorentzian(i*de, 0.0, 2.0*ep.gamma(lambda_VV, 1000)) for i in range(len(e))])
kerVCo300 = np.array([ep.Lorentzian(i*de, 0.0, 2.0*ep.gamma(lambda_VCo, 300)) for i in range(len(e))])
kerVCo1000 = np.array([ep.Lorentzian(i*de, 0.0, 2.0*ep.gamma(lambda_VCo, 1000)) for i in range(len(e))])

# convolve the e-dos with the electronic broadening:
dVVep300 = convolveSymArr(eVV, dVV, kerVV300)
dVVep1000 = convolveSymArr(eVV, dVV, kerVV1000)

dVCoep300 = convolveSymArr(eVCo, dVCo, kerVCo300)
dVCoep1000 = convolveSymArr(eVCo, dVCo, kerVCo1000)

# calculate the chemical potentials:
# integrate the dos:
intdVV = antiderivate(eVV, dVV)
intdVVep300 = antiderivate(eVV, dVVep300)
intdVVep1000 = antiderivate(eVV, dVVep1000)

intdVCo = antiderivate(eVCo, dVCo)
intdVCoep300 = antiderivate(eVCo, dVCoep300)
intdVCoep1000 = antiderivate(eVCo, dVCoep1000)

# calc mu
NelVV = 5.0
NelVCo = 5.25

# mu for the non-smeared DOS should be close to zero
# but calculate it for consistency
muVV = findmu(eVV, intdVV, NelVV)
muVV300 = findmu(eVV, intdVVep300, NelVV)
muVV1000 = findmu(eVV, intdVVep1000, NelVV)

muVCo = findmu(eVCo, intdVCo, NelVCo)
muVCo300 = findmu(eVCo, intdVCoep300, NelVCo)
muVCo1000 = findmu(eVCo, intdVCoep1000, NelVCo)

energyVV = np.array(eVV)
energyVCo = np.array(eVCo)

energyVV300 = energyVV - muVV300
energyVV1000 = energyVV - muVV1000

energyVCo300 = energyVCo - muVCo300
energyVCo1000 = energyVCo - muVCo1000

# plot results:

fig_width_pt = 600.0  # Get this from LaTeX using \showthe\columnwidth
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


pl.figure(1)
pl.axes([0.125,0.15,0.95-0.125,0.95-0.15])
pl.plot(eVV - muVV, dVV, 'b-', label='V 0K')
pl.plot(eVV - muVV300, dVVep300, 'g-', label='V 300K')
pl.plot(eVV - muVV1000, dVVep1000, 'r-', label='V 1000K')
pl.legend()
pl.plot(eVV, dVVep300, 'g--', label=None)
pl.plot(eVV, dVVep1000, 'r--', label=None)
pl.xlabel('$E-\mu(T)$ (eV)')
pl.ylabel('(states/eV/at.)')
pl.xlim((-0.7,0.7))
pl.ylim((0,3.0))
pl.grid(True)


pl.figure(2)
pl.axes([0.125,0.15,0.95-0.125,0.95-0.15])
pl.ylim((0,3.0))
pl.plot(eVCo - muVCo, dVCo, 'b-', label='V$_{15}$Co$_1$ 0K')
pl.plot(eVCo - muVCo300, dVCoep300, 'g-', label='V$_{15}$Co$_1$ 300K')
pl.plot(eVCo - muVCo1000, dVCoep1000, 'r-', label='V$_{15}$Co$_1$ 1000K')
pl.legend()
pl.plot(eVCo, dVCoep300, 'g--', label=None)
pl.plot(eVCo, dVCoep1000, 'r--', label=None)
pl.xlabel('$E-\mu(T)$ (eV)')
pl.ylabel('(states/eV/at.)')
pl.xlim((-0.7,0.7))
pl.ylim((0,3.0))
pl.grid(True)

# calculate the change in density at the Fermi level:
nefVV = findNEF(eVV, dVV, muVV)
nefVV1000 = findNEF(eVV, dVVep1000, muVV1000)

nefVCo = findNEF(eVCo, dVCo, muVCo)
nefVCo1000 = findNEF(eVCo, dVCoep1000, muVCo1000)

# entropy stuff:
Ts = np.arange(0, 4200, 200)
Ts[0] = 10

Sbare = np.array([Selbare(eVV,dVV, T) for T in Ts])
Sep = np.array([Selph(eVV, dVV, 5.0, 0.7, T) for T in Ts])

pl.plot(Ts, Sbare, 'b+-', Ts, Sep, 'r+-')

Cbare = Ts[0:-1] * derivate(Ts, Sbare)
Cep = Ts[0:-1] * derivate(Ts, Sep)

pl.plot(Ts[0:-1], Cbare, 'b+-', Ts[0:-1], Cep, 'r+-')

