from scipy import *
from scipy.optimize import *
from pylab import *
from dos_utils import *
from matdyn import *

polyOrder = 2
def fitFunc(p, x):
    return p[0] * x + p[1] * x * x 
#+ p[2] * x * x * x + p[3] * x * x * x * x

def fitOmega(prcntVolume, prcntOmega):
#    fitfunc = lambda p, x: p[0] * x + p[1] * x * x
    errFunc = lambda p, x, y: (y - fitFunc(p, x))
    # initial parameters:
    pinit = [prcntOmega[1]/prcntVolume[1], 0.0]
#    pinit = [1,0]
    v, success = leastsq(errFunc, pinit, args=(prcntVolume, prcntOmega) )
    if success < 1 or success > 4:
        print success
        print prcntOmega
        print v
 #   assert success != 1, "fitOmega: Fitting was not successful"
    return v
    
    
def getQuasiFreqs(prcntVol, coeffOmega, Omega0):
    nPoints = shape(coeffOmega)[0]
    nModes = shape(coeffOmega)[1]
    fittedFreqs = zeros(shape = (nPoints, nModes) )
    for i in range( nPoints ):
        for j in range( nModes ):
            p = coeffOmega[i,j]
            fittedFreqs[i,j] = p[0]*prcntVol
    return (fittedFreqs+1.0)*Omega0        
        
def getFittedFreqs(prcntVol, coeffOmega, Omega0):
    nPoints = shape(coeffOmega)[0]
    nModes = shape(coeffOmega)[1]
    fittedFreqs = zeros(shape = (nPoints, nModes) )
    for i in range( nPoints ):
        for j in range( nModes ):
            p = coeffOmega[i,j]
            fittedFreqs[i,j] = fitFunc(p, prcntVol)
    return (fittedFreqs+1.0)*Omega0
    
def loadPhonons(indexRange, fPrefix):
    fname = fPrefix + str(0) + '.modes'
    Pol, Omega, qPoints = matdyn(fname)
    volOmega = zeros(shape=(len(indexRange), shape(Omega)[0], shape(Omega)[1]  ) )
    volPol = zeros(shape=(len(indexRange), shape(Pol)[0], shape(Pol)[1], shape(Pol)[2], shape(Pol)[3] ) )
    volPol[0] = Pol
    volOmega[0] = Omega
    for i in range(1,len(indexRange)):
        fname = fPrefix + str(indexRange[i]) + '.modes'
        Pol, Omega, qPoints = matdyn(fname)
        volPol[i] = Pol
        volOmega[i] = Omega
    return volPol, volOmega, qPoints        
        
def fitMatdyn(indexRange, prefix ):
    volPol, volOmega, qPoints = loadPhonons(indexRange, prefix)
    prcntVol = array(indexRange)/1000.0
    # percent change in Omegas relative to equilibrium    
    volPrcntOmega = zeros(shape(volOmega))
    # introduce Omega0 to get rid of 0s in denominator
    Omega0 = copy(volOmega[0])
    for i in range(shape(Omega0)[0]):
        for j in range(shape(Omega0)[1]):
            if Omega0[i,j] == 0.:
                Omega0[i, j] = 1             
    for i in range(len(indexRange)):
        volPrcntOmega[i] = (volOmega[i] - volOmega[0])/Omega0
    coeffOmega = zeros( shape = ( shape(Omega0)[0], shape(Omega0)[1], polyOrder) )
    nonlinearity = 0
    for i in range( shape(Omega0)[0] ):
        for j in range( shape(Omega0)[1] ):
            coeffOmega[i,j] = fitOmega(prcntVol, volPrcntOmega[:,i,j])
            nonlinearity = nonlinearity + abs(coeffOmega[i,j][1])

#    print volPrcntOmega[:,1000,3]
#    print fitOmega(prcntVol, volPrcntOmega[:,1000,3])
#    plot(volPrcntOmega[:,1000,3])
#    show()

    return prcntVol, coeffOmega, volOmega, nonlinearity
    

if __name__ == '__main__':
    indexRange1 = range(0,5,2)
    print indexRange1
    prcntVol1, coeffOmega1, volOmega1, nonlinearity1 = fitMatdyn(indexRange1, './matdyn')
    indexRange2 = range(0,21,2)
    print indexRange2
    prcntVol2, coeffOmega2, volOmega2, nonlinearity2 = fitMatdyn(indexRange2, './matdyn')
    print 'Change in nonlinearuty (nonlinearity2/nonlinearity1): ', nonlinearity2/nonlinearity1
    lastPrcntVolume1 = prcntVol1[len(indexRange1) - 1]
    lastPrcntVolume2 = prcntVol2[len(indexRange2) - 1]    
    quasiOmega1 = getQuasiFreqs(lastPrcntVolume1, coeffOmega1, volOmega1[0])
    quasiOmega2 = getQuasiFreqs(lastPrcntVolume2, coeffOmega2, volOmega2[0])
    
    deltaOmega = 0.5
    minOmega = 0.0
    maxOmega = 110
    axisOmega = linspace(minOmega,maxOmega,int((maxOmega - minOmega)/deltaOmega))
    quasiDos1 = calcPhonDos(quasiOmega1,minOmega, maxOmega, deltaOmega)
    quasiDos2 = calcPhonDos(quasiOmega2,minOmega, maxOmega, deltaOmega)
    Dos2 = calcPhonDos(volOmega2[len(indexRange2) - 1],minOmega, maxOmega, deltaOmega)
    Dos1 = calcPhonDos(volOmega1[len(indexRange1) - 1],minOmega, maxOmega, deltaOmega)    
    fittedOmega2 = getFittedFreqs(lastPrcntVolume2, coeffOmega2, volOmega2[0])
    fittedDos2 = calcPhonDos(fittedOmega2,minOmega, maxOmega, deltaOmega)

    volPol, volOmega, qPoints = loadPhonons(indexRange2, './matdyn')
    histGamma = calcCoeffPartPhonDos(coeffOmega2, 0, volOmega2[len(indexRange2) - 1], volPol[len(indexRange2) - 1], 0, minOmega, maxOmega, deltaOmega)
    print histGamma
    plot(axisOmega,histGamma, label="Coefficients")
    axis([minOmega, maxOmega, -1, 3])
    show()
#    quasiOmega = getQuasiFreqs(lastPrcntVolume, coeffOmega, volOmega[0])
#    fittedDos = calcPhonDos(fittedOmega,0., 110., 0.5
#    quasiDos = calcPhonDos(quasiOmega,0., 110., 0.5)
#    Dos = calcPhonDos(volOmega[len(indexRange) - 1],0., 110., 0.5)
#    plot(axisOmega,Dos1, label="original DOS")
#    axis([minOmega, maxOmega, 0, 300])
#    plot(axisOmega,quasiDos1, label="Quasi harmonic DOS")    
#    axis([minOmega, maxOmega, 0, 300])
    plot(axisOmega,fittedDos2, label="original DOS")
    axis([minOmega, maxOmega, 0, .04])
    plot(axisOmega,quasiDos2, label="Quasi harmonic DOS")    
    axis([minOmega, maxOmega, 0, .04])
    show()
#    print coeffOmega
#    print volPrcntOmega[:,1000,1]
#    print fitOmega(prcntVol, volPrcntOmega[:,1000,1])
#    plot(volPrcntOmega[:,1000,1])
#    axis([0, 22, 0, -0.04])
#    show()
