from numpy import *

def calcPhonDos(Omega, minOmega = 0., maxOmega = 110., deltaOmega = 1.0):
    histOmega = zeros(int((maxOmega-minOmega)/deltaOmega))
    norm = 0.0
    for modes in Omega:
        for omega in modes:
            idx = int( (abs(omega) - minOmega)/deltaOmega )
            if idx < len(histOmega):
                histOmega[idx] = histOmega[idx] + 1.0
                norm = norm + 1.0
    return histOmega/norm

def calcPartPhonDos(Omega, PolVec, iAtom, minOmega = 0., maxOmega = 110., deltaOmega = 1.0):
    histPartOmega = zeros(int((maxOmega-minOmega)/deltaOmega))
    norm = 0.0
    for modes, vectors in zip(Omega, PolVec):
        for omega, vector in zip(modes, vectors[:,iAtom,:]):
            idx = int( (abs(omega) - minOmega)/deltaOmega )
            if idx < len(histPartOmega):
                weight = (real(vector*vector.conjugate())).sum()
                histPartOmega[idx] = histPartOmega[idx] + weight
                norm = norm + weight
    return histPartOmega/norm

def  calcCoeffPartPhonDos(coeffOmega, iOrder, Omega, PolVec, iAtom, minOmega = 0., maxOmega = 110., deltaOmega = 1.0):
    histPartOmega = zeros(int((maxOmega-minOmega)/deltaOmega))
    histPartCoeff = zeros(int((maxOmega-minOmega)/deltaOmega))
    norm = 0.0
    for modes, vectors, coeffs in zip(Omega, PolVec, coeffOmega[:,:,iOrder]):
        for omega, vector, coeff in zip(modes, vectors[:,iAtom,:], coeffs):
            idx = int( (abs(omega) - minOmega)/deltaOmega )
            if idx < len(histPartCoeff):
                weight = (real(vector*vector.conjugate())).sum()
                histPartCoeff[idx] = histPartCoeff[idx] - weight*coeff
                histPartOmega[idx] = histPartOmega[idx] + weight
                norm = norm + weight*coeff
    for idx in range(len(histPartCoeff)):
        if abs(histPartOmega[idx]) > 0.0 :
           histPartCoeff[idx] = histPartCoeff[idx]/histPartOmega[idx]
        else:
            histPartCoeff[idx] = 0.0 
    return histPartCoeff
    

if __name__ == '__main__':
     from matdyn import *
     from pylab import *
     Pol, Omega, qPoints = matdyn( 'matdyn0.modes' )
     deltaOmega = 0.5
     minOmega = 0.0
     maxOmega = 110.
     axisOmega = linspace(minOmega,maxOmega,int((maxOmega - minOmega)/deltaOmega))
     histPartOmega0 = calcPartPhonDos(Omega, Pol, 0, minOmega, maxOmega, deltaOmega)/3.
     histPartOmega1 = calcPartPhonDos(Omega, Pol, 1, minOmega, maxOmega, deltaOmega)/3.
     histPartOmega2 = calcPartPhonDos(Omega, Pol, 2, minOmega, maxOmega, deltaOmega)/3.
     histPartOmegaB = histPartOmega1 + histPartOmega2
#     histTotalOmega = histPartOmega0 + histPartOmega1
     totalDos = calcPhonDos(Omega, minOmega, maxOmega, deltaOmega)
     plot(axisOmega,histPartOmega0)
     axis([minOmega, maxOmega, 0, 0.04])
     plot(axisOmega, histPartOmegaB)     
     axis([minOmega, maxOmega, 0, 0.04])
#     plot(histTotalOmega)
#     plot(axisOmega, totalDos)
#     axis([minOmega, maxOmega, 0, 0.15])
     show()
#     print totalDos
#     print histPartOmega0

