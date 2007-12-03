f=file('OUTCAR')
#content=f.read()

#import re
##first get the text part we're interested in
#eigPart=re.compile('Eigenvectors and eigenvalues of the dynamical matrix').search(content)#.*Eigenvectors after division by SQRT(mass)',content)
#print eigPart
#now get the values we're intereseted in
#for line in eigPart:
    
    
lines=f.readlines()
for i in range(len(lines)):
    if lines[i]==' Eigenvectors and eigenvalues of the dynamical matrix\n':
        start=i; continue
    if lines[i]==' Eigenvectors after division by SQRT(mass)\n':
        stop=i; break
       
for line in lines[start:stop]:
    print line