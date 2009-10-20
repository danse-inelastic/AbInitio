from volutherm import VoluTherm
import numpy
import os

# AlB2:
a_range = numpy.array([5.67151771874, 5.6746747706760061, 5.6776204820139551, 5.6806360863588239, 5.683679971382686, 5.686754414520701, 5.6898242129665277, 5.695890566218945, 5.6989053376288545, 5.7018340788532464, 5.7049184321499897, 5.7078757068584496, 5.7107870032685142])
c_a_range = numpy.array([1.09052436402, 1.0908826780484495, 1.0913596389405691, 1.0917930570942032, 1.0922069499496534, 1.0926000248067018, 1.0929925903638811, 1.0938405493690089, 1.0942595677867497, 1.0947249745271406, 1.0950976015945779, 1.0955402735584503, 1.0960062816752161])

e_total = numpy.array([-18.01301094, -18.01300702, -18.01299464, -18.01297406, -18.01294522, -18.01290826, -18.01286360, -18.01274950, -18.01268107, -18.01260483, -18.01252076, -18.01242918, -18.01233015])



indexRange = range(0,9,2)
prcntVol = 2.0*numpy.array(indexRange)/1000.0

c_range = c_a_range*a_range
voluTherm = VoluTherm('config.ini', indexRange, prcntVol, a_range[indexRange], \
                      c_range[indexRange], e_total[indexRange])

#voluTherm.totalFreeEnergy(0.0,300.0)
#voluTherm.minimizeFreeEnergy(temperature = 300.0)
line = ''
for v,i in zip(prcntVol, indexRange):
    Cv, total, electronic, phonon = voluTherm.totalFreeEnergy(v,300.0)
    #v, total, electronic, phonon = (0,0,0,0)
    line = line + '%f    %f    %f    %f    %f\n'%(v, total, electronic, phonon, Cv)
    print line
        # save dos:
    for atom in voluTherm.voluPhon.phonQHA.structure.diffpy():
        #for atom in mp.structure.diffpy()
        oldfname = 'projected_DOS.'+ atom.element
        fname = str(i) + '_' + oldfname
        os.system('cp ' + oldfname + ' ' + fname)
file = open('free_energy.out', 'w')
file.write(line)
file.close()

#for pV in prcntVol:
#    print voluTherm.voluPhon.a.fittedValue(pV)
#    print voluTherm.voluPhon.c.fittedValue(pV)
#    print voluTherm.voluPhon.energy.fittedValue(pV)
