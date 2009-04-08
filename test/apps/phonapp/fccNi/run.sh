vaspapp.py -name=vasp-fccNi -ecutoff=240.0 -xcf=PAW-PBE -mpmesh=2,2,2 -unitcell=fccNi.xyz -generateInputsOnly=False
phonapp.py -name=phon-fccNi -unitcell=fccNi.xyz -supersize=4,4,4 -qgridsize=20,20,20 -dosaxis=0,20 -amplitude=0.01
