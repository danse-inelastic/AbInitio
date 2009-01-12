vaspapp.py -name=vasp-FeAl -ecutoff=240.0 -xcf=PAW-PBE -mpmesh=2,2,2 -unitcell=FeAl.xyz -generateInputsOnly=False
phonapp.py -name=phon-FeAl -unitcell=FeAl.xyz -supersize=2,2,2 -qgridsize=20,20,20 -dosaxis=0,20 -amplitude=0.01
