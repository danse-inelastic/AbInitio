vaspapp.py -name=vasp-Si -ecutoff=240.0 -xcf=PAW-PBE -mpmesh=2,2,2 -unitcell=Si.xyz -generateInputsOnly=False
phonapp.py -name=phon-Si -unitcell=Si.xyz -supersize=2,2,2 -qgridsize=20,20,20 -dosaxis=0,20 -amplitude=0.01
