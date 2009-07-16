from qe.models import *
import utils
import datetime
from django.http import HttpResponse,  HttpResponseRedirect
from django.template.loader import get_template
from django.shortcuts import render_to_response
from forms import *

HOME = '/vini/'

electron_config = """
 &control
    calculation='scf'
    restart_mode='from_scratch',
    tprnfor = .true.
    prefix='ni',
    pseudo_dir = './',
    outdir='temp/'
 /
 &system    
    ibrav=2, 
    celldm(1) =6.65, 
    nat=  1, 
    ntyp= 1,
    nspin=2,
    starting_magnetization(1)=0.5,
    degauss=0.02,
    smearing='gauss',
    occupations='smearing',
    ecutwfc =27.0
    ecutrho =300.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00 
K_POINTS AUTOMATIC
4 4 4 1 1 1
"""

phonon_config = """
phonons of Ni at gamma
 &inputph
  tr2_ph=1.0d-16,
  prefix='ni',
  ldisp=.true.,
  nq1=2,
  nq2=2,
  nq3=2,
  amass(1)=58.6934,
  outdir='temp/',
  fildyn='ni.dyn',
 /
"""

def welcome(request):
    return render_to_response('welcome.html')

def set_electron_params(request):
    electron = {'ibrav':  2,  
                    'ecutwfc': 27.0, 
                    'occupations': 'smearing', 
                    'smearing': 'gauss', 
                    'degauss': 0.02, 
                    'atomic_species': 'Ni  26.98', 
                    'atomic_positions': 'Ni 0.00 0.00 0.00', 
                    'kpoints': '4 4 4 1 1 1'
                    }
    
    if request.method == 'POST':
        form = ElectronsForm(request.POST)
        if form.is_valid():
            ibrav           = form.cleaned_data['ibrav']
            ecutwfc      = form.cleaned_data['ecutwfc']
            occupations = form.cleaned_data['occupations']
            smearing    = form.cleaned_data['smearing']
            degauss     = form.cleaned_data['degauss']
            atomic_species = form.cleaned_data['atomic_species']
            atomic_positions = form.cleaned_data['atomic_positions']
            kpoints       = form.cleaned_data['kpoints']

            return HttpResponseRedirect('/edit-electron-params/')
    else:
        form = ElectronsForm(initial=electron)
        
    return render_to_response('set_electron_params.html',  {'form': form})

def edit_electron_params(request):
    if request.method == 'POST':
        form = ElectronsConfigForm(request.POST)
        if form.is_valid():
            config = form.cleaned_data['config']
            job = Job(type='electron',  status='configured',  created=datetime.datetime.now(),  config=config)
            job.save()
            #job_id = 1
            return HttpResponseRedirect('/saved-electron-params/%s/' % job.id)#job_id) 
    else:
        form = ElectronsConfigForm(initial= {'config': electron_config})

    return render_to_response('edit_electron_params.html',  {'form': form})

def saved_electron_params(request,  job_id):
    job = Job.objects.get(id=job_id)
    if request.method == 'POST':
        return HttpResponseRedirect('/run-electron-simulation/%s/' % job.id) 
        
    config = job.config #.replace('\n',  '<br />')
    return render_to_response('saved_electron_params.html',  {'config': config})
    

def run_electron_simulation(request,  job_id):
    # Run simulation, update, job's status, then redirect to Jobs
    job = Job.objects.get(id=job_id)
    job.status = "running"
    job.save()
    utils.run_pw_simulation(infile   = "ni.scf.in",  outfile = "ni.scf.out")
    utils.run_pw_dos(infile = "ni.scf.dos.in")
    utils.create_pw_plot(infile = "ni.scf.dos.out",  imagefile = "ni_scf_dos")
    job.status = "finished"
    job.save()

    return HttpResponseRedirect('/electron-jobs/%s/' % job_id) 

def electron_jobs(request,  job_id):
    job = Job.objects.get(id=job_id)
    
    return render_to_response('electron_jobs.html',  {'job': job})
    
def electron_dos(request):
    return render_to_response('electron_dos.html')

def set_phonon_params(request):
    if request.method == 'POST':
        return HttpResponseRedirect('/edit-phonon-params/')
        
    return render_to_response('set_phonon_params.html')

def edit_phonon_params(request):
    if request.method == 'POST':
        form = PhononsConfigForm(request.POST)
        if form.is_valid():
            config = form.cleaned_data['config']
            job = Job(type='phonon',  status='configured',  created=datetime.datetime.now(),  config=config)
            job.save()
            #job_id = 2
            return HttpResponseRedirect('/saved-phonon-params/%s/' % job.id) #job_id) 
    else:
        form = PhononsConfigForm(initial= {'config': phonon_config})

    return render_to_response('edit_phonon_params.html',  {'form': form})

def saved_phonon_params(request,  job_id):
    job = Job.objects.get(id=job_id)
    if request.method == 'POST':
        return HttpResponseRedirect('/run-phonon-simulation/%s/' % job_id) 
        
    config = job.config #.replace('\n',  '<br />')
    return render_to_response('saved_phonon_params.html',  {'config': config}) 

def run_phonon_simulation(request,  job_id):
    # Run simulation, update, job's status, then redirect to Jobs
    job = Job.objects.get(id=job_id)
    job.status = "running"
    job.save()
    #utils.run_pw_simulation(infile   = "ni.scf.in",  outfile = "ni.scf.out")
    #utils.run_pw_dos(infile = "ni.scf.dos.in")
    utils.create_ph_plot(infile = "ni.ph.dos.out",  imagefile = "ni_ph_dos")
    job.status = "finished"
    job.save()

    return HttpResponseRedirect('/phonon-jobs/%s/' % job_id)

def phonon_jobs(request,  job_id):
    job = Job.objects.get(id=job_id)
    return render_to_response('phonon_jobs.html',  {'job': job})
    
def phonon_dos(request):
    return render_to_response('phonon_dos.html')


