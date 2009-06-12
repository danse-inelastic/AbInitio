from django.db import models

class Job(models.Model):
    type            = models.CharField(max_length = 50)
    status          = models.CharField(max_length = 50)
    created        = models.DateTimeField() 
    config          = models.TextField()

"""
class Electrons(models.Model):
    ibrav           = models.IntegerField()
    ecutwfc       = models.FloatField()
    occupations = models.CharField(max_length = 100)
    smearing     = models.CharField(max_length = 100)
    degauss       = models.FloatField()
    atomic_species = models.TextField()
    atomic_positions = models.TextField()
    kpoints         = models.TextField()
    configuration = models.TextField()
    job               = models.ForeignKey(Job)
    
class Phonons(models.Model):
    nq1             = models.IntegerField()
    nq2             = models.IntegerField()
    nq3             = models.IntegerField()
    configuration = models.TextField()
    job               = models.ForeignKey(Job)
"""

class Plot(models.Model):
    filename      = models.CharField(max_length = 100)
    size             = models.FloatField() # in KB
    created       = models.DateTimeField()
    job               = models.ForeignKey(Job) # Many plots can be produces during the single Job
    
