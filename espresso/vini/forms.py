from django import forms
#from qe.models import Electrons

class ElectronsForm(forms.Form):
    ibrav           = forms.IntegerField(widget=forms.TextInput(attrs={'size':5}))
    ecutwfc       = forms.FloatField(widget=forms.TextInput(attrs={'size':5}))
    occupations = forms.CharField(max_length = 100,  widget=forms.TextInput(attrs={'size':10}))
    smearing     = forms.CharField(max_length = 100,  widget=forms.TextInput(attrs={'size':10}))
    degauss       = forms.FloatField(widget=forms.TextInput(attrs={'size':5}))
    atomic_species = forms.CharField(widget=forms.TextInput(attrs={'size':20}))
    atomic_positions = forms.CharField(widget=forms.Textarea(attrs={'rows': 2,  'cols': 23}))
    kpoints         = forms.CharField(widget=forms.Textarea(attrs={'rows': 2,  'cols': 23}))
    
class ElectronsConfigForm(forms.Form):
    config = forms.CharField(widget=forms.Textarea(attrs={'rows': 30,  'cols': 50}))
    
class PhononsConfigForm(forms.Form):
    config = forms.CharField(widget=forms.Textarea(attrs={'rows': 15,  'cols': 50}))

