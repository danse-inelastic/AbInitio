import models
from django.conf.urls.defaults import *

HOME = '/vini/'

urlpatterns = patterns('',
    (r'^$',  'vini.views.welcome'), 
    (r'^static/(?P<path>.*)$', 'django.views.static.serve',
    {'document_root': '/home/dexity/Documents/Work/CalTech/Tasks/espresso/vini'}), # change, to have the 'static' reference work
    (r'^set-electron-params/$',  'vini.views.set_electron_params'),
    (r'^edit-electron-params/$',  'vini.views.edit_electron_params'),
    (r'^saved-electron-params/(\d+)/$',  'vini.views.saved_electron_params'),
    (r'^run-electron-simulation/(\d+)/$',  'vini.views.run_electron_simulation'),
    (r'^electron-jobs/(\d+)/$',  'vini.views.electron_jobs'),
    (r'^electron-dos/$',  'vini.views.electron_dos'),
    (r'^set-phonon-params/$',  'vini.views.set_phonon_params'),
    (r'^edit-phonon-params/$',  'vini.views.edit_phonon_params'),
    (r'^saved-phonon-params/(\d+)/$',  'vini.views.saved_phonon_params'),
    (r'^run-phonon-simulation/(\d+)/$',  'vini.views.run_phonon_simulation'),
    (r'^phonon-jobs/(\d+)/$',  'vini.views.phonon_jobs'),
    (r'^phonon-dos/$',  'vini.views.phonon_dos'),
)

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

# Example:
# (r'^vini/', include('vini.foo.urls')),

# Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
# to INSTALLED_APPS to enable admin documentation:
# (r'^admin/doc/', include('django.contrib.admindocs.urls')),

# Uncomment the next line to enable the admin:
# (r'^admin/(.*)', admin.site.root),
