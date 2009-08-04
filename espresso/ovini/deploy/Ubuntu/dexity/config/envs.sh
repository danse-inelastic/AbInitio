EXPORT_ROOT=/home/dexity/exports/ovini
EXPORT_SOURCE=/home/dexity/danse-workspace/AbInitio/espresso/ovini

#Need to explicitly specify HOME for .matplotlib (for web app)
# www-data user doesn't set HOME variable
export HOME=/tmp
export ESPRESSO_PATH=/home/dexity/distribs/espresso-4.0.5/bin
export PATH=$ESPRESSO_PATH:$PATH

export PYRE_DIR=/home/dexity/exports/pythia
export PATH=$EXPORT_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$EXPORT_ROOT/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$EXPORT_ROOT/lib:$DYLD_LIBRARY_PATH
export PYTHONPATH=$PYRE_DIR/modules:$PYTHONPATH
export PYTHONPATH=$EXPORT_ROOT/modules:$PYTHONPATH
