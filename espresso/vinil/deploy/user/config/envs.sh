export EXPORT_ROOT=/home/user/exports/vinil
export EXPORT_SOURCE=/home/user/danse-workspace/AbInitio/espresso/vinil
export PYRE_DIR=/home/user/exports/pythia
export LUBAN_DIR=/home/user/exports/luban

#Need to explicitly specify HOME for .matplotlib (for web app)
# www-data user doesn't set HOME variable
export HOME=/tmp
export ESPRESSO_PATH=/home/user/distribs/espresso-4.0.5/bin
export PATH=$ESPRESSO_PATH:$PATH
export PATH=$EXPORT_ROOT/bin:$PATH

export LD_LIBRARY_PATH=$EXPORT_ROOT/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$EXPORT_ROOT/lib:$DYLD_LIBRARY_PATH

# Setting PYTHONPATH
export PYTHONPATH=$PYRE_DIR/modules:$PYTHONPATH
export PYTHONPATH=$LUBAN_DIR/modules:$PYTHONPATH
export PYTHONPATH=$EXPORT_ROOT/modules:$PYTHONPATH
