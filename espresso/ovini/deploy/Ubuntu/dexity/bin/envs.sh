ROOT=/home/dexity/exports/ovini
EXPORT_SOURCE=/home/dexity/danse-workspace/AbInitio/espresso/ovini

export PYRE_DIR=/home/dexity/exports/pythia
export PATH=$ROOT/bin:$PATH
export LD_LIBRARY_PATH=$ROOT/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$ROOT/lib:$DYLD_LIBRARY_PATH
export PYTHONPATH=$PYRE_DIR/modules:$PYTHONPATH

export PYTHONPATH=$EXPORT_SOURCE:$PYTHONPATH
