#!/usr/bin/env bash
EXPORT_ROOT=/home/dexity/exports/ovini
source $EXPORT_ROOT/bin/envs.sh
cd $EXPORT_ROOT/cgi && python main.py $@