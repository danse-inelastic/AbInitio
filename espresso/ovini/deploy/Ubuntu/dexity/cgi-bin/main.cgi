#!/usr/bin/env bash
EXPORT_ROOT=/home/dexity/exports/ovini
source $EXPORT_ROOT/config/envs.sh
cd $EXPORT_ROOT/cgi && python main.py $@