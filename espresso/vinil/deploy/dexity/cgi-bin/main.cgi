#!/usr/bin/env bash
EXPORT_ROOT=/home/dexity/exports/vinil
source $EXPORT_ROOT/config/envs.sh
cd $EXPORT_ROOT/cgi-bin && webmain.py $@