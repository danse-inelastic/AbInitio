#!/usr/bin/env bash

# Tests should be run from tests directory

# Add parser directory to PYTHONPATH first
export PWD=`pwd`
export PYTHONPATH="${PWD}/../parser"

python tests.py