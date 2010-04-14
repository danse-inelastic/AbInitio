#!/usr/bin/env bash

# Set env var if it is not set
if [ ! -n "$CASSANDRA" ]; then
    source envs.sh
fi

export MASTER_DIR=$CASSANDRA/config/taskmaster

# Start daemons
cassandra start $MASTER_DIR
