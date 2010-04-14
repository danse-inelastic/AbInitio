#!/usr/bin/env bash

# Set env var if it is not set
if [ ! -n "$CASSANDRA" ]; then
    source envs.sh
fi

export WORKER_DIR=$CASSANDRA/config/worker

# Start daemons
cassandra start $WORKER_DIR
