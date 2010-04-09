#!/usr/bin/env bash

# Set env var if it is not set
if [ ! -n "$CASSANDRA" ]; then
    source envs.sh
fi

export MASTER_DIR=$CASSANDRA/config/taskmaster
export WORKER_DIR=$CASSANDRA/config/worker

# Start daemons
cassandra start $WORKER_DIR
cassandra start $MASTER_DIR