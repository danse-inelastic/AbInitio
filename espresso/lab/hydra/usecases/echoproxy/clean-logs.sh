#!/usr/bin/env bash

# Set env var if it is not set
if [ ! -n "$CASSANDRA" ]; then
    source envs.sh
fi


export LOGFILE="twistd.log"
export MASTER_DIR=$CASSANDRA/config/taskmaster
export WORKER_DIR=$CASSANDRA/config/worker

rm $MASTER_DIR/$LOGFILE
rm $WORKER_DIR/$LOGFILE
