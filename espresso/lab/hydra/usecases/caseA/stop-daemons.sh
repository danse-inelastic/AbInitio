#!/usr/bin/env bash

if [ -f twistd.pid ]; then
    kill `cat twistd.pid`
else
    echo "File twistd.pid does not exist!"
fi