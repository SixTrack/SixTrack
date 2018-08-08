#!/usr/bin/env bash

rsync -a $1/doc/user_manual .
cd user_manual
make clean
make TEXFLAGS=-interaction=nonstopmode > build.log

EXSTAT=$?

if [ $EXSTAT -ne 0 ]; then
  echo "ERROR Failed to build user manual"
fi
exit $EXSTAT
