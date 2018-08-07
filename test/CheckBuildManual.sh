#!/usr/bin/env bash

cd ../doc/user_manual
make clean
make TEXFLAGS=-interaction=nonstopmode

EXSTAT=$?

if [ $EXSTAT -ne 0 ]; then
  echo "ERROR Failed to build user manual"
fi
exit $EXSTAT
