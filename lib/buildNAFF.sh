#!/usr/bin/env bash

set -e #Exit on error

echo ""
echo " Initialising NAFF"
echo "==================="
echo ""

cd ..
git submodule init lib/NAFFlib
git submodule update lib/NAFFlib
cd lib

echo ""
