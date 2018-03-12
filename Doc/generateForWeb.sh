#!/bin/bash

#
# This script generates the files necessary for the /SixTrack/web/docs folder on the website.
# Requires the latexml package to be installed.
# Written by Veronica Berglyd Olsen, Feb 2018
#

CURR=$(pwd)
MUSER=$CURR/user_manual
MPHYS=$CURR/physics_manual
MBUILD=$CURR/building_sixtrack
OUSERF=$CURR/html/user_full
#OUSERS=$CURR/html/user_split
OPHYS=$CURR/html/physics_full
OBUILD=$CURR/html/build_full

# LaTeXML Options
MATHJAX='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML'
FORMAT=html5

mkdir -pv $OUSERF
#mkdir -pv $OUSERS
mkdir -pv $OPHYS
mkdir -pv $OBUILD

echo ""
echo "*******************************"
echo "* Generating User Manual HTML *"
echo "*******************************"
echo ""

cd $MUSER
make
cp $MUSER/six.pdf $CURR/html/user_manual.pdf
latexml six.tex | latexmlpost --dest=$OUSERF/manual.html --format=$FORMAT --javascript=$MATHJAX -
#latexml six.tex | latexmlpost --dest=$OUSERS/manual.html --format=$FORMAT --javascript=$MATHJAX --splitat=chapter -
$CURR/cleanupHTML.py $OUSERF
rm -v $OUSERF/*.html
echo "<?php header('Location: manual.php'); ?>" > $OUSERF/index.php

echo ""
echo "**********************************"
echo "* Generating Physics Manual HTML *"
echo "**********************************"
echo ""

cd $MPHYS
make
cp $MPHYS/sixphys.pdf $CURR/html/physics_manual.pdf
latexml sixphys.tex | latexmlpost --dest=$OPHYS/manual.html --format=$FORMAT --javascript=$MATHJAX -
$CURR/cleanupHTML.py $OPHYS
rm -v $OPHYS/*.html
echo "<?php header('Location: manual.php'); ?>" > $OPHYS/index.php

echo ""
echo "********************************"
echo "* Generating Build Manual HTML *"
echo "********************************"
echo ""

cd $MBUILD
make
cp $MBUILD/building_sixtrack.pdf $CURR/html/building_sixtrack.pdf
latexml building_sixtrack.tex | latexmlpost --dest=$OBUILD/manual.html --format=$FORMAT --javascript=$MATHJAX -
$CURR/cleanupHTML.py $OBUILD
rm -v $OBUILD/*.html
echo "<?php header('Location: manual.php'); ?>" > $OBUILD/index.php

echo ""
echo "**********"
echo "*  DONE  *"
echo "**********"
echo ""
echo "The content of the folder 'html' can now be uploaded to /afs/cern.ch/project/sixtrack/web/docs/"
echo ""
