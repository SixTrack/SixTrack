#!/bin/bash

#
# This script generates the files necessary for the /SixTrack/web/docs folder on the website.
# Requires the latexml package to be installed as well as image magick for figure conversion.
# Written by Veronica Berglyd Olsen, Feb 2018. Updated 2019-04-08
# Note: Does not seem to work with latexml 8.3. Please use 8.2.
#

CURR=$(pwd)
MUSER=$CURR/user_manual
TUSER=$CURR/user_manual_temp
MPHYS=$CURR/physics_manual
MBUILD=$CURR/building_sixtrack
OUSERF=$CURR/html/user_full
OPHYS=$CURR/html/physics_full
OBUILD=$CURR/html/build_full

# LaTeXML Options
MATHJAX='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML'
FORMAT=html5

mkdir -pv $OUSERF
mkdir -pv $OPHYS
mkdir -pv $OBUILD

rm -r html
mkdir html

echo ""
echo "*****************************************"
echo "* Generating User Manual LaTeX and HTML *"
echo "*****************************************"
echo ""

# Make temp directory
rm -rfv $TUSER
rsync -avPh $MUSER/ $TUSER
cd $TUSER

# Build PDF
make clean
make | tee ../latexUserManual.log
cp $TUSER/six.pdf $CURR/html/user_manual.pdf

# Remove unsupported stuff for LaTeX source for HTML build
for FILE in *.tex; do
  echo "Cleaning up file $FILE"
  sed -i 's/\\begin{cverbatim}/\\begin{verbatim}/g' $FILE
  sed -i 's/\\end{cverbatim}/\\end{verbatim}/g' $FILE
  sed -i 's/\\begin{ctverbatim}/\\begin{verbatim}/g' $FILE
  sed -i 's/\\end{ctverbatim}/\\end{verbatim}/g' $FILE
  sed -i 's/\\begin{longtabu}/\\begin{tabular}/g' $FILE
  sed -i 's/\\end{longtabu}/\\end{tabular}/g' $FILE
  sed -i 's/\\arraybackslash//g' $FILE
  sed -i '/\\todo/d' $FILE
  sed -i '/\\pdfbookmark/d' $FILE
  sed -i '/\\tabulinesep/d' $FILE
done

# Build HTML
latexml six.tex | latexmlpost --dest=$OUSERF/manual.html --format=$FORMAT --javascript=$MATHJAX - | tee ../htmlUserManual.log
HTMEX=$?
$CURR/cleanupHTML.py $OUSERF
rm -v $OUSERF/*.html
echo "<?php header('Location: manual.php'); ?>" > $OUSERF/index.php
rm -rfv $TUSER

echo ""
echo "***********************************"
echo "* Generating Physics Manual LaTeX *"
echo "***********************************"
echo ""

cd $MPHYS
make | tee ../latexPhysicsManual.log
cp $MPHYS/sixphys.pdf $CURR/html/physics_manual.pdf

echo ""
echo "*********************************"
echo "* Generating Build Manual LaTeX *"
echo "*********************************"
echo ""

cd $MBUILD
make | tee ../latexBuildManual.log
cp $MBUILD/building_sixtrack.pdf $CURR/html/building_sixtrack.pdf

echo ""
echo "**********"
echo "*  DONE  *"
echo "**********"

if [ $HTMEX != 0 ]; then
  echo ""
  echo "ERROR during HTML conversion. It may have failed."
  echo "If the error is in converting figures to png, make sure the line"
  echo "  <policy domain=\"coder\" rights=\"read|write\" pattern=\"PDF\" />"
  echo "is set to read|write in /etc/ImageMagick/policy.xml"
fi

echo ""
echo "The content of the folder 'html' can now be uploaded to /afs/cern.ch/project/sixtrack/web/docs/"
echo "RUN: rsync -rvPh html/ /afs/cern.ch/project/sixtrack/docs/"
echo ""
