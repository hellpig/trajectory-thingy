#!/bin/bash
# Creates files: code, data.dat, data2.dat
# Makes plots (final plot has a "pause" statement)
# All units are MKS.
# Run via...
#   ./run.sh
# If permissions are a problem, run...
#   chmod 777 run.sh
# If dependencies are a problem, run...
#   sudo apt-get install gnuplot build-essential

# To run on macOS...
# To allow 3D graph to be rotated,
#   I needed to do "set term x11 enhanced", and I need XQuartz installed.
# To get the dependencies, I use MacPorts...
#   sudo port install gnuplot gcc9
# Note: on default macOS, g++ is a link for clang, which works,
#   but I replace g++ with g++-mp-9 after installing gcc9 via MacPorts.

# To run on Windows using Cygwin64...
#  - install:  https://cygwin.com/install.html
#       making sure to select gcc-g++ and bash
#  - install gnuplot (cygwin's didn't work):  https://sourceforge.net/projects/gnuplot/
#  - add the following bin folders to Path in Windows (your folders may differ)...
#       C:\Program Files\gnuplot\bin
#       C:\cygwin64\bin
#  - in this run.sh, add .exe to executables (maybe not necessary): g++, code, and gnuplot
#  - in this code and plot.gpl, remove "set term x11 enhanced" stuff (or change to "set term windows")
#  - in this run.sh, replace   echo \"0 0\"   with   ECHO 0 0
#  - in PowerShell (or Command Prompt), go to directory then type:  bash run.sh




g++ -O3 -fwhole-program -march=native -fno-stack-protector code.cpp -o code

# run the code to create the data files
./code


gnuplot -p -e "set term x11 enhanced; plot 'data2.dat' u 4:5 title 'deviations (y is North)', '<echo \"0 0\"' u 1:2 title 'no deviation'"

gnuplot -p -e "set term x11 enhanced; plot 'data.dat' u 9:8 with lines title 'altitude(great-circle distance)', 100000"

gnuplot plot.gpl
