#!/bin/bash

dir=$PWD
parentdir="$(dirname "$dir")"
rsync -t $parentdir/tough/t2/*.py $dir/
rsync -t $parentdir/vesa/VESA_02_01_13/*.py $dir/
rm CO2Properties_1_0.py
rm __init__.py
