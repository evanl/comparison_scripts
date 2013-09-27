#!/bin/bash
dir=$PWD
parentdir="$(dirname "$dir")"
rsync -t *t2* $parentdir/tough/t2/ 
rsync -t *eclipse* $parentdir/tough/t2/ 

rsync -t vesa* $parentdir/vesa/VESA_02_01_13/ 
rsync -t *eclipse* $parentdir/vesa/VESA_02_01_13/
