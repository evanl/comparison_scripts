#Author - Evan Leister
import tough.t2.process_t2_output as pt2
import vesa.VESA_02_01_13.vesa_output_read as pv
from time import time, clock
import sys
import os
from subprocess import call

def compare(vesa_simtitle, tough_simtitle, \
        vesa_hydro_folder = 'unif_hydro'):
    print "------------------------------------------------------------\n"
    print "Comparing "
    print "  VESA simulation: " + vesa_simtitle
    print "with"
    print "TOUGH2 simulation: " + tough_simtitle
    print "------------------------------------------------------------\n"
    current_directory = os.getcwd()
    vesa_dirname =  current_directory + "/vesa/VESA_02_01_13/" + vesa_simtitle
    tough_dirname = current_directory + "/tough/t2/" + tough_simtitle + "_dir/"

    os.chdir(tough_dirname)
    t2_grid, t2_timesteps = pt2.process_t2_output(tough_simtitle)
    # this reads the output for a given unit, but will not work for 
    # multiple units. the positions will have to be combined for that later
    os.chdir(vesa_dirname)
    vesa_cells, vesa_timesteps = pv.read_output_data(layer = vesa_simtitle)

    os.chdir(current_directory)
    print "yay"
    return 0

if __name__ == '__main__': 
    if len(sys.argv) != 3:
        sys.exit("Proper usage for compare.py is \n" +
                "python compare.py <vesa_simtitle> <tough_simtitle>\n")
    vesa_simtitle = sys.argv[1]
    tough_simtitle = sys.argv[2]
    compare(vesa_simtitle, tough_simtitle)
