#Author - Evan Leister
import vesa_reading_functions as vr
from subprocess import call
import sys
from time import time, clock
import glob

def move_files(fmt, sim_title):
    print ('moving all ' + fmt + ' files to the directory ' + sim_title)
    movestring = "cp " + "*" + fmt + " " + sim_title
    call(movestring, shell=True)
    call (["cp","massBalanceError.txt",sim_title])
    call (["cp",layer + ".txt",sim_title])
    call (["cp","System.txt",sim_title])
    call (["cp","InjWells.txt",sim_title])
    call (["cp","p.csv",sim_title])
    call (["cp","scbar.csv",sim_title])
    call (["cp","MassBalance.csv",sim_title])
    call (["cp","thickness.txt",sim_title])
    return 0

if __name__ == '__main__':
    if len(sys.argv) == 2:
        layer = sys.argv[1]
        sim_title = sys.argv[1]
    elif len(sys.argv) == 3:
        layer = sys.argv[2]
        sim_title = sys.argv[1]
    else:
        print "please run vesa_output_read in the following way \n"
        print "$ python vesa_output_read.py <sim_title> (optional: <layer_id>)"
        print "\n"

    fmt = 'png'
    
    cells, time_steps = vr.read_output_data(layer = layer)
    t_read = clock()

    uniform = True
    sleipner = False
    if uniform == True:
        nx = 25
        ny = 25
    else:
        nx = 65
        ny = 119
    vr.plot_cross_sections(cells, time_steps, nx, axis = 1, index = nx/2, \
            fmt = 'png')
    vr.mass_balance_read_print()
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="saturation", fmt = fmt)
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="pressure", fmt = fmt, )
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="delta_p", fmt = fmt, )

    print "Creating directory: " + sim_title 
    call(["mkdir",sim_title])

    hydro_folder = "unif_hydro"
    hydro_layer_name = "unif_hydro"
    vr.plot_wellhead_pressure(cells, time_steps, hydro_folder, hydro_layer_name,\
            fmt = fmt, sleipner = False )

    move_files(fmt, sim_title)
    print "total time processing vesa"
    print clock() - t_read
   
