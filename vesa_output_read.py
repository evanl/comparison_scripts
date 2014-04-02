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

    uniform = False
    sleipner = True
    yearwise = True
    if uniform == True:
        nx = 25
        ny = 25
    else:
        nx = 65
        ny = 119
    vr.plot_cross_sections(cells, time_steps, nx, axis = 1, index = nx/2, \
            fmt = 'png', yearwise = yearwise)
    vr.plot_cross_sections(cells, time_steps, nx, axis = 0, index = ny/3 + 2, \
            fmt = 'png', yearwise = yearwise)
    vr.mass_balance_read_print()
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="saturation", fmt = fmt, yearwise = yearwise)
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="pressure", fmt = fmt, yearwise = yearwise)
    vr.plot_vesa_timesteps(cells, time_steps, nx, ny, \
            valtype="delta_p", fmt = fmt, yearwise = yearwise)

    print "Creating directory: " + sim_title 
    call(["mkdir",sim_title])

    hydro_folder = "sl_hydro_700"

    hydro_layer_name = "sl_hydro_700"
    vr.plot_wellhead_pressure(cells, time_steps, hydro_folder, hydro_layer_name,\
            fmt = fmt, sleipner = False )

    move_files(fmt, sim_title)
    print "total time processing vesa"
    print clock() - t_read
   
