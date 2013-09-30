#Author - Evan Leister
import vesa_reading_functions as vr
from subprocess import call
from time import time, clock
import glob

def move_files(fmt, sim_title):
    print ('moving all ' + fmt + ' files to the directory ' + sim_title)
    call(movestring, shell=True)
    movestring = "cp " + "*" + fmt + " " + sim_title
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
    layer = 'sl_const_iso'
    cells, time_steps = vr.read_output_data(layer = layer)
    t_read = clock()

    uniform = False
    sleipner = True
    if uniform == True:
        nx = 25
        ny = 25
    else:
        nx = 65
        ny = 119
    vr.plot_cross_sections(cells, time_steps, axis = 1, index = nx/2, fmt = 'png',
            uniform = uniform)

    vr.mass_balance_read_print()

    sim_title = 'sl_const_iso'
    fmt = 'png'
    vr.plot_vesa_saturation(cells, time_steps, nx, ny, \
            plottype="saturation", fmt = fmt, \
            cbar = True, frame = True, 
            thickness = False, sleipner = sleipner)
    vr.plot_vesa_pressure(cells, time_steps, nx, ny, \
            plottype="pressure", fmt = fmt, \
            cbar = True, frame = True)
    print "Creating directory: " + sim_title 
    call(["mkdir",sim_title])

    hydro_folder = "unif_hydro"
    hydro_layer_name = 'unif_hydro'
    vr.plot_wellhead_pressure(cells, time_steps, hydro_folder, hydro_layer_name,\
            fmt = fmt, sleipner = False )

    move_files(fmt, sim_title)
    print "total time processing vesa"
    print clock() - t_read
   
