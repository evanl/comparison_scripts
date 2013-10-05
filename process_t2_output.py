#Author - Evan Leister
import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import output_t2_funcs as t2o

def process_t2_output(sim_title ):
    # this must be changed to reflect the input data
    # create an output grid object
    grid = t2o.T2grid()
    # read the coordinate data from MESH
    grid.read_mesh()
    # count the number of outputs
    num_outputs = 0 
    f = open(sim_title + '.out')
    a = f.read()
    b = a.split()
    for char in b:
        if char == 'OUTPUT':
            num_outputs +=1
    print "The number of timesteps output: " + str(num_outputs)
    f.close()

    # read the file for the timestep data
    f = open(sim_title + '.out')
    grid.read_initial_mass(f)

    time_steps = []
    for i in range(num_outputs):
        timetest = t2o.T2Timestep()
        timetest.readOutput(f, grid)
        timetest.get_co2_density_viscosity()
        time_steps.append(timetest)
    f.close()

    return grid, time_steps

if __name__ == '__main__':
    print "PROCESSING TOUGH2 OUTPUT FOR " + str(sys.argv[1])
    if len(sys.argv) != 2:
        sys.exit( "Please specify a simulation title")
    sim_title = sys.argv[1]
    hydro = False
    two_d = False
    sleipner = False
    section = False
    shale = True
    # creates a T2Grid object and a list of T2Timestep objects
    grid, time_steps = process_t2_output(sim_title)
    # choose format for plots.
    fmt = 'png'
    if two_d == True:
        eleme = 'aH732'
        k_layer = 0
        xind = 65/2
        yind = 119/2
    else:
        # sleipner
        if sleipner == True:
            eleme = 'IH732'
            k_layer = 3
            xind = 32
            yind = 77
            if section == True:
                eleme = 'AE717'
                xind = 17
                yind = 47
                k_layer = 26
        else:
            eleme = 'mB212'
            k_layer = 24
            xind = 12
            yind = 12
        #modified for center
        #eleme = 'cA2 2'
        #k_layer = 2
        # column
        #eleme = 'dA0 0'
        #k_layer = 0
        # box

    for basis in range(1,4):
        t2o.plot_planar_contours(grid, time_steps, sim_title, fmt, \
                two_d = two_d, sleipner = sleipner, section = section,\
                shale = shale, axis = basis,\
                i_in = xind, j_in = yind, k_in = k_layer)

    t2o.plot_wellhead_pressure(grid, time_steps, \
            well_cell_id = eleme)

    if hydro == True:
        t2o.plot_incon_change_vs_index(grid, time_steps, vs_elev = True)
        t2o.check_3d_hydro_pressure(grid, time_steps)

    t2o.plot_mass_balance(grid, time_steps)

    t2o.write_viscosity(grid, time_steps)

    sim_title_dir = sim_title + '_dir'
    call(["mkdir",sim_title_dir])
    movestring = "mv " + "*" + fmt +" " + sim_title_dir
    call(movestring, shell=True)
    call (["mv",sim_title,sim_title_dir])
    call (["mv",sim_title+'.out',sim_title_dir])
    call (["mv",'INCON' ,sim_title_dir])
    call (["mv",'MESH' ,sim_title_dir])
    call (["mv",'SAVE' ,sim_title_dir])
    call (["mv",'rho_visc' ,sim_title_dir])
    call (["cp",'CO2TAB' ,sim_title_dir])
