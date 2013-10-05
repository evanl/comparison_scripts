#Author - Evan Leister
import sys
import string
from time import clock
import read_eclipse as re
import eclipse_cells as ec
import input_t2_funcs as it2f

def create_t2_input(sim_title, two_d = False, uniform = False,\
        sleipner = False, hydro = False, hydro_directory = False,\
        num_steps = 11, days_per_step = 365.25, fs = [],\
        bc_type = 1, column = False, linear_rp = False,
        shale = True):
    print 'CREATING TOUGH2 INPUT FILE'
    if two_d == True:
        eleme = 'aH732'
    else:
        eleme = 'IH732'
        if fs_sec != []:
            eleme = 'AE717'

    if uniform == True:
        #eleme = 'cA2 2'
        #eleme = 'dA0 0'
        eleme = 'yB212'
        nx = 25
        ny = 25
        nz = 25
        dx = 50
        dy = 50
        dz = 0.2
    else:
        e_cel, nx, ny, nz = re.read_eclipse()

    print "Read time"
    t_read = clock()
    print t_read

    output_day_list = []
    #output_day_list.append(30.)
    for i in range(num_steps):
        output_day_list.append(days_per_step * (i+1))

    print "nx, ny, nz"
    print nx, ny, nz

    print 'Creating T2 input files for simulation: ' + sim_title
    f = open(sim_title,'w')
    f.write('*'+ sim_title +'* test, 2d-case\n')

    #input parameters
    name = 'sands'
    sand_density = 2600.
    if sleipner == True:
        porosity = 0.35 
        xperm = 2.e-12
        yperm = 2.e-12
        zperm = 2.e-12
    else:
        porosity = 0.35
        xperm = 2.e-13
        yperm = xperm
        zperm = xperm / 100.

    cap = 'vanGenuchten'
    if linear_rp == False:
        rel_perm = 'vanGenuchten'
    else:
        rel_perm = 'linear'

    if cap == 'vanGenuchten':
        cp_vals = [0.4, 0.2, 1.61e-5, 1.e7, 0.999]
    else:
        cp_vals = [2.0e4, 0.0, 1.0]# linear

    if rel_perm == 'vanGenuchten':
        rp_vals = [0.4, 0.2, 1.0, 0.05]
    else:
        rp_vals = [0.2, 0., 1., 1.]# linear

    it2f.plot_relperm_cap(rp_vals, cp_vals, fmt = 'png',\
            rp = rel_perm, cp = cap)

    it2f.write_separator(f, 'ROCKS')
    it2f.write_rocks(f, name, sand_density, porosity, xperm, yperm, zperm, \
        cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
        cap = cap, rel_perm = rel_perm)
    name = 'shale'
    shale_density = 2600.
    porosity = 0.35
    xperm = 1.e-18
    yperm = 1.e-18
    zperm = 1.e-18
    it2f.write_rocks(f, name, shale_density, porosity, xperm, yperm, zperm, \
        rp_vals, cp_vals, thermk = 2.51, specheat = 920., \
        cap = cap, rel_perm = rel_perm,\
        end = True)

    # Cavanagh test case
    # output_day_list.append(365.25 * 50)

    it2f.write_multi(f)
    it2f.write_selec(f)
    it2f.write_start(f)
    it2f.write_param(f, tmax = output_day_list[-1] * 24 * 3600)
    it2f.write_solvr(f)

    phase = 'co2'
    if hydro == True:
        mass_rate = 0.0
    else:
        if sleipner == True:
            mass_rate = 4.475
        else:
            mass_rate = 3.17
    # in megatons per year
    # estimated flow rate.
    massinflow = [0.0198, 0.0405, 0.0437, 0.0540, 0.0740, 0.1030, \
                                0.1390, 0.1830, 0.2370, 0.2960, 0.370]
    kg_per_sec_factor = 31.71 # kg/s per mt/yr
    kgInflow = [ x * kg_per_sec_factor for x in massinflow]

    it2f.write_gener(f, eleme, phase = phase, mass_rate = mass_rate, \
            column = column)
            #, kg_inflow = kgInflow, times = output_day_list )

    it2f.write_times(f, output_day_list)
    it2f.write_foft(f)
    it2f.write_coft(f)
    it2f.write_goft(f)
    it2f.write_separator(f, 'ENDCY')
    f.close()

    # f = open('mrad1','w')
    # f.write('*mrad1* Make mesh for rad1layer\n')
    # it2f.write_meshmaker(f , flat = True, nx = 5, Ny = 5, nz = 3)
    # it2f.write_separator(f, 'ENDFI')
    # f.close()

    # create an input grid instance
    tg = it2f.T2InputGrid(nx, ny, nz)

    solubility = 0.454104e-3
    #solubility = 0.0

    if uniform == True:
        brine_density = 1016.4
        altered_cell = 'eA2 2'
        altered_cell = 'none'
        tg.fill_uniform_grid(porosity, dx, dy, dz, density = brine_density,\
                solubility_limit = solubility,\
                altered_cell = altered_cell)
        e_cel = 'uniform'
        tg.write_mesh(e_cel, two_d = two_d, uniform = uniform,\
                boundary_type = bc_type, shale = shale)
    else:
        brine_density = 1019.35
        tg.fill_3d_grid(e_cel, temperature = 32., density = brine_density,\
                two_d = two_d, solubility = solubility,\
                five_section = fs, shale = shale)
        tg.write_mesh(e_cel, two_d = two_d, uniform = uniform,\
                boundary_type = bc_type, shale = shale)


    if hydro_directory == False:
        tg.write_incon(porosity)
    else:
        tg.use_old_incon(hydro_directory)

    print "Write Time"
    print clock() - t_read
    tg.plot_cells(show = False, two_d = two_d)
    tg.plot_cells_slice(direc = 2, ind = 0, show = False)
    return "create_t2_input COMPLETE"

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit( "Please specify a simulation title")
    sim_title = str(sys.argv[1])

    # for creating column of injection for VESA comparison
    # if column is not needed, just make it false
    endstr = 'B212'
    letters = string.lowercase
    column = []
    for i in range(25):
        column.append(letters[i] + endstr)
    column = False

    # start ijk values for the section of sleipner data
    i_start = 15
    j_start = 45 
    k_start = 24
    fs_sec = [i_start, j_start, k_start]
    # if section is desired, comment this line
    fs_sec = []

    # if hd == False, it means that the initial condition will be generated 
    # using density estimates. If a directory is specified, the 
    # initial pressures and dissolved fractions will be taken from
    # the 'hd' + '_dir/'
    sleipner = False
    shale = True
    hydro = False
    uniform = True
    
    # sanity check
    if hydro == False:
        hd = False
    else:
        hd = 'unif_hydro'
    if uniform == True:
        shale = True
        sleipner = False
    print create_t2_input(sim_title, two_d = False, uniform = uniform, \
            sleipner = sleipner, hydro = hydro, hydro_directory = hd, \
            num_steps = 24, days_per_step = 5, fs = fs_sec,\
            bc_type = 1, column = column, linear_rp = False,\
            shale = shale)

