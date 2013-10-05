#Author - Evan Leister
import vesa_writing_functions as wf
import read_eclipse as re
import eclipse_cells as ec
from time import clock
import sys

def vesa_make_input(layer_id, uniform = False, hydro = False,\
        timestep_days = 1., output_days = 365, simtime_years = 11):
    """  
    """
    print "creating input files"
    output_control = 'all massbalance nowellpressure'
    layers = [layer_id]

    #injwells
    wellid = 1
    ratio = 1.
    layer_id = 1
    l_type = 1
    #massinflow = [[0.0198, 0.0405, 0.0437, 0.0540, 0.0740, 0.1030, \
                  #0.1390, 0.1830, 0.2370, 0.2960, 0.370]]
    massinflow = [0.1418]
    if uniform == True:
        massinflow = [0.1]
    if hydro == True:
        massinflow = [0.]

    simtime_days = [simtime_years * 365]

    # density of C02 [kg/m^3]
    co2_rho = 700.
    #density of brine [kg/m^3]
    brine_rho = 1020.
    #viscosity of CO2 [Pa s]
    co2_mu = 0.0000590
    #viscosity of brine [Pa s ]
    brine_mu = 0.000690
    #residual saturation of C02
    sc_res = 0.00
    #residual saturation of brine
    sb_res = 0.2
    # compressibility of CO2 [1/Pa]
    c_co2 = 0.00000000
    #compressibility of brine
    c_brine = 0.00000000
    # compressibility of rock
    c_rock = 0.00000000

    # capillary pressure stuff: 0 = sharp interface
    cap_rp_id = 0

    if uniform == True:
        nx = 25
        ny = 25
        nz = 1
        dx = 50.
        dy = 50.
        dz = 5.
        center_depth = 862.
        phi = 0.35
        k = 200
        unit = wf.Layer(layers[0], l_type, layer_id, co2_rho, brine_rho, co2_mu,\
                brine_mu, sc_res, sb_res, c_co2, c_brine, c_rock, cap_rp_id,\
                nx, ny, nz = nz, gradient = 10. )
        unit.fill_uniform_grid(dx, dy, dz, center_depth, phi, k)
        xwell = 625.0
        ywell = 625.0
        t_read = clock()
    else:
        e_cells, nx, ny, nz = re.read_eclipse()
        t_read = clock()
        print "Read time"
        print t_read
        unit = wf.Layer(layers[0], l_type, layer_id, co2_rho, brine_rho, co2_mu,\
                brine_mu, sc_res, sb_res, c_co2, c_brine, c_rock, cap_rp_id,\
                nx, ny, nz = nz, gradient = 10. )
        unit.fill_nonuniform_grid(e_cells)
        xwell = 1600.
        ywell = 2057.75

    wf.write_system( timestep_days, output_days, simtime_years, \
            output_control, layers)
    injs = []
    injs.append(wf.Injector(wellid, xwell, ywell, ratio, layer_id, simtime_days, \
            massinflow))
    wf.write_injwells(injs)
    unit.write_layer()

    print "Write Time"
    print clock() - t_read
    print "vesa_make_input COMPLETE"
    return 0

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Please run vesa_make_input in the following way: \n"
        print "$ python vesa_make_input.py <layer_id/simtitle>"
    layer_id = sys.argv[1]
    timestep_days = 0.5
    output_days = 5
    simtime_years = 120 / 365.
    vesa_make_input(layer_id, uniform = True, hydro = False,\
            timestep_days = timestep_days, output_days = output_days,\
            simtime_years = simtime_years)

