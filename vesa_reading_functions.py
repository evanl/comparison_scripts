import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import read_eclipse as re
import eclipse_cells as ec
from matplotlib import cm

class Cell:
    def __init__(self, x, y, z_top, z_bot, p_init):
        self.x_center = x
        self.y_center = y
        self.z_top = z_top
        self.z_bot = z_bot
        self.p_init = p_init
        self.thickness = self.z_top - self.z_bot
        self.saturation = []
        self.pressure = []
        self.delta_p = []
        self.sat_thickness = []

    def __str__(self):
        a = str(self.get_x()) + ' '
        b = str(self.get_y()) + ' '
        c = str(self.get_z_top()) + ' ' 
        d = str(self.get_z_bot()) + ' '
        return a + b + c + d

    def get_x(self):
        return self.x_center

    def get_y(self):
        return self.y_center

    def get_z_top(self):
        return self.z_top

    def get_thickness(self):
        return self.z_top - self.z_bot

    def get_z_bot(self):
        return self.z_bot

    def get_item_list(self, valtype):
        if valtype == 'pressure':
            cvals = self.pressure
        elif valtype == 'saturation':
            cvals = self.saturation
        elif valtype == 'thickness':
            cvals = self.sat_thickness
        elif valtype == 'delta_p':
            cvals = self.delta_p
        else: 
            print "please specify a valid list"
            return 1
        return cvals

def read_output_data(layer = 'SleipnerL9'):
    print "Reading VESA output from layer: " + layer
    # reads the pressure file
    r = open('p.csv','r')

    # line 1 - x values
    line = r.readline()
    tempx = line.split(', ')
    tempx.remove('')
    tempx[-1] = tempx[-1].replace(' \n','')

    # line 2 - y values
    line = r.readline()
    tempy = line.split(', ')
    tempy.remove('')
    tempy[-1] = tempy[-1].replace(' \n','')

    #get boundaries, initial pressure
    z_top = []
    z_bot = []
    p_init = []
    fin = open(layer + '.txt','r')
    for i in range(14):
        line = fin.readline()
    while line:
        s = line.split(',')
        if s[-1] == '\n':
            s[-1] = s[-1].replace( '\n','')
        if s[-1] == '':
            s.remove('')
        z_bot.append(s[8])
        z_top.append(s[9])
        p_init.append(s[12])
        line = fin.readline()

    # initialize all cells
    cells = []
    for i in range(0,len(tempx)):
        c_in = Cell(float(tempx[i]), float(tempy[i]), \
                float(z_top[i]), float(z_bot[i]), float(p_init[i]))
        cells.append(c_in)
    # line 3 - layer ID
    line = r.readline()
    # line 4 to end - pressures
    line = r.readline()
    templine = line.split(',')
    baseline = templine
    baseline[-1] = templine[-1].replace(' \n','')
    time_steps = []
    while line:
        time_steps.append(float(templine[0]))
        templine[-1] = templine[-1].replace(' \n','')
        for i in range(1,len(templine)):
            cells[i-1].pressure.append(float(templine[i]))
            cells[i-1].delta_p.append(float(templine[i]) - cells[i-1].p_init)
        line = r.readline()
        templine = line.split(',')
    r.close()
    fin.close()
    # reads in the saturation values
    r = open('scbar.csv','r')
    # skip the first 3 lines of the text file.
    # information has already been read
    for i in range(0,4):
        line = r.readline()
    while line:
        templine = line.split(',')
        del templine[0]
        templine[-1] = templine[-1].replace('\n','')
        for i in range(0,len(templine)):
            cells[i].saturation.append(float(templine[i]))
            cells[i].sat_thickness.append(\
                    float(templine[i]) * cells[i].get_thickness() /\
                    (1 - 0.2))
        line = r.readline()

    r.close()
    return cells , time_steps

def mass_balance_read_print():
    r = open("MassBalance.csv",'r')
    line = r.readline()
    line = r.readline()
    years_mass = []
    injected_mass = []
    actual_mass = []
    boundary_mass= []
    while line:
        line = line.split(',')
        years_mass.append(float(line[0]))
        actual_mass.append(float(line[1]))
        injected_mass.append(float(line[2]))
        boundary_mass.append(float(line[3]))
        line = r.readline()
    perdiff_mass = []
    for i in range(len(injected_mass)):
        if injected_mass[i] != 0.:
            perdiff_mass.append(\
                (actual_mass[i] - injected_mass[i]+ boundary_mass[i]) \
                / injected_mass[i] * 100)
        else: 
            perdiff_mass.append(actual_mass[i] - injected_mass[i] + \
                    boundary_mass[i])
    print 'year  | mass balance error'
    for j in range(0,len(perdiff_mass)):
        print j+1, '    | '+str(perdiff_mass[j])

    f = open("massBalanceError.txt",'w')
    f.write("Year       | Mass Balance Error \n" )
    for j in range(0,len(perdiff_mass)):
        f.write(str(j+1) + "     | " + str(perdiff_mass[j]) + "\n")
    f.close()
        
def val_bounds(cells, valtype):
    # calculates upper and lower bounds for pressure to ensure
    # comparable colorscales across contour plots
    if valtype == 'saturation':
        v_min = 0.0
        v_max = 0.8
    elif valtype == 'thickness':
        v_min = 0.0
        v_max = 15.
    else:
        val_bound = []
        for c in cells:
            cvals = c.get_item_list(valtype)
            for el in cvals:
                val_bound.append(el)
        v_min = min(val_bound)
        v_max = max(val_bound)
    return v_min, v_max

def make_plot_grid(cells, time_index, nx, ny, valtype):
    val = []
    x_list = []
    y_list = []
    temp_val = []
    counter = 0
    for cel in cells:
        if counter < nx: 
            x_list.append(cel.get_x())
        temp_val.append(cel.get_item_list(valtype)[time_index])
        counter +=1 
        if (counter % nx ) == 0:
            y_list.append(cel.get_y())
            val.append(temp_val)
            temp_val = [] 
    xs = np.array(x_list)
    ys = np.array(y_list)
    x, y = np.meshgrid(xs, ys)
    zval = np.asarray(val)
    return x, y, zval

def plot_timestep_contour(x, y, zval, time, i, valtype, v_val, fmt,\
        yearwise = False):
    yrstring = '{:4d}'.format(int(time/365 + 1998))
    f_val = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    ax_val = f_val.add_subplot(111)
    ax_val.set_xlabel('x-direction [m]')
    ax_val.set_ylabel('y-direction [m]')
    cs_val = ax_val.contourf(x,y,zval,v_val) 
    ax_val.set_aspect('equal')
    if valtype =='pressure':
        clab = valtype + " [Pa]"
        cfm = '%.3e'
    elif valtype == 'delta_p':
        clab = valtype + " [Pa]"
        cfm = '%.3e'
    elif valtype == 'thickness':
        clab = "CO2 Plume Thickness [m]"
        cfm = '%.1f'
    elif valtype == 'saturation':
        clab = "CO2 Saturation []"
        cfm = '%.2f'
    cb_val = plt.colorbar(cs_val, shrink = 0.8, \
            extend = 'both', ticks = v_val, format=cfm)
    cb_val.set_label(clab)
    if yearwise == True:
        val_str = valtype + '_' + yrstring
    else:
        val_str = valtype + '_' + '{:02d}'.format(i+1)
    f_val.suptitle(val_str)
    f_val.savefig(val_str + "." + fmt ,bbox_inches='tight', \
            format = fmt)
    plt.clf()
    plt.close()
    return 0

def plot_vesa_timesteps(cells, time_steps, nx, ny,  \
        valtype = 'pressure', fmt = 'eps', yearwise = False):
    font = { 'size' : '12'}
    matplotlib.rc('font', **font)

    n_levels = 21
    v_min, v_max = val_bounds(cells, valtype)
    v_val = np.linspace(v_min, v_max, num = n_levels)
    
    for time_index in range(0,len(cells[0].get_item_list(valtype))):
        print "Plotting " + valtype +  " timestep: " + str(time_index)
        x, y, zval = make_plot_grid(cells, time_index, nx, ny, valtype)
        plot_timestep_contour(x, y, zval, time_steps[time_index], \
                time_index, valtype, v_val, fmt, yearwise = yearwise)

def plot_wellhead_pressure(cells, time_steps, hydro_directory, hydro_layer_name, \
        x_well = 1600., y_well = 2057.75,\
         fmt = 'png', sleipner = True):
    print "Plotting Wellhead Pressure..."
    # find cell with that is the wellhead cell
    if sleipner == True:
        well_head_index = 2697
    else:
        well_head_index = 313
    os.chdir(hydro_directory + '/')
    hydro_cells, hydro_time_steps = read_output_data(layer = hydro_layer_name)
    os.chdir('../')

    font = { 'size' : '16'}
    matplotlib.rc('font', **font)
    # for that cell, get the pressure over time
    if well_head_index != 0:
        pres_list = cells[well_head_index].get_item_list('pressure')
        # add hydrostatic initial conditions
        pres_list.insert(0,\
                hydro_cells[well_head_index].get_item_list('pressure')[0])
        print pres_list
        time_steps.insert(0,0.)
        pres = np.asarray(pres_list)
        pres = pres/pow(10.,3)
        time_ar = np.asarray(time_steps)

        f = plt.figure(num=None , dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        ax = f.add_subplot(111)
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Wellhead Pressure [kPa]')
        p = plt.plot(time_ar, pres)
        f.savefig('wellhead_pressure' + '.' + fmt)
        plt.clf()
        plt.close()

def make_cross_sections(cells, time_index, axis, index, nx):
    plume = []
    top = []
    bot = []
    y_list = []
    tempsat = []
    counter = 0
    if axis == 1:
        for cel in cells:
            if (counter % nx ) == 0:
                y_list.append(cel.get_y())
                counter = 0
            if counter == index:
                plume.append(cel.get_item_list('thickness')[time_index])
                top.append(cel.get_z_top())
                bot.append(cel.get_z_bot())
            counter +=1 
    elif axis == 0:
        for cel in cells:
            if counter / nx == 0:
                y_list.append(cel.get_x())
            if counter / nx == index:
                plume.append(cel.get_item_list('thickness')[time_index])
                top.append(cel.get_z_top())
                bot.append(cel.get_z_bot())
            counter +=1

    zb = np.asarray(bot)
    zt = np.asarray(top)
    ys = np.array(y_list)
    zsat = np.asarray(plume)
    plume = zt - zsat
    if len(ys) != len(zt) != len(zb) != len(zsat):
        print "NONEQUAL LENGTH ARRAYS"
        return 1

    return ys, zb, zt, plume

def plot_cross_sections(cells, time_steps, nx, axis = 2, index = 32,\
        fmt = 'png', yearwise = False):
    # make a plot of a vertical slice:
    # include top and bottom boundaries
    # include sharp interface saturation thickness
    print "making cross section with axis: " + str(axis) + ", index: " +\
            str(index)
    for i in range(0,len(cells[0].get_item_list('saturation'))):
        print "Plotting Cross Section ..." + str(i)
        font = { 'size' : '16'}
        matplotlib.rc('font', **font)
        ys, zb, zt, plume = make_cross_sections(cells, i, axis, index, nx)
        f = plt.figure(num=None , dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        #title_string = \
            #'Cross section of formation: Axis {0}, Index {1}: '.format(axis, index)
        title_string = ''
        yrstring = '{:4d}'.format(int(time_steps[i]/365 + 1998))
        if yearwise == True:
            title_string += ' in ' + yrstring
        else:
            title_string += 'Time t = {0} days'.format(time_steps[i])
        f.suptitle(title_string)
        ax = f.add_subplot(111)
        ax.set_xlabel('Lateral Distance [m]')
        ax.set_ylabel('Elevation [m]')
        p0 = plt.plot(ys, plume, label = "CO2 Thickness")
        p1 = plt.plot(ys, zb, label = "Bottom Boundary")
        p2 = plt.plot(ys, zt, label = "Top Boundary")
        plt.legend(loc=4)
        sect_str = '{:02d}'.format(i+1)
        f.savefig('cross_section_' + str(axis) + '_' + str(index) + '_' \
                + sect_str + '.' + fmt)
        plt.close()
    return 0
