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
        for sat in self.saturation:
            self.sat_thickness.append(sat * self.thickness)

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

def plot_vesa_saturation(cells, time_steps, nx, ny, \
        valtype = "saturation", fmt = 'eps', sleipner = False):
    font = { 'size' : '12'}
    matplotlib.rc('font', **font)
    for i in range(0,len(cells[0].get_item_list(valtype))):
        print "Plotting Saturation timestep: " + str(i)
        x, y, zsat = make_plot_grid(cells, i, nx, ny, valtype)
        yrstring = '{:4d}'.format(int(time_steps[i]/365 + 1998))
        if sleipner == True:
            f_sat = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
                facecolor = 'w', edgecolor = 'k')
            f_sat.suptitle("CO2 " + valtype + " in " + yrstring )
        else:
            f_sat = plt.figure(num=None, figsize=(10,10), dpi = 480, \
                facecolor = 'w', edgecolor = 'k')
            f_sat.suptitle("CO2 "+valtype +  " on day: " +\
                    str(time_steps[i]))
        ax_sat = f_sat.add_subplot(111)
        ax_sat.set_xlabel('x-direction [m]')
        ax_sat.set_ylabel('y-direction [m]')
        plt.xticks(np.arange(0,3050,1000))
        n_levels = 21
        if sleipner == True:
            h = 15.
        else:
            h = 5.
        if valtype == 'thickness':
            v_val = np.linspace(0.,h,num=n_levels)
        else: 
            v_val = np.linspace(0.,0.8,num=n_levels)
        cs_sat = ax_sat.contourf(x,y,zsat,v_val) 
        ax_sat.set_aspect('equal')
        cb_sat = plt.colorbar(cs_sat, shrink = 0.8, \
                extend = 'both', ticks =v_val )
        if valtype == 'thickness':
            cb_sat.set_label("CO2 Plume Thickness [m]")
        else:
            cb_sat.set_label("CO2 Saturation []")
        figstr = 'sat_' + '{:02d}'.format(i+1)
        f_sat.savefig(figstr + "." + fmt ,bbox_inches='tight',\
                format = fmt)
        plt.clf()
        plt.close()
        
def val_bounds(cells, valtype):
    # calculates upper and lower bounds for pressure to ensure
    # comparable colorscales across contour plots
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

def plot_timestep_contour(x, y, zval, time, i, valtype, v_val, sleipner, fmt):
    yrstring = '{:4d}'.format(int(time/365 + 1998))
    if sleipner == True:
        f_val = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        f_val.suptitle(valtype + " in " + yrstring )
        plt.xticks(np.arange(0,3050,1000))
    else:
        f_val = plt.figure(num=None, figsize=(10,10), dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        f_val.suptitle(valtype + " on day" + str(time))
    ax_val = f_val.add_subplot(111)
    ax_val.set_xlabel('x-direction [m]')
    ax_val.set_ylabel('y-direction [m]')
    cs_val = ax_val.contourf(x,y,zval,v_val) 
    ax_val.set_aspect('equal')
    cb_val = plt.colorbar(cs_val, shrink = 0.8, \
            extend = 'both', ticks = v_val, format='%.3e' )
    cb_val.set_label(valtype + " [Pa]")
    val_str = valtype + '_' + '{:02d}'.format(i+1)
    f_val.savefig(val_str + "." + fmt ,bbox_inches='tight', \
            format = fmt)
    plt.clf()
    plt.close()

    return 0

def plot_vesa_timesteps(cells, time_steps, nx, ny,  \
        valtype = 'pressure', fmt = 'eps', \
        sleipner = False):
    font = { 'size' : '12'}
    matplotlib.rc('font', **font)

    n_levels = 21
    v_min, v_max = val_bounds(cells, valtype)
    v_val = np.linspace(v_min, v_max, num = n_levels)
    
    for i in range(0,len(cells[0].get_item_list(valtype))):
        print "Plotting " + valtype +  " timestep: " + str(i)
        x, y, zval = make_plot_grid(cells, i, nx, ny, valtype)
        plot_timestep_contour(x, y, zval, time_steps[i], i, valtype, \
             v_val, sleipner, fmt)

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

    # for that cell, get the pressure over time
    if well_head_index != 0:
        pres_list = cells[well_head_index].get_item_list('pressure')
        print pres_list
        print hydro_cells[well_head_index].get_item_list('pressure')[0]
        # add hydrostatic initial conditions
        pres_list.insert(0,\
                hydro_cells[well_head_index].get_item_list('pressure')[0])
        time_steps.insert(0,0.)
        pres = np.asarray(pres_list)/pow(10.,3.)
        time_ar = np.asarray(time_steps)

        f = plt.figure(num=None , dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        f.suptitle('wellhead pressure vs time')
        ax = f.add_subplot(111)
        ax.set_xlabel('time[days]')
        ax.set_ylabel('pressure[kPa]')
        p = plt.plot(time_ar, pres)
        f.savefig('wellhead_pressure' + '.' + fmt)
        plt.clf()
        plt.close()

def plot_cross_sections(cells, time_steps, axis = 1, index = 35, fmt = 'png', \
        uniform = False):
    # make a plot of a vertical slice:
    # include top and bottom boundaries
    # include sharp interface saturation thickness

    if uniform == False:
        nx = 65
    else: 
        nx = 25
    for i in range(0,len(cells[0].get_item_list('saturation'))):
        print "Plotting Cross Section ..." + str(i)
        plume = []
        top = []
        bot = []
        y_list = []
        tempsat = []
        counter = 0
        for cel in cells:
            if (counter % nx ) == 0:
                y_list.append(cel.get_y())
                counter = 0
            if counter == index:
                plume.append(cel.saturation[i] * cel.get_thickness())
                top.append(cel.get_z_top())
                bot.append(cel.get_z_bot())
            counter +=1 

        zb = np.asarray(bot)
        zt = np.asarray(top)
        ys = np.array(y_list)
        zsat = np.asarray(plume)
        plume = zt - zsat
        print len(ys), len(zt), len(zb), len(zsat)

        f = plt.figure(num=None , dpi = 480, \
            facecolor = 'w', edgecolor = 'k')
        title_string = 'Cross section of formation: Axis {0}, Index {1}'.format(axis, index)
        f.suptitle(title_string)
        ax = f.add_subplot(111)
        ax.set_xlabel('distance [m]')
        ax.set_ylabel('[m]')
        p0 = plt.plot(ys, plume, label = "CO2 Plume")
        p1 = plt.plot(ys, zb, label = "Bottom Boundary")
        p2 = plt.plot(ys, zt, label = "Top Boundary")
        plt.legend(loc=4)
        sect_str = '{:02d}'.format(i+1)
        f.savefig('cross_section' + sect_str + '.' + fmt)
        plt.close()
    return 0
