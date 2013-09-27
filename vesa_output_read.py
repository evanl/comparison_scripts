# This script reads VESA output files 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import read_eclipse as re
import eclipse_cells as ec
from matplotlib import cm
from subprocess import call
import os
import glob
from time import time, clock

class Cell:
    def __init__(self,x,y,z_top,z_bot):
        self.x_center = x
        self.y_center = y
        self.z_top = z_top
        self.z_bot = z_bot
        self.thickness = self.z_top - self.z_bot
        self.saturation = []
        self.pressure = []
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

    #get top and bottom boundaries
    z_top = []
    z_bot = []
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
        line = fin.readline()


    # initialize all cells
    cells = []
    for i in range(0,len(tempx)):
        c_in = Cell(float(tempx[i]), float(tempy[i]), \
                float(z_top[i]), float(z_bot[i]))
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
            # check if the hydrostatic case works
            # if float(templine[i]) != float(baseline[i]):
            #     print i, templine[i], baseline[i]
        line = r.readline()
        templine = line.split(',')

    # reads in the saturation values

    r.close()
    fin.close()

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
        plottype = "Contour", fmt = 'eps', cbar = True, frame = False,\
        thickness = True, sleipner = False):
    font = { 'size' : '12'}
    matplotlib.rc('font', **font)
    for i in range(0,len(cells[0].get_item_list('saturation'))):
        sat = []
        sat_thick = []
        x_list = []
        y_list = []
        tempsat = []
        temptt = []
        counter = 0
        for cel in cells:
            if counter < nx: 
                x_list.append(cel.get_x())
            tempsat.append(cel.saturation[i] )
            temptt.append(cel.saturation[i] * cel.get_thickness())
            counter +=1 
            if (counter % nx ) == 0:
                y_list.append(cel.get_y())
                sat.append(tempsat)
                sat_thick.append(temptt)
                tempsat = [] 
                temptt = []
        xs = np.array(x_list)
        ys = np.array(y_list)
        X, Y = np.meshgrid(xs, ys)
        if thickness == True:
            zsat = np.asarray(sat_thick)
        else:
            zsat = np.asarray(sat)
        if frame == True:
            print "Frame"
            yrstring = '{:4d}'.format(int(time_steps[i]/365 + 1998))
            if sleipner == True:
                f_sat = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
                    facecolor = 'w', edgecolor = 'k')
                if thickness == True:
                    f_sat.suptitle("CO2 Plume thickness in " + yrstring )
                else:
                    f_sat.suptitle("CO2 Saturation in " + yrstring )
            else:
                f_sat = plt.figure(num=None, figsize=(10,10), dpi = 480, \
                    facecolor = 'w', edgecolor = 'k')
                if thickness == True:
                    f_sat.suptitle("CO2 Plume Thickness on day: " +\
                            str(time_steps[i]))
                else:
                    f_sat.suptitle("CO2 Saturation on day: " +\
                            str(time_steps[i]))
            ax_sat = f_sat.add_subplot(111)
            ax_sat.set_xlabel('x-direction [m]')
            ax_sat.set_ylabel('y-direction [m]')
            plt.xticks(np.arange(0,3050,1000))
            # plt.yticks(rotation=90)
            plotframe = 0
            
        else:
            print "notframe"
            f_sat = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
                facecolor = 'w', edgecolor = 'k', frameon = False)
            ax_sat = f_sat.add_subplot(111)
            ax_sat.set_axis_off()
            plotframe = 'tight'
        print "Plotting Saturation timestep: " + str(i)
        n_levels = 21
        if sleipner == True:
            h = 15.
        else:
            h = 5.
        if thickness == True:
            v_sat = np.linspace(0.,h,num=n_levels)
        else: 
            v_sat = np.linspace(0.,0.8,num=n_levels)
        cs_sat = ax_sat.contourf(X,Y,zsat,v_sat) 
        ax_sat.set_aspect('equal')
        if cbar == True:
            cb_sat = plt.colorbar(cs_sat, shrink = 0.8, \
                    extend = 'both', ticks =v_sat )
            if thickness == True:
                cb_sat.set_label("CO2 Plume Thickness [m]")
            else:
                cb_sat.set_label("CO2 Saturation []")
        figstr = 'sat_' + '{:02d}'.format(i+1)
        f_sat.savefig(figstr + "." + fmt ,bbox_inches=plotframe,\
                format = fmt)
        plt.clf()
        plt.close()


def plot_vesa_pressure(cells, time_steps, nx, ny,  \
        plottype = 'pressure', fmt = 'eps', cbar = True, frame = False, \
        sleipner = False):
    font = { 'size' : '12'}
    matplotlib.rc('font', **font)

    # calculates upper and lower bounds for pressure to ensure
    # comparable colorscales across contour plots
    val_bound = []
    for c in cells:
        cvals = c.get_item_list('pressure')
        for el in cvals:
            val_bound.append(el)

    p_min = min(val_bound)/pow(10.,6.)
    p_max = max(val_bound)/pow(10.,6.)
    n_levels = 21

    for i in range(0,len(cells[0].get_item_list('pressure'))):
        val = []
        x_list = []
        y_list = []
        temp_val = []
        counter = 0
        for cel in cells:
            if counter < nx: 
                x_list.append(cel.get_x())
            temp_val.append(cel.pressure[i]/pow(10.,6.) )
            counter +=1 
            if (counter % nx ) == 0:
                y_list.append(cel.get_y())
                val.append(temp_val)
                temp_val = [] 
        xs = np.array(x_list)
        ys = np.array(y_list)
        X, Y = np.meshgrid(xs, ys)
        zval = np.asarray(val)
        if frame == True:
            print "Frame"
            yrstring = '{:4d}'.format(int(time_steps[i]/365 + 1998))
            if sleipner == True:
                f_val = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
                    facecolor = 'w', edgecolor = 'k')
                f_val.suptitle("Bottom Pressure in " + yrstring )
            else:
                f_val = plt.figure(num=None, figsize=(10,10), dpi = 480, \
                    facecolor = 'w', edgecolor = 'k')
                f_val.suptitle("Bottom Pressure on day" + str(time_steps[i]))
            ax_val = f_val.add_subplot(111)
            ax_val.set_xlabel('x-direction [m]')
            ax_val.set_ylabel('y-direction [m]')
            plt.xticks(np.arange(0,3050,1000))
            # plt.yticks(rotation=90)
            plotframe = 0
        else:
            print "notframe"
            f_val = plt.figure(num=None, figsize=(7.5,10), dpi = 480, \
                facecolor = 'w', edgecolor = 'k', frameon = False)
            ax_val = f_val.add_subplot(111)
            ax_val.set_axis_off()
            plotframe = 'tight'
        print "Plotting Pressure timestep: " + str(i)
        v_val = np.linspace(p_min, p_max, num = n_levels)
        cs_val = ax_val.contourf(X,Y,zval,v_val) 
        ax_val.set_aspect('equal')
        if cbar == True:
            cb_val = plt.colorbar(cs_val, shrink = 0.8, \
                    extend = 'both', ticks = v_val )
            cb_val.set_label("Pressure [MPa]")
        val_str = 'pres_' + '{:02d}'.format(i+1)
        f_val.savefig(val_str + "." + fmt ,bbox_inches=plotframe, \
                format = fmt)
        plt.clf()
        plt.close()

def plot_wellhead_pressure(cells, time_steps, hydro_directory, hydro_layer_name, \
        x_well = 1600., y_well = 2057.75,\
         fmt = 'png', sleipner = True):
    print "Plotting Wellhead Pressure..."
    # find cell with that is the wellhead cell
    if sleipner == True:
        well_head_index = 2697
    else:
        well_head_index = 313
    # for i in range(len(cells)):
    #     if (cells[i].get_x() == x_well) and \
    #             (cells[i].get_y() == y_well):
    #         print i, cells[i].get_x(), cells[i].get_y()
    #         well_head_index = i

    # get initial timestep from a hydrostatic dataset
    os.chdir(hydro_directory + '/')
    hydro_cells, hydro_time_steps = read_output_data(layer = hydro_layer_name)
    os.chdir('../')

    # for that cell, get the pressure over time
    if well_head_index != 0:
        pres_list = cells[well_head_index].get_item_list('pressure')
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
        sat = []
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
                sat.append(cel.saturation[i] * cel.get_thickness())
                top.append(cel.get_z_top())
                bot.append(cel.get_z_bot())
            counter +=1 

        zb = np.asarray(bot)
        zt = np.asarray(top)
        ys = np.array(y_list)
        zsat = np.asarray(sat)
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

if __name__ == '__main__':
    layer = 'sl_const'
    cells, time_steps = read_output_data(layer = layer)
    t_read = clock()

    uniform = False
    sleipner = True
    if uniform == True:
        nx = 25
        ny = 25
    else:
        nx = 65
        ny = 119
    plot_cross_sections(cells, time_steps, axis = 1, index = nx/2, fmt = 'png',
            uniform = uniform)

    mass_balance_read_print()

    sim_title = 'sl_const'
    fmt = 'png'
    plot_vesa_saturation(cells, time_steps, nx, ny, \
            plottype="saturation", fmt = fmt, \
            cbar = True, frame = True, 
            thickness = False, sleipner = sleipner)
    plot_vesa_pressure(cells, time_steps, nx, ny, \
            plottype="pressure", fmt = fmt, \
            cbar = True, frame = True)
    print "Creating directory: " + sim_title 
    call(["mkdir",sim_title])
    movestring = "cp " + "*" + fmt + " " + sim_title

    hydro_folder = "unif_hydro"
    hydro_layer_name = 'unif_hydro'
    plot_wellhead_pressure(cells, time_steps, hydro_folder, hydro_layer_name,\
            fmt = fmt, sleipner = False )

    print ('moving all ' + fmt + ' files to the directory ' + sim_title)
    call(movestring, shell=True)
    call (["cp","massBalanceError.txt",sim_title])
    call (["cp",layer + ".txt",sim_title])
    call (["cp","System.txt",sim_title])
    call (["cp","InjWells.txt",sim_title])
    call (["cp","p.csv",sim_title])
    call (["cp","scbar.csv",sim_title])
    call (["cp","MassBalance.csv",sim_title])
    call (["cp","thickness.txt",sim_title])
    print "moved files in "
    print clock() - t_read
   
