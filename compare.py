#Author - Evan Leister
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from time import time, clock
import sys
import os
from subprocess import call
import tough.t2.process_t2_output as pt2
import vesa.VESA_02_01_13.vesa_reading_functions as vr

class Compare:
    def __init__(self, sim_titles, num_vesa_sims, 
            sleipner, section,\
            vesa_hydro_folder = 'unif_hydro'):
        self.num_vesa_sims = num_vesa_sims
        self.simulations = []
        print "------------------------------------------------------------\n"
        print "Comparing VESA Simulations:"
        for i in range(num_vesa_sims):
            print 20*" " + sim_titles[i]
        print "with"
        print "          TOUGH2 Simulations:"
        for i in range(num_vesa_sims, len(sim_titles)):
            print  20*" "+ tough_simtitle
        print "------------------------------------------------------------\n"
        current_directory = os.getcwd()
        for i in range(num_vesa_sims):
            print "reading VESA Simulation: " + sim_titles[i]
            vesa_dirname =  current_directory + "/vesa/VESA_02_01_13/" + sim_titles[i]
            os.chdir(vesa_dirname)
            vesa_cells, vesa_timesteps = vr.read_output_data(layer = vesa_simtitle)
            self.simulations.append(Simulation(sim_titles[i], 'vesa', \
                    vesa_cells, vesa_timesteps))
        for i in range(num_vesa_sims, len(sim_titles)):
            print "reading TOUGH2 Simulation: " + sim_titles[i]
            tough_dirname = current_directory + "/tough/t2/" + sim_titles[i] + "_dir/"
            os.chdir(tough_dirname)
            t2_grid, t2_timesteps = pt2.process_t2_output(tough_simtitle)
            self.simulations.append(Simulation(sim_titles[i], 'tough', \
                    t2_grid, t2_timesteps))
        # this reads the output for a given unit, but will not work for 
        # multiple units. the positions will have to be combined for that later

        os.chdir(current_directory)
        print "YAY"

        self.sleipner = sleipner
        self.section = section

    def create_blank_figure(self, fontsize = 10, compare_type = 'none'):
        self.set_font_size(size = fontsize)
        if compare_type == 'section':
            self.fig = plt.figure(figsize=(12.0,10.0), dpi=480)
        elif compare_type == 'horizontal':
            self.fig = plt.figure(figsize=(8.0,10.0), dpi=480)
        else:
            self.fig = plt.figure()
        return 0

    def set_font_size(self, size = 'small'):
        #font = {'family': 'scalable', 'weight' : 'bold', 'size' : size}
        font = {'size' : size}
        matplotlib.rc('font', **font)
        return 0

    def plot_contour(self, x, y, z, position = 111, label = False,\
            title = 'saturation', xlab = 'x [m]', ylab = 'z [m]'):
        """ takes in numpy arrays of same shape
        """
        ax_c = self.fig.add_subplot(position)
        if label == 'saturation []':
            n_levels = 8
            v = np.linspace(0., 0.8, num = n_levels)
            CS = ax_c.contourf(x,y,z,v)
        else:
            CS = ax_c.contourf(x,y,z)
        cbaxes = self.fig.add_axes([0.8,0.1,0.03,0.8])
        CB = plt.colorbar(CS, shrink = 1.0, pad=0.02, fraction = 0.07,\
                extend = 'both', format='%.2f')
        #if label != False:
            #CB.set_label(label)

        ax_c.set_xlabel(xlab)
        ax_c.set_ylabel(ylab)
        ax_c.set_title(title)
        return 0

    def plot_graph(self, x, y, position = 111, \
            title = 'saturation', xlab = 'x [m]', ylab = 'z [m]'):
        """ takes in two 1d numpy arrays
        """
        ax_g = self.fig.add_subplot(position)
        ax_g.plot(x,y )
        ax_g.set_xlabel(xlab)
        ax_g.set_ylabel(ylab)
        ax_g.set_title(title)
        return 0

    def add_t2_contours(self, sim_index, axis, space_index, valtype,\
            time_indices = [0,1]):
        pos = 220
        for time_index in time_indices:
            self.simulations[sim_index].time_steps[time_index].make_plot_grid(\
                    self.simulations[sim_index].grid_cells, \
                    axis, space_index,\
                    valtype)
            x, y, z = self.simulations[sim_index].time_steps[time_index].format_plot_grid(\
                    self.simulations[sim_index].grid_cells, \
                    axis, self.sleipner, self.section)
            pos +=1
            label_str = self.year_index(time_index)
            if axis == 3:
                xlab = 'x [m]'
                ylab = 'y [m]'
            elif axis == 2:
                xlab = 'y [m]'
                ylab = 'z [m]'
            elif axis == 1:
                xlab = 'x [m]'
                ylab = 'z [m]'
            self.plot_contour(x,y,z, position = pos, \
                    xlab = xlab, ylab = ylab,\
                    title = label_str, label = 'saturation []')
        return 0

    def add_vesa_sections(self, sim_index, axis, space_index, nx,\
            time_indices=[0,1]):
        pos = 222
        for time_index in time_indices:
            pos +=1
            ys, zb, zt, plume = vr.make_cross_sections(\
                    self.simulations[sim_index].grid_cells, \
                    time_index, axis, space_index, nx)
            title_string = self.year_index(time_index)
            self.plot_graph(ys, plume, pos)
            self.plot_graph(ys, zb, pos)
            self.plot_graph(ys, zt, pos, title = title_string)
        return 0

    def add_vesa_contours(self, sim_index, nx, ny, valtype, \
            time_indices=[0,1]):
        pos = 222
        for time_index in time_indices:
            pos +=1
            x, y, zval = vr.make_plot_grid(\
                    self.simulations[sim_index].grid_cells, time_index, \
                    nx, ny, valtype)
            self.plot_contour(x, y, zval, position = pos,\
                    xlab = 'x [m]', ylab = 'y [m]',\
                    label = valtype)
        return 0

    def create_cross_section_comparison(self, title, time_indices=[0,1]):
        self.create_blank_figure(compare_type ='section')
        fmt = 'png'
        nx = 65
        ny = 119
        x_ind = nx/2
        y_ind = ny/2
        v_index = 0
        t_index = 1
        self.add_t2_contours(t_index, 2, x_ind, 'saturation',\
                time_indices = time_indices)
        self.add_vesa_sections(v_index, 2, x_ind, nx,\
                time_indices = time_indices)
        self.fig.tight_layout()
        self.fig.savefig(title + '.' + fmt, bbox_inches='tight',format=fmt)
        self.fig.clf()
        return 0

    def create_contour_comparison(self, title, time_indices=[0,1]):
        nx = 65
        ny = 119
        z_ind = 3
        self.create_blank_figure(fontsize = 10, compare_type = 'horizontal')
        fmt = 'png'
        valtype = 'saturation'
        v_index = 0
        t_index = 1
        self.add_t2_contours(t_index, 3, z_ind, valtype,\
                time_indices = time_indices)
        self.add_vesa_contours(v_index,nx, ny, valtype,\
                time_indices = time_indices)
        self.fig.tight_layout()
        self.fig.savefig(title + '.' + fmt, bbox_inches='tight', format=fmt)
        self.fig.clf()
        return 0

    def year_index(self, time_index):
        yrstring = '{:4d}'.format(int(1999 + time_index))
        return yrstring

class Simulation:
    def __init__(self, sim_title, sim_type, grid_cells, time_steps):
        self.sim_title = sim_title
        self.sim_type = sim_type
        self.grid_cells = grid_cells
        self.time_steps = time_steps


if __name__ == '__main__': 
    if len(sys.argv) != 3:
        sys.exit("Proper usage for compare.py is \n" +
                "python compare.py <vesa_simtitle> <tough_simtitle>\n")
    sim_titles = []
    for i in range(1,len(sys.argv)):
        sim_titles.append(sys.argv[i])
    vesa_simtitle = sys.argv[1]
    tough_simtitle = sys.argv[2]
    num_vesa_sims = 1
    sleipner = True
    section = False
    c = Compare( sim_titles, num_vesa_sims, sleipner, section)
    temp = '32'
    yr_index = [0,2,4,6,8]
    for i in yr_index:
        time_indices = [i, i+1]
        if i == 0:
            yri = 99
        else:
            yri = i-1
        title = 'sleipner_' + temp + '_section_' + str(yr_index.index(i)) + \
                '_' + '{:2d}'.format(yri) + '_' + \
                '{:2d}'.format(i)
        c.create_cross_section_comparison(title, time_indices = time_indices)
        contour_title = 'sleipner_' + temp + '_contour_' +\
                str(yr_index.index(i)) + \
                '_' + '{:2d}'.format(i-1) + '_' + \
                '{:2d}'.format(i)
        c.create_contour_comparison(contour_title, time_indices = time_indices)
