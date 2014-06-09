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
            sleipner, section, split = 0,\
            vesa_hydro_folder = 'unif_hydro', parallel = False, \
            double_balance = False,\
            t2_read = True, vesa_read = True):
        self.num_vesa_sims = num_vesa_sims
        self.simulations = []
        self.part2 = parallel
        self.split = split
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
        if vesa_read == True:
            for i in range(num_vesa_sims):
                print "reading VESA Simulation: " + sim_titles[i]
                vesa_dirname =  current_directory + "/vesa/VESA_02_01_13/" + \
                        sim_titles[i]
                os.chdir(vesa_dirname)
                vesa_cells, vesa_timesteps = vr.read_output_data(layer = \
                        vesa_simtitle)
                self.simulations.append(Simulation(sim_titles[i], 'vesa', \
                        vesa_cells, vesa_timesteps))
        if t2_read == True:
            for i in range(num_vesa_sims, len(sim_titles)):
                print "reading TOUGH2 Simulation: " + sim_titles[i]
                tough_dirname = current_directory + "/tough/t2/" + sim_titles[i] + "_dir/"
                os.chdir(tough_dirname)
                t2_grid, t2_timesteps = pt2.process_t2_output(tough_simtitle, \
                        parallel = self.part2, split = self.split,\
                        double_balance = double_balance)
                self.simulations.append(Simulation(sim_titles[i], 'tough', \
                        t2_grid, t2_timesteps))
        # this reads the output for a given unit, but will not work for 
        # multiple units. the positions will have to be combined for that later

        os.chdir(current_directory)
        print "YAY"

        self.sleipner = sleipner
        self.section = section

    def create_blank_figure(self, fontsize = 14, compare_type = 'none'):
        self.set_font_size(size = fontsize)
        if compare_type == 'section':
            self.fig = plt.figure(figsize=(20.0,6.0), dpi=960)
        elif compare_type == 'horizontal':
            self.fig = plt.figure(figsize=(8.0,10.0), dpi=960)
        elif compare_type == 'plume':
            self.fig = plt.figure(figsize=(10.0,2.5), dpi=960)
        else:
            self.fig = plt.figure()
        return 0

    def set_font_size(self, size = 'small'):
        #font = {'family': 'scalable', 'weight' : 'bold', 'size' : size}
        font = {'size' : size}
        matplotlib.rc('font', **font)
        return 0

    def plot_contour(self, x, y, z, position = 111, ax_label = False,\
            contour_label = True, label = 'saturation []',\
            title = 'saturation', xlab = 'x [m]', ylab = 'z [m]',\
            section = False):
        """ takes in numpy arrays of same shape
        """
        ax_c = self.fig.add_subplot(position)
        plt.tick_params(which='major', length=3, color = 'w')
        if label == 'saturation []':
            n_levels = 8
            v = np.linspace(0., 0.8, num = n_levels)
            CS = ax_c.contourf(x,y,z,v)
        else:
            CS = ax_c.contourf(x,y,z)
        if contour_label == True:
            self.fig.subplots_adjust(right=0.84)
            cb_axes = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
            plt.tick_params(which='major', length=3, color = 'k')
            self.fig.colorbar(CS, cax = cb_axes, format = '%.2f')
            #CB = plt.colorbar(CS, shrink = 1.0, pad=0.02, fraction = 0.07,\
                    #extend = 'both', format='%.2f')
        if ax_label != True:
            ax_c.set_yticklabels([])
        if section == True:
            ax_c.xaxis.set_ticks(np.arange(0,6000,2500))
        else:
            ax_c.xaxis.set_ticks(np.arange(0,3500,1000))
            ax_c.set_xticklabels([])


        #ax_c.set_xlabel(xlab)
        #ax_c.set_aspect('equal')
        ax_c.set_title(title)
        return 0

    def plot_graph(self, x, y, position = 111, \
            title = 'saturation', xlab = 'x [m]', ylab = 'z [m]',\
            legend_label = 'y', legend = False,\
            ax_label = True):
        """ takes in two 1d numpy arrays
        """
        ax_g = self.fig.add_subplot(position)
        ax_g.plot(x,y, label = legend_label)
        ax_g.set_xlabel(xlab)
        ax_g.set_ylabel(ylab)
        ax_g.set_title(title)
        ax_g.xaxis.set_ticks(np.arange(0,6000,2500))
        if ax_label != True:
            ax_g.set_yticklabels([])
        #if legend == True:
            #leg_axes = self.fig.add_axes([0.85,0.15,0.05,0.7])
            #leg = plt.legend(loc=1)
            ##ax_g.legend(loc=1)
        return 0

    def add_t2_contours(self, sim_index, axis, space_index, valtype,\
            time_indices = [0,1], section = False):
        pos = 130
        for time_index in time_indices:
            contour_label = False
            ax_label = False
            if time_indices.index(time_index) == 0:
                ax_label = True
            elif time_indices.index(time_index) == 2:
                contour_label = True
            label_str = self.year_index(time_index)
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
                    contour_label = contour_label, ax_label = ax_label,\
                    xlab = xlab, ylab = ylab,\
                    title = label_str, label = 'saturation []', 
                    section = section)
        return 0

    def add_vesa_sections(self, sim_index, axis, space_index, nx,\
            time_indices=[0,1]):
        pos = 130
        for time_index in time_indices:
            legend = False
            ax_label = False
            if time_indices.index(time_index) == 0:
                ax_label = True
            elif time_indices.index(time_index) == 2:
                legend = True
            pos +=1
            ys, zb, zt, plume = vr.make_cross_sections(\
                    self.simulations[sim_index].grid_cells, \
                    time_index, axis, space_index, nx)
            title_string = self.year_index(time_index)
            self.plot_graph(ys, plume, pos, legend_label='CO2 Thickness')
            self.plot_graph(ys, zb, pos, legend_label='Bottom Shale')
            self.plot_graph(ys, zt, pos, title = title_string, \
                    ax_label = ax_label,\
                    legend_label='Caprock', legend = legend)
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

    def create_cross_section_comparison(self, title, sec_type = 'tough',\
            time_indices=[0,1], fmt = 'png'):
        self.create_blank_figure(fontsize = 16, compare_type ='section')
        nx = 65
        ny = 119
        x_ind = nx/2
        y_ind = ny/2
        v_index = 0
        t_index = 1
        valtype = 'saturation'
        if sec_type == 'tough':
            self.add_t2_contours(t_index, 2, x_ind, valtype,\
                    time_indices = time_indices, section = True)
        elif sec_type == 'vesa':
            self.add_vesa_sections(v_index, 1, x_ind, nx,\
                    time_indices = time_indices)
        #self.fig.tight_layout()
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

    def create_vesa_plume_match(self, title, fmt):
        pos = 160
        nx = 65
        ny = 119
        self.create_blank_figure(fontsize = 14, compare_type = 'plume')
        time_indices = [0, 2, 3, 5, 7, 9]
        sim_index = 0
        for time_index in time_indices:
            pos +=1
            contour_label = False
            ax_label = False
            if time_indices.index(time_index) == 0:
                ax_label = True
            elif time_indices.index(time_index) == 5:
                contour_label = True
            label_str = self.year_index(time_index)
            x, y, zval = vr.make_plot_grid(\
                    self.simulations[sim_index].grid_cells, time_index, \
                    nx, ny, 'saturation')
            label_str = self.year_index(time_index)
            self.plot_contour(x, y, zval, position = pos,\
                    xlab = 'x [m]', ylab = 'y [m]',\
                    ax_label = ax_label, contour_label = contour_label, \
                    title = label_str, label = 'saturation []')
        #self.fig.tight_layout()
        self.fig.savefig(title + '.' + fmt, bbox_inches='tight',format=fmt)
        self.fig.clf()
        return 0

    def create_tough_plume_match(self, title, fmt, thickness = True):
        pos = 160
        self.create_blank_figure(fontsize = 14, compare_type = 'plume')
        time_indices = [0, 2, 3, 5, 7, 9]
        sim_index = 1
        axis = 3
        space_index = 3
        if thickness == True:
            thick = 'thickness'
        else:
            thick = 'saturation'
        for time_index in time_indices:
            pos +=1
            contour_label = False
            ax_label = False
            if time_indices.index(time_index) == 0:
                ax_label = True
            elif time_indices.index(time_index) == 5:
                contour_label = True
            label_str = self.year_index(time_index)
            self.simulations[sim_index].time_steps[time_index].make_plot_grid(\
                    self.simulations[sim_index].grid_cells, \
                    axis, space_index, thick)
            x, y, z = self.simulations[sim_index].time_steps[time_index].format_plot_grid(\
                    self.simulations[sim_index].grid_cells, \
                    axis, self.sleipner, self.section)
            label_str = self.year_index(time_index)
            self.plot_contour(x,y,z, position = pos, \
                    xlab = 'x [m]', ylab = 'y [m]',\
                    ax_label = ax_label, contour_label = contour_label, \
                    title = label_str, label = 'saturation []')
        #self.fig.tight_layout()
        self.fig.savefig(title + '.' + fmt, bbox_inches='tight',format=fmt)
        self.fig.clf()
        return 0

    def year_index(self, time_index):
        yrstring = '{:4d}'.format(int(1999 + time_index))
        return yrstring

    def make_2yr_contours_sections(self, temp):
        yr_index = [0,2,4,6,8]
        for i in yr_index:
            time_indices = [i, i+1]
            if i == 0:
                yri = 99
            else:
                yri = i-1
            title = 'sleipner_' + temp + '_section_' + \
                    str(yr_index.index(i)) + \
                    '_' + '{:2d}'.format(yri) + '_' + \
                    '{:2d}'.format(i)
            self.create_cross_section_comparison(title, \
                    time_indices = time_indices)
            contour_title = 'sleipner_' + temp + '_contour_' +\
                    str(yr_index.index(i)) + \
                    '_' + '{:2d}'.format(i-1) + '_' + \
                    '{:2d}'.format(i)
            self.create_contour_comparison(contour_title, \
                    time_indices = time_indices)
        return 0

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
    parallelt2 = True
    double_balance = True
    split = 0
    fmt = 'pdf'
    t2_read = True
    vesa_read = True
    c = Compare( sim_titles, num_vesa_sims, sleipner, section, \
            parallel = parallelt2, split = split,\
            double_balance = double_balance, t2_read = t2_read,\
            vesa_read = vesa_read)

    #TOUGH2
    if t2_read == True:
        toughid = '42_nocap'
        sec_type = 'tough'
        title ='t' + toughid + '_' + sec_type + '_section' 
        c.create_cross_section_comparison(title, sec_type = sec_type,\
                time_indices = [1, 4, 7], fmt = fmt)
        tpl = 't' + toughid + '_tough_plume'
        c.create_tough_plume_match(tpl, fmt, thickness = True)

    #VESA
    if vesa_read == True:
        vesa_id = '42'
        sec_type = 'vesa'
        title ='v' + vesa_id + '_' + sec_type + '_section' 
        c.create_cross_section_comparison(title, sec_type = sec_type,\
                time_indices = [1, 4, 7], fmt = fmt)
        vpl = 'v' + vesa_id + '_vesa_plume_sharp'
        c.create_vesa_plume_match(vpl, fmt)

