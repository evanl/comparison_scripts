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
    def __init__(self, vesa_simtitle, tough_simtitle,\
            sleipner, section,\
            vesa_hydro_folder = 'unif_hydro'):
        print "------------------------------------------------------------\n"
        print "Comparing "
        print "  VESA simulation: " + vesa_simtitle
        print "with"
        print "TOUGH2 simulation: " + tough_simtitle
        print "------------------------------------------------------------\n"
        current_directory = os.getcwd()
        vesa_dirname =  current_directory + "/vesa/VESA_02_01_13/" + vesa_simtitle
        tough_dirname = current_directory + "/tough/t2/" + tough_simtitle + "_dir/"

        os.chdir(tough_dirname)
        t2_grid, t2_timesteps = pt2.process_t2_output(tough_simtitle)
        # this reads the output for a given unit, but will not work for 
        # multiple units. the positions will have to be combined for that later
        os.chdir(vesa_dirname)
        vesa_cells, vesa_timesteps = vr.read_output_data(layer = vesa_simtitle)

        os.chdir(current_directory)
        print "yay"

        self.t2_grid = t2_grid
        self.t2_timesteps = t2_timesteps
        self.vesa_cells = vesa_cells
        self.vesa_timesteps = vesa_timesteps
        self.sleipner = sleipner
        self.section = section

    def create_blank_figure(self, compare_type = 'none'):
        if compare_type == 'section':
            self.f = plt.figure(figsize=(12.0,10.0), dpi=480)
        elif compare_type == 'horizontal':
            self.f = plt.figure(figsize=(10.0,12.0), dpi=480)
        else:
            self.f = plt.figure()

    def set_font_size(self, size = 'small'):
        #font = {'family': 'scalable', 'weight' : 'bold', 'size' : size}
        font = {'size' : size}
        matplotlib.rc('font', **font)

    def plot_contour(self, x, y, z, position = 111, label = 'values [units]',\
            xlab = 'x [m]', ylab = 'z [m]'):
        """ takes in numpy arrays of same shape
        """
        ax_c = self.f.add_subplot(position)
        if label == 'saturation []':
            n_levels = 8
            v = np.linspace(0., 0.8, num = n_levels)
            CS = ax_c.contourf(x,y,z,v)
        else:
            CS = ax_c.contourf(x,y,z)
        CB = plt.colorbar(CS, shrink = 1.0, pad=0.02, fraction = 0.07,\
                extend = 'both', format='%.2f')
        CB.set_label(label)
        ax_c.set_xlabel(xlab)
        ax_c.set_ylabel(ylab)
    def plot_graph(self, x, y, position = 111, \
            xlab = 'x [m]', ylab = 'z [m]'):
        """ takes in two 1d numpy arrays
        """
        ax_g = self.f.add_subplot(position)
        ax_g.plot(x,y )
        ax_g.set_xlabel(xlab)
        ax_g.set_ylabel(ylab)
    def add_t2_contours(self, axis, index, valtype):
        pos = 320
        for j in range(2):
            for i in range(2):
                self.t2_timesteps[i].make_plot_grid(self.t2_grid, axis, index,\
                        valtype)
                print "Self.sleipner" + str(self.sleipner)
                x, y, z = self.t2_timesteps[i].format_plot_grid(self.t2_grid, \
                        axis, self.sleipner, self.section)
                pos +=1
                self.plot_contour(x,y,z, position = pos, \
                        label = 'saturation []')
        return 0
    def add_vesa_sections(self, axis, index, nx):
        pos = 324
        for i in range(2):
            pos +=1
            ys, zb, zt, plume = vr.make_cross_sections(self.vesa_cells, \
                    i, axis, index, nx)
            self.plot_graph(ys, plume, pos)
            self.plot_graph(ys, zb, pos)
            self.plot_graph(ys, zt, pos)
        return 0
    def create_cross_section_comparison(self, title):
        self.create_blank_figure(compare_type ='section')
        fmt = 'eps'
        self.set_font_size(size = 8)
        nx = 65
        nx = 25
        x_ind = nx/2
        y_ind = 77
        self.add_t2_contours(2, x_ind, 'saturation')
        self.add_vesa_sections(1, nx/2, nx)
        #self.f.savefig('x_section_saturation.'+fmt, bbox_inches='tight',format=fmt)
        self.f.savefig(title + '.' + fmt, bbox_inches='tight',format=fmt)
        self.f.clf()


if __name__ == '__main__': 
    if len(sys.argv) != 3:
        sys.exit("Proper usage for compare.py is \n" +
                "python compare.py <vesa_simtitle> <tough_simtitle>\n")
    vesa_simtitle = sys.argv[1]
    tough_simtitle = sys.argv[2]
    sleipner = False
    section = False
    c = Compare(vesa_simtitle, tough_simtitle, sleipner, section)

    cl = [[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]
    z = np.asarray(cl)
    z = z / z.max()
    xl = [1.,2.,3.]
    yl = [4.,5.,6.]
    yl1 = [4.5, 5., 5.5]
    x, y = np.meshgrid(xl, yl)
    xg = np.asarray(xl)
    yg = np.asarray(yl1)

    title = 'unif_saturation_linear'
    c.create_cross_section_comparison(title)

