#Author - Evan Leister
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import string
from operator import itemgetter
import CO2Properties_1_0 as co2

class T2grid(object):
    """ This class represents an overarching structure to post-process a 
            TOUGH2 simulation. 
            The aim is that this will store all of the relevant information
            in dictionaries with a list of elements that can access them all.
    """
    def __init__(self):
        self.elements = []
        self.pres_init = self.get_initial_pressure()
        self.x = {}
        self.y = {}
        self.z = {}
        self.i = {}
        self.j = {}
        self.k = {}
        self.x_vals = []
        self.y_vals = []
        self.z_vals = []
        self.aq_mass_total = 0.
        self.aq_mass_water = 0.
        self.aq_mass_co2 = 0.
        self.aq_mass_nacl = 0.

    def get_initial_pressure(self):
        initial_pressure_dict = {}
        f = open('INCON','r')
        for i in range(2):
            line = f.readline()
        s = line.split()
        while s != []:
            if len(s[0]) != 5:
                s[0] = s[0] + " " + s[1]
            key = s[0]
            line = f.readline()
            s = line.split()
            initial_pressure_dict[key] = float(s[0])
            line = f.readline()
            s = line.split()
        return initial_pressure_dict

    def read_mesh(self):
        """ reads the MESH file in the directory and populates the 
                Elements list
                X,Y,Z coordinates
        """
        f = open('MESH','r')
        line = f.readline() #Should be ELEME
        line = f.readline()
        s = line.split()
        while  line != ('     \n') and s != []:
            if s[0] != 'ina':
                if len(s[0]) == 3:
                    eleme = s[0] + ' ' + s[1]
                elif len(s[0]) == 5 :
                    eleme = s[0]
                else :
                    print "error in read"
                self.elements.append(eleme)
                if len(s[-1]) == 30:
                    self.x[eleme] = float(s[-1][0:10])
                    self.y[eleme] = float(s[-1][10:20])
                    self.z[eleme] = float(s[-1][20:30])
                elif len(s[-1]) == 19:
                    self.x[eleme] = float(s[-2])
                    self.y[eleme] = float(s[-1][:9])
                    self.z[eleme] = float(s[-1][9:])
                else : 
                    self.x[eleme] = float(s[-3])
                    self.y[eleme] = float(s[-2])
                    self.z[eleme] = float(s[-1])

            line = f.readline()
            s = line.split()
        self.make_x_list()
        self.make_y_list()
        self.make_z_list()
        self.make_index_lists()
        return 0

    def make_x_list(self):
        print "making X list"
        for el in self.elements:
            flag = True
            if self.x_vals == []:
                self.x_vals.append(self.x[el])
            else:
                for entry in self.x_vals:
                    if entry == self.x[el]:
                        flag = False
                if flag:
                    self.x_vals.append(self.x[el])
        self.x_vals.sort()
        print "x has " + str(len(self.x_vals)) + " unique values"
        return 0

    def make_y_list(self):
        print "making Y list"
        for el in self.elements:
            flag = True
            if self.y_vals == []:
                self.y_vals.append(self.y[el])
            else:
                for entry in self.y_vals:
                    if entry == self.y[el]:
                        flag = False
                if flag:
                    self.y_vals.append(self.y[el])
        self.y_vals.sort()
        print "y has " + str(len(self.y_vals)) + " unique values"
        return 0

    def make_z_list(self):
        print "making Z list"
        for el in self.elements:
            flag = True
            if self.z_vals == []:
                self.z_vals.append(self.z[el])
            else:
                for entry in self.z_vals:
                    if entry == self.z[el]:
                        flag = False
                if flag:
                    self.z_vals.append(self.z[el]) 
        self.z_vals.sort()
        print "z has " + str(len(self.z_vals)) + " unique values"
        return 0

    def make_index_lists(self):
        i_map, j_map, k_map = get_index_maps()
        for el in self.elements:
            l = list(el)
            il = l[3] + l[4]
            jl = l[1] + l[2]
            kl = l[0]
            self.i[el] = i_map[il]
            self.j[el] = j_map[jl]
            self.k[el] = k_map[kl]
        return 0

    def read_initial_mass(self, f):
        print "Reading Initial Mass"
        line = f.readline()
        s = line.split()
        
        # read through lines until 'VOL.' is hit, 
        # get initial mass balances for each phase
        #  VOL. (m^3) * 0.00000000E+00 0.62971493E+07 0.00000000E+00
        while s == [] or s[0]!= 'VOL.':
            line = f.readline()
            s = line.split()
        self.gas_mass_water = float(s[9])
        self.gas_mass_nacl =  float(s[10])
        self.gas_mass_co2 =   float(s[11])
        
        #reads one line to get to masses and aqueous phase
        # MASS (kg)  * 0.00000000E+00 0.64182984E+10 0.00000000E+00 \
        # AQUEOUS    * 0.62100916E+10
        line = f.readline()
        s = line.split()
        self.gas_mass_total = float(s[3])
        self.aq_mass_total =float(s[4])
        self.aq_mass_water =float(s[8])
        self.aq_mass_nacl = float(s[9])
        self.aq_mass_co2 =  float(s[10])
        # Makes sure to not double read sand/system balances
        while s == [] or s[0]!= 'VOL.':
            line = f.readline()
            s = line.split()
        return "gas mass" + str(self.gas_mass_total )
        return f

class T2Timestep(object):
    """ 
    this class stores a grid for a given set of time output
    """
    def __init__(self):
        # injected mass
        self.step_time = 0. # days
        self.rate = 0. # kg/s
        self.mass_injected = 0. # kg

        # total mass inside system
        self.aq_mass_total = 0.
        self.aq_mass_water = 0.
        self.aq_mass_co2 = 0.
        self.aq_mass_nacl = 0.
        self.gas_mass_total = 0.
        self.gas_mass_water = 0. 
        self.gas_mass_co2 = 0.
        self.gas_mass_nacl = 0.

        # element-wise parameters
        self.elements = []
        self.pressure = {} 
        self.pres_diff = {}
        self.temperature = {}
        self.sat_gas = {}
        self.x_co2 = {}
        self.rho_gas = {}
        self.plot_grid = []

        # these two parameters are obtained from reading CO2TAB
        self.co2_viscosity = {}
        self.co2_density = {}
    
    def get_plot_value(self, element, valtype):
        if valtype == 'saturation':
            v = self.sat_gas[element]
        elif valtype == 'pressure':
            v = self.pressure[element]
        elif valtype == 'delta_p':
            v = self.pres_diff[element]
        elif valtype == 'temperature':
            v = self.temperature[element]
        elif valtype == 'internal_density':
            v = self.rho_gas[element]
        elif valtype == 'external_density':
            v = self.co2_viscosity[element]
        elif valtype == 'viscosity':
            v = self.co2_viscosity[element]
        else:
            print "please specify a valid element type: \n"
            print "saturation, pressure, delta_p, temperature"
            print "internal_density, external_density, viscosity"
            return 1
        return v

    def readOutput(self, f, g, grid, parallel = False, year_add = 0.,
            double_read = False):
        """
        reads through the already opened t2run.out file until it finds the 
        next block of output data. 
        Populates the pressure, saturation, density and time values for the 
        timestep
        """
        print "reading timestep"
        line = f.readline()
        s = line.split()

        # reads lines until the next OUTPUT block is encountered
        while s == [] or s[0] != 'OUTPUT':
            line = f.readline()
            s = line.split()

        # when the line with OUTPUT is encountered, the time is read
        # from that same line
        days_added = year_add * 365.25 
        self.step_time = float(s[-2]) + days_added
        print "step time, year_add"
        print self.step_time, year_add
        
        # reads through the block that looks like the following:
        # @@@@@@@@@@@@@@@ ...
        #                 ...
        #  TOTAL TIME     ...
        # 0.157790E+09    ...
        # @@@@@@@@@@@@@@@ ...
        #                 ...
        # in order to avoid the '@' symbols
        for i in range(6):
            line = f.readline()
            s = line.split()

        # reads blank lines or skips them until the end of the output block
        # signified by the line of '@' symbols
        while s == [] or line[1] != '@':
            # skips the lines that are blank
            if s == []:
                line = f.readline()
                s = line.split()
            # if a line contains the element header:
            # ELEM. INDEX   P      T       SG ...
            #             (Pa)   (deg.C)      ...
            # then this block skips it
            if s[0] == 'ELEM.':
                for i in range(3):
                    line = f.readline()
                s = line.split()


            # reads a subblock of output data that looks like:
            #  aJ3 1  25 0.83487E+07  32.00 ...
            #  aJ2 1  26 0.83584E+07  32.00 ...
            while s != [] and line [1] != '@' :
                # if a line contains the element header:
                # ELEM. INDEX   P      T       SG ...
                #             (Pa)   (deg.C)      ...
                # then this block skips it
                if s[0] == 'ELEM.':
                    for i in range(2):
                        line = f.readline()
                    s = line.split()


                # if the cell is not a dummy boundary cell
                # adds the element and its values to the dictionaries
                if s[0][:3] != 'ina' :

                    # handles the case for blank spaces 
                    # in between final characters
                    if len(s[0]) == 3:
                        s[0] = s[0] + ' ' + s[1]
                        del s[1]
                    # handles the case for full-length strings
                    # NOTE
                    # This will likely have to change once 
                    # the full 3d version of 
                    # tough2 is implemented and the indices change from 2 to 3d
                    if len(s[0]) == 9:
                        l = list(s[0])
                        el = ''.join(l[0:5])
                        ind = ''.join(l[5:9])
                        del s[0]
                        s.insert(0,ind)
                        s.insert(0,el)

                    eleme = s[0]
                    self.elements.append(eleme)
                    self.pressure[eleme] = float(s[2])
                    self.pres_diff[eleme] = \
                            self.pressure[eleme] - grid.pres_init[eleme]
                    self.temperature[eleme] = float(s[3])
                    self.sat_gas[eleme] = float(s[4])
                    self.x_co2[eleme] = float(s[8])
                    self.rho_gas[eleme] = float(s[11])
                line = f.readline()
                s = line.split()
            # end while for reading sub-block
        # end while for reading full output block
        # the current line is the '@@@@@@@@' line 
        # from the end of the output block

        if parallel == False:
            # for a single cell injection rate:
            while s == [] or s[0] != 'ELEMENT':
                line = f.readline()
                s = line.split()

            # once the line is element 
            for i in range(3):
                line = f.readline()
                s = line.split()
            if len(s[0]) == 3:
                self.rate = float(s[5])
            else: 
                self.rate = float(s[4])
            # calculates injected mass up until that point 
            # NOTE: 
            # FOR CONSTANT RATE ONLY. 
            self.mass_injected = self.rate * self.step_time * 24. * 3600.

            f = self.read_balance_block(f)
        else:
            g = self.read_balance_block(g, double_read = double_read)

        print "timestep read"
        return f

    def read_balance_block(self, f, double_read = False):
        line = f.readline()
        s = line.split()
        # read through lines until 'VOL.' is hit, get 
        # initial mass balances for each phase
        # VOL. (m^3) * 0.00000000E+00 0.62971493E+07 0.00000000E+00 \
        # GAS PHASE  * 0.16497908E+06 0.00000000E+00 0.11330593E+09 0.65583271E+14
        while s == [] or s[0]!= 'VOL.':
            line = f.readline()
            s = line.split()
        self.gas_mass_water = float(s[9])
        self.gas_mass_nacl =  float(s[10])
        self.gas_mass_co2 =   float(s[11])
        
        # PHASES     *      GAS          AQUEOUS        SOLID       \
        # COMPONENTS *     WATER          SALT             CO2           HEAT
        #reads one line to get to masses and aqueous phase
        # MASS (kg)  * 0.00000000E+00 0.64182984E+10 0.00000000E+00 \
        # AQUEOUS    * 0.62100916E+10
        line = f.readline()
        s = line.split()
        self.gas_mass_total = float(s[3])
        self.aq_mass_total =float(s[4])
        self.aq_mass_water =float(s[8])
        self.aq_mass_nacl = float(s[9])
        self.aq_mass_co2 =  float(s[10])
        # Makes sure to not double read sand/system balances
        if double_read == True:
            while s == [] or s[0]!= 'VOL.':
                line = f.readline()
                s = line.split()
        return f

    def get_co2_density_viscosity(self):
        NPK, NTK, PKV, TKV, Density, Viscosity, satLine = \
                co2.ReadTabFromFile()
        for el in self.elements:
            P = self.pressure[el]
            T = self.temperature[el]
            rho, mu = co2.GetCO2DensityViscosity(P, T, NPK, NTK, PKV, TKV, Density,\
                    Viscosity, satLine)
            self.co2_viscosity[el] = float(mu)
            self.co2_density[el] = float(rho)
        return 0 

    def make_plot_grid(self, grid, axis, index, valtype = False ):
        print "making plot grid"
        print "axis: " + str(axis)
        print "index: " + str(index)
        print "valtype: " + str(valtype)
        vals = []
        if axis == 1:
            for x in grid.x_vals:
                tempvals = []
                for el in self.elements:
                    if x == grid.x[el]:
                        if grid.j[el] == index:
                            v = self.get_plot_value(el, valtype)
                            tempvals.append(\
                                    (grid.x[el], grid.z[el], v))
                tempvals = sorted(tempvals, key = itemgetter(1))
                vals.append(tempvals)
        elif axis == 2:
            for y in grid.y_vals:
                tempvals = []
                for el in self.elements:
                    if y == grid.y[el]:
                        if grid.i[el] == index:
                            v = self.get_plot_value(el, valtype)
                            tempvals.append(\
                                    (grid.y[el], grid.z[el], v))
                tempvals = sorted(tempvals, key = itemgetter(1))
                vals.append(tempvals)
        elif axis == 3:                
            if valtype == 'thickness':
                vtype = 'saturation'
                nzmax = max(grid.k.iteritems(), key=itemgetter(1))[1]
                nzmin = min(grid.k.iteritems(), key=itemgetter(1))[1] 
                print "nzmin ", nzmin, "nzmax ", nzmax
                for i in range(len(grid.x_vals)):
                    tempvals = []
                    for j in range(len(grid.y_vals)):
                        num_cells = 0
                        sat_sum = 0.
                        for k in range(nzmin, nzmax):
                            el = get_element_chars(i, j, k)
                            v = self.get_plot_value(el, vtype)
                            num_cells +=1
                            sat_sum += v
                        sat_frac = sat_sum / num_cells
                        tempvals.append(\
                            (grid.x[el], grid.y[el], sat_frac))
                    tempvals = sorted(tempvals,key=itemgetter(1))
                    vals.append(tempvals)
            else:
                for x in grid.x_vals:
                    tempvals = []
                    for el in self.elements:
                        if x == grid.x[el]:
                            if grid.k[el] == index:
                                v = self.get_plot_value(el, valtype)
                                tempvals.append(\
                                    (grid.x[el], grid.y[el], v))
                    tempvals = sorted(tempvals,key=itemgetter(1))
                    vals.append(tempvals)
        else:
            print "make_plot_grid: Please Specify Valid Axis\n\n "
            return 1
        self.plot_grid = vals
        return 0

    def format_plot_grid(self, grid, axis, sleipner, section, shale = False):
        """ spits out 3 numpy arrays"""
        if axis == 1:
            nxp = len(grid.x_vals)
            if sleipner == True:
                if section == True:
                    nyp = 5
                else:
                    if shale == True:
                        nyp = 43
                    else:
                        nyp = 34
            else:
                nyp = len(grid.z_vals)
        elif axis == 2:
            nxp = len(grid.y_vals)
            if sleipner == True:
                if section == True:
                    nyp = 5
                else:
                    if shale == True:
                        nyp = 43
                    else:
                        nyp = 34
            else:
                nyp = len(grid.z_vals)
        elif axis == 3:
            nxp = len(grid.x_vals)
            nyp = len(grid.y_vals)

        xpl = np.zeros((nxp, nyp))
        ypl = np.zeros((nxp, nyp))
        val = np.zeros((nxp, nyp))

        if self.plot_grid != []:
            for i in range(nxp):
                for j in range(nyp):
                    #print i, j, len(self.plot_grid), len(self.plot_grid[0])
                    xpl[i][j] = float(self.plot_grid[i][j][0])
                    ypl[i][j] = float(self.plot_grid[i][j][1])
                    val[i][j] = float(self.plot_grid[i][j][2])
        else: 
            print "run T2grid.make_plot_grid() before using this function"
            return 1
        return xpl, ypl, val

    def plot_planar_timestep(self, grid, axis, index, valtype, \
            name = 'test',fmt='png', sleipner = True, \
            section = False, shale = True):
        valstr = "Valtype = " + str(valtype) + '\n'
        axstr = "Axis =  " + str(axis) + '\n'
        indstr = "Index = " + str(index) 
        print "Plotting Planar Timestep: " + '\n' + valstr + axstr + indstr
        xpl, ypl, val = self.format_plot_grid(grid, axis, sleipner, section, shale)
        font = { 'size' : '16'}
        matplotlib.rc('font', **font)
        f = plt.figure(num=None, dpi=480, facecolor= 'w',\
                #figsize=(7.5,10), 
            edgecolor ='k')
        ax = f.add_subplot(111)
        if sleipner == True:
            title_time = '{:.0f}'.format(self.step_time/ 365.25 + 1998)
            title_time = '{:.0f} days'.format(self.step_time)
        else:
            title_time = '{:.0f} days'.format(self.step_time)
        title_time = title_time.zfill(2)

        Nlevels = 21
        if valtype == 'pressure':
            CS = ax.contourf(xpl,ypl,val) 
            CB = plt.colorbar(CS, shrink = 0.8, extend = 'both')
            CB.set_label("Pressure [Pa]")
            #f.suptitle("T2slice, simulation: "+ name + ": + " + \
                    #"\n Pressure in  " + title_time + \
                    #axstr + indstr)
            f.suptitle("Time = " + title_time)
        elif valtype == 'saturation':
            V = np.linspace(0.,0.8,num=Nlevels)
            CS = ax.contourf(xpl,ypl,val,V) 
            CB = plt.colorbar(CS, shrink = 0.8, extend = 'both', ticks =V )
            CB.set_label("Saturation []")
            #f.suptitle("T2slice, simulation: "+ name + ": + " + \
                    #"\n Saturation in  " + title_time +\
                    #axstr + indstr)
            f.suptitle("Time = " + title_time)
        elif valtype == 'delta_p':
            CS = ax.contourf(xpl,ypl,val) 
            CB = plt.colorbar(CS, shrink = 0.8, extend = 'both')
            CB.set_label("Pressure Difference P - Po [Pa]")
            #f.suptitle("T2slice, simulation: "+ name + ": + " + \
                    #"\n Pressure Difference in  " + title_time+\
                    #axstr + indstr)
            f.suptitle("Time = " + title_time)
        else: 
            print "pressure or saturation ? "
            return 1

        if axis == 1:
            ax.set_xlabel('x-direction [m]')
            ax.set_ylabel('elevation [m]')
        elif axis == 2:
            ax.set_xlabel('y-direction [m]')
            ax.set_ylabel('elevation [m]')
        elif axis == 3:
            ax.set_xlabel('x-direction [m]')
            ax.set_ylabel('y-direction [m]')
            ax.set_aspect('equal')

        pltstring = '_'.join(['T2_slice', str(axis), str(index), name, valtype,\
                title_time ])
        plt.savefig(pltstring +\
                 '.'+ fmt, bbox_inches = 0, format = fmt)
        plt.close()
        plt.clf()
        return 0

def plot_mass_balance(grid, timestepList ):
    """ This function plots the mass balances of differenc components in the aqueous 
    phases and displays them in chart form.
    """
    matplotlib.rcParams.update({'font.size':12})
    
    # create total mass balance chart
    aq_m_tot = []
    aq_m_wat = []
    aq_m_co2 = []
    a1_m_nacl = []
    gas_m_tot = []
    gas_m_wat = []
    gas_m_co2 = []
    gas_m_nacl = []

    mass_injected = []
    total = []
    time_plot = []

    #for time t = 0
    aq_m_tot.append(grid.aq_mass_total)
    aq_m_wat.append(grid.aq_mass_water)
    aq_m_co2.append(grid.aq_mass_co2)
    a1_m_nacl.append(grid.aq_mass_nacl)

    gas_m_tot.append(0.)
    gas_m_wat.append(0.)
    gas_m_co2.append(0.)
    gas_m_nacl.append(0.)

    mass_injected.append(0.)
    total.append(grid.aq_mass_total)
    time_plot.append(0.)

    for step in timestepList:
        aq_m_tot.append(step.aq_mass_total)
        aq_m_wat.append(step.aq_mass_water)
        aq_m_co2.append(step.aq_mass_co2)
        a1_m_nacl.append(step.aq_mass_nacl)

        gas_m_tot.append(step.gas_mass_total)
        gas_m_wat.append(step.gas_mass_water)
        gas_m_co2.append(step.gas_mass_co2)
        gas_m_nacl.append(step.gas_mass_nacl)

        mass_injected.append(step.mass_injected)
        total.append(step.aq_mass_total + step.gas_mass_total)
        time_plot.append(step.step_time/(365 ))

    # create total phase curves
    aqtot = np.asarray(aq_m_tot) - aq_m_tot[0]
    gastot = np.asarray(gas_m_tot)
    inj = np.asarray(mass_injected)
    tot = np.asarray(total) - total[0]
    t = np.asarray(time_plot)
    f = plt.figure(num=None, figsize=(12,10), dpi=480, facecolor= 'w',\
        edgecolor ='k')
    # f.suptitle("Mass Components")
    ax = f.add_subplot(111)
    ax.plot(t,inj,     label = 'injected')
    ax.plot(t,aqtot,   label = 'delta aqueous')
    ax.plot(t,gastot,  label = 'gaseous')
    ax.plot(t,tot,     label = 'delta total')
    ax.set_xlabel('time [years]')
    ax.set_ylabel('phase mass [kg]')
    plt.legend(bbox_to_anchor = (0., 1.02, 1., 0.102), loc = 2,
            ncol= 4, mode="expand", borderaxespad=0.)
    plt.savefig('massbalance.png')
    plt.clf()
    plt.close()
    
    return 0

def plot_wellhead_pressure(grid, time_steps, \
        well_cell_id = 'bA1 1', fmt = 'png'):
    t_list = []
    p_list = []
    for step in time_steps:
        t_list.append(step.step_time)
        p_list.append(step.pressure[well_cell_id])

    # hydrostatic part
    f = open('INCON','r')
    line = f.readline()
    s = line.split()
    while s[0] != well_cell_id:
        line = f.readline()
        s = line.split()
        if len(s[0]) != 5:
            s[0] = s[0] + " " + s[1]
    line = f.readline()
    s = line.split()
    hydro_well_pressure = float(s[0])

    t_list.insert(0,0.)
    p_list.insert(0, hydro_well_pressure)
    time_ar = np.asarray(t_list)
    pres_ar = np.asarray(p_list) / 1000.
    f = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    f.suptitle('wellhead pressure vs time')
    ax = f.add_subplot(111)
    ax.set_xlabel('time[days]')
    ax.set_ylabel('pressure[kPa]')
    p = plt.plot(time_ar, pres_ar)
    f.savefig('wellhead_pressure' + '.' + fmt)
    plt.clf()
    plt.close()
    return 0

def plot_incon_change_vs_index(grid, time_steps, \
        show = False, fmt = 'png', vs_elev = False):
    pressure_difference =[]
    elev = []
    for i in range(len(grid.elements)):
        pressure_difference.append(time_steps[0].pres_diff[grid.elements[i]])
        elev.append(grid.z[grid.elements[i]])
    el_diff = np.asarray(elev)
    p_diff = np.asarray(pressure_difference)
    f = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    f.suptitle('pressure difference vs. elevation')
    ax = f.add_subplot(111)
    ax.set_ylabel('pressure difference [Pa]')
    if vs_elev == True:
        p = plt.scatter(el_diff, p_diff)
        ax.set_xlabel('elev. [m]')
    else:
        p = plt.plot(p_diff)
        ax.set_xlabel('index []')
    f.savefig('p_diff_vs_index' + '.' + fmt)
    plt.clf()
    if show == True:
        plt.show()
    plt.close()
    return 0

def plot_planar_contours(grid, time_steps, sim_title, fmt='png',\
        two_d = True, sleipner = True, section = False, shale = True,\
        axis = 3,\
        i_in = False, j_in = False, k_in = False):
    """ if axis is 3, standard 2d contour plots will be generated with 
        k as the z-direction index
        if axis is 1, a cross section will be made along the x axis
        with j-index
        if axis is 2, a cross section will be made along the x axis with 
        i index.
    """
    for i in range(len(time_steps)):
        #print "Plotting timestep " + str(i)
        #print "time = " + str(float(time_steps[i].step_time)/365.) + " years" 

        if axis == 1:
            index = j_in
        elif axis == 2:
            index = i_in
        elif axis == 3:
            index = k_in

        if two_d == True:
            index = 0
            axis = 3
        # create pressure difference 
        time_steps[i].make_plot_grid(grid, axis = axis, index = index, \
                valtype = 'delta_p')
        time_steps[i].plot_planar_timestep(grid, axis = axis, index = index,\
                valtype = 'delta_p',\
                name = sim_title, fmt = fmt, sleipner = sleipner,\
                section = section, shale = shale)
        # create pressure contours
        time_steps[i].make_plot_grid(grid, axis = axis, index = index, \
                valtype = 'pressure')
        time_steps[i].plot_planar_timestep(grid, axis = axis, index = index, \
                valtype = 'pressure',\
                name = sim_title, fmt = fmt, sleipner = sleipner,\
                section = section, shale = shale)
        # create saturation
        time_steps[i].make_plot_grid(grid, axis = axis, index = index, \
                valtype = 'saturation')
        time_steps[i].plot_planar_timestep(grid, axis = axis, index = index,\
                valtype = 'saturation',\
                name = sim_title, fmt = fmt, sleipner = sleipner,\
                section = section, shale = shale)
    return 0

def check_3d_hydro_pressure(grid, time_steps):
    print "Checking Hydrostatic Equilibrium"
    print len(time_steps)
    for i in range(len(time_steps)):
        if i != 0:
            count = 0
            for el in grid.elements:
                p1 = time_steps[i].pressure[el]
                p2 = time_steps[i-1].pressure[el]
                x = grid.x[el]
                y = grid.y[el]
                z = grid.z[el]
                if p1 != p2:
                    count +=1
                    ## printout used for checkijk
                    #print i, el, p1, p2, abs(p1 - p2), x, y, z
            print "step " + str(i)+ " compared to step " + str(i-1)
            print "number of inconsistent elements" + str(count)
    return 0

def get_x_char( i):
    if i > 99 or i < 0 :
        print "change your x indexing system"
        return 1
    else : 
        c1 = str(i)
        if len(c1) == 1:
            char = ' ' + c1
        elif len(c1) == 2:
            char = c1
        else:
            print "I'm pretty sure you probably put an index less than 1"
        return char

def get_y_char(j):
    if j > 259 or j < 0:
        print "change y indexing system"
        return 1
    else : 
        letters = string.uppercase
        c1 = letters[ j/10]
        c2 = str( j - 10 * (j/10))
        char = c1 + c2
        return char

def get_z_char(k):
    if k > 51 or k < 0:
        print 'fix your z indecisis'
        return 1
    else:
        letters =  string.lowercase + string.uppercase
        char = letters[k]
        return char

def get_element_chars(i, j, k):
    c1 = get_x_char( i)
    c2 = get_y_char( j)
    c3 = get_z_char( k)
    el = c3 + c2 + c1
    if len(el) != 5:
        print "i, j, k"
        print i, j, k 
        print "c1, c2, c3"
        print c1, c2, c3
        print "failed to give correct character length"
        return 1
    else:
        return el

def get_index_maps():
    # populate list of i j k element chars
    i_map = {}
    j_map = {}
    k_map = {}
    for i in range(90):
        c = get_x_char(i)
        i_map[c] = i
    for j in range(150):
        c = get_y_char(j)
        j_map[c] = j
    for k in range(45):
        c = get_z_char(k)
        k_map[c] = k
    return i_map, j_map, k_map

def write_viscosity(grid, time_steps):
    f= open('rho_visc','w')
    for i in range(len(time_steps)):
        stepno = '{:d}'.format(i)
        f.write('output step ' + stepno + '\n')
        count = 0
        internal_count = 0
        internal_rho = 0.
        sum_mu = 0.
        sum_rho = 0.
        for el in time_steps[i].elements:
            x = '{:.0f}'.format(grid.x[el])
            y = '{:.0f}'.format(grid.y[el])
            z = '{:.0f}'.format(grid.z[el])
            mu = '{:.0f}'.format(time_steps[i].co2_viscosity[el])
            #f.write(el + x + y + z + mu + '\n')
            count += 1
            sum_mu += time_steps[i].co2_viscosity[el]
            sum_rho += time_steps[i].co2_density[el]
            rho_gas = time_steps[i].rho_gas[el]
            if rho_gas > 10.:
                internal_count += 1
                internal_rho += rho_gas
        f.write('Average CamelotViscosity for time_step ' + stepno + ": " + \
                str(sum_mu / float(count)) + '\n')
        f.write('Average CamelotDensity for time_step ' + stepno + ": " + \
                str(sum_rho / float(count)) + '\n')
        if internal_count != 0:
            f.write('Average OUTPUT Density for time_step ' + stepno + ": " + \
                    str(internal_rho / float(internal_count)) + '\n')
        else: 
            f.write('NO CO2 in DOMAIN for this simulation')
    f.close()
    return 0

def read_coft():
    try:
        with open('COFT','r'):
            f = open('COFT','r')
            line = f.readline()
            s = line.split(',')
            time = []
            flux_1 = []
            flux_2 = []
            flux_3 = []
            while line:
                time.append(float(s[1]))
                flux_1.append(float(s[3]))
                flux_2.append(float(s[4]))
                flux_3.append(float(s[5]))
                line = f.readline()
                s = line.split(',')

            t = np.asarray(time)
            f1 = np.asarray(flux_1)
            f2 = np.asarray(flux_2)
            f3 = np.asarray(flux_3)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(t, f1, label='gas flow')
            ax.plot(t, f2, label='liquid flow')
            ax.plot(t, f3, label='total flow')
            ax.legend()
            plt.show()
    except IOError:
        print "COFT WAS NOT GENERATED"
    return 0

def plot_timesteps(filename):
    f = open(filename,'r')
    line = f.readline()
    s = line.split()
    time = []
    timesteps = []
    while line:
        if len(s) > 8 and s[6] == 'DT':
            time.append(float(s[5])/(3600.))
            timesteps.append(float(s[8])/3600.)
        elif len(s) > 8 and s[5] == 'DT':
            time.append(float(s[4])/(3600.))
            timesteps.append(float(s[7])/3600.)
        line = f.readline()
        s = line.split()
    f.close()
    steps = np.asarray(timesteps)
    time_vals = np.asarray(time)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time_vals, steps)
    ax.set_xlabel('time[hr]')
    ax.set_ylabel('timestep[hr]')
    plt.show()
    
