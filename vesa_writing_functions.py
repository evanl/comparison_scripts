#Author - Evan Leister
import eclipse_cells as ec
import numpy as np
import matplotlib.pyplot as plt

class Injector(object):
    def __init__(self, index, x, y, ratio, layer_id, end_days, mass_rate):
        self.index = index
        self.x = x
        self.y = y
        self.ratio = ratio
        self.layer_id = layer_id
        self.end_days = end_days # list of intervals in integer days
        self.mass_rate = mass_rate # list of mass rate is given in Mt/yr
        self.radius = 1.0
        if len(mass_rate) != len(end_days):
            print "Mass inflow and interval ends must match"
            return 1
    def write_injector(self, f):
        print "mass rate is"
        print self.mass_rate
        f.write(', '.join([str(self.index), str(self.x), str(self.y),\
                str(self.ratio), str(self.layer_id),\
                str(self.radius), str(len(self.mass_rate))]))
        f.write('\n')
        for i in range(len(self.mass_rate)):
            print len(self.mass_rate), len(self.end_days)
            f.write(', '.join([str(self.end_days[i]), \
                    str(self.mass_rate[i])]))
            f.write('\n')
        return f

def write_injwells(injectors):
    f = open("InjWells.txt","w")
    print "writing Injwells.txt"
    for inj in injectors:
        inj.write_injector(f)
    f.close()
    return 0

def write_system(timestep_days, output_days, simtime_years, output_control, layers):
    f = open("System.txt", "w")
    print "writing System.txt"
    f.write(''.join([str(timestep_days),'\n']))
    f.write(''.join([str(output_days),'\n']))
    years = simtime_years
    f.write(''.join([str(simtime_years),'\n']))
    f.write(''.join([output_control,'\n']))
    # mass balance output 'string' massbalance nomassbalance
    f.write('massbalance \n')
    # number of layers in the model [int]
    f.write(''.join([str(len(layers)),'\n']))
    # file names that contain the grid data
    for s in layers:
        f.write(s)
    f.close()
    return 0

class Layer(object):
    def __init__(self, layer_name, l_type, l_id, l_co2_rho, l_bri_rho, l_co2_mu,\
            l_bri_mu, sc_res, sb_res, c_co2, c_bri, c_roc, cap_rp_id, \
            nx, ny, nz = 1, gradient = 10.,\
            homogeneous = False, permval = 2000., poroval = False):
        """
        self.layer_name : Name of Text File
        l_type : Layer Type
        l_id : Layer ID
        l_co2_rho : CO2 density
        self.l_bri_rho : Brine Density
        self.l_co2_mu : CO2 viscosity
        self.l_bri_mu : brine viscosity
        self.sc_res : CO2 residual saturation
        self.sb_res : brine residual saturation
        self.c_co2 : CO2 compressibility
        self.c_bri : Brine Compressibility
        self.c_roc : Rock Compressibility
        self.cap_rp_id : Capillary - Rel/perm ID
        self.gradient : Pressure gradient
        """
        self.layer_name = layer_name
        self.l_type = l_type
        self.l_id = l_id
        self.l_co2_rho = l_co2_rho
        self.l_bri_rho = l_bri_rho
        self.l_co2_mu = l_co2_mu
        self.l_bri_mu = l_bri_mu
        self.sc_res = sc_res
        self.sb_res = sb_res
        self.c_co2 = c_co2
        self.c_bri = c_bri
        self.c_roc = c_roc
        self.cap_rp_id = cap_rp_id
        # list of GridCell objects
        self.grid_cells = []
        self.gradient = gradient #[MPa/km]
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.homogeneous = homogeneous
        self.permval = permval
        self.poroval = poroval


    def fill_uniform_grid(self, dx, dy, dz, center_depth, phi, k):
        for j in range(self.ny):
            for i in range(self.nx):
                x = (dx/2. + dx * i)
                y = (dy/2. + dy * j)
                top_b = -center_depth + dz/2.
                bottom_b = -center_depth - dz/2.

                if i == (self.nx - 1):
                    east_bc = 3
                else:
                    east_bc = 1

                if j == (self.ny - 1):
                    north_bc = 3
                else:
                    north_bc = 1

                if i == 0:
                    west_bc = 3
                else: 
                    west_bc = 1
                    
                if j == 0:
                    south_bc = 3
                else:
                    south_bc = 1

                pressure = -self.gradient * bottom_b * 1000
                gc = GridCell(top_b, bottom_b, x, y, dx, dy, phi, k,\
                        west_bc, east_bc, south_bc, north_bc, pressure)
                self.grid_cells.append(gc)
            

        return 0

    def plot_perm_data(self, e_cells):
        cell_ind = np.zeros(len(e_cells))
        anis = np.zeros(len(e_cells))
        depth = np.zeros(len(e_cells))
        perm = np.zeros(len(e_cells))
        poro = np.zeros(len(e_cells))
        for i in range(len(e_cells)):
            cell_ind[i] = i
            perm[i] = e_cells[i].getXPermeability()
            anis[i] = e_cells[i].getZPermeability() / \
                        e_cells[i].getXPermeability()
            depth[i] = e_cells[i].getTopZ()
            poro[i] = e_cells[i].getPorosity()
        print "plotting anisotropy ratio"
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        a = ax1.plot(cell_ind, anis)
        ax1.set_xlabel('cell_index []')
        ax1.set_ylabel('anisotropy kz/kx')
        plt.savefig('ec_anis_cells.png')
        plt.close()
        print "plotting permeabilityvsdepth"
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        b = ax2.scatter(depth, perm)
        ax2.set_xlabel('depth [m]')
        ax2.set_ylabel('permeability [md]')
        plt.savefig('ec_perm_depth.png')
        plt.close()
        print "plotting porosity vsdepth"
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        b = ax3.scatter(depth, poro)
        ax3.set_xlabel('depth [m]')
        ax3.set_ylabel('porosity []')
        plt.savefig('ec_poro_depth.png')
        plt.close()

        return 0

    def fill_nonuniform_grid(self, e_cells):
        print "Filling nonuniform grid" 
        count = 0
        k = 0
        columnz = []
        columnk = []
        columnphi = []
        check_col = self.nz
        check_plane = self.nx*self.ny
        for j in range(self.ny-1,-1,-1):
            for i in range(0,self.nx):
                for k in range(0,self.nz):
                    ind = (i + self.nx *j) + check_plane*k
                    if i == 32 and j == 77:
                        print k, e_cells[ind].getZPermeability()
                    if e_cells[ind].getXPermeability() > 1:
                        columnz.append(e_cells[ind].getTopZ())
                        columnk.append(e_cells[ind].getXPermeability())
                        columnphi.append(e_cells[ind].getPorosity())
                # spits out the averages after the column index is filled.     
                kmean, kvar = stats(columnk)
                if self.homogeneous == True:
                    k_write = self.permval
                else:
                    k_write = kmean
                if self.poroval == False:
                    phimean , phivar = stats(columnphi)
                else:
                    phimean = self.poroval

                top_b = -columnz[0]
                bottom_b = -columnz[-1]

                x = e_cells[ i + self.nx * j].getCenterX()
                y = e_cells[ i + self.nx * j].getCenterY()

                # get correct dx and dy
                if i != (self.nx-1):
                    x_1 = e_cells[i+1 + self.nx*j].getCenterX()
                    x_0 = e_cells[ i + self.nx * j].getCenterX()
                    dx = (x_1 - x_0)
                else:
                    x_1 = e_cells[i + self.nx*j].getCenterX()
                    x_0 = e_cells[ i-1 + self.nx * j].getCenterX()
                    dx = (x_1 - x_0)
                if j != 0:
                    y_1 = e_cells[i + self.nx *(j-1)].getCenterY()
                    y_0 = e_cells[i + self.nx * j].getCenterY()
                    dy = y_1 - y_0
                else: 
                    y_1 = e_cells[i + self.nx *j].getCenterY()
                    y_0 = e_cells[i + self.nx *(j+1)].getCenterY()
                    dy = y_1 - y_0

                #boundary condition key [integers]
                # 1 = internal
                # 3 = constant pressure
                # 4 = no flow
                if i == (self.nx - 1):
                    east_bc = 3
                else:
                    east_bc = 1

                if j == 0:
                    north_bc = 3
                else:
                    north_bc = 1

                if i == 0:
                    west_bc = 3
                else: 
                    west_bc = 1
                    
                if j == (self.ny - 1):
                    south_bc = 3
                else:
                    south_bc = 1

                pressure = -self.gradient * bottom_b * 1000
                gc = GridCell(top_b, bottom_b, x, y, dx, dy, phimean, k_write,\
                        west_bc, east_bc, south_bc, north_bc, pressure)
                self.grid_cells.append(gc)

                # increment loops
                count += 1
                columnz = []
                columnk = []
                columnphi = []
                k = 0
        return 0

    def write_layer(self):
        print "writing layer " + self.layer_name
        f = open("".join([self.layer_name,'.txt']),"w")
        g = open("thickness.txt","w")
        # layer type
        f.write(''.join([str(self.l_type) + '\n']))
        # layer id
        f.write(''.join([str(self.l_id) + '\n']))
        # fluid parameters
        f.write(''.join([str(self.l_co2_rho),'\n']))
        f.write(''.join([str(self.l_bri_rho),'\n']))
        f.write(''.join([str(self.l_co2_mu),'\n']))
        f.write(''.join([str(self.l_bri_mu),'\n']))
        f.write(''.join([str(self.sc_res),'\n']))
        f.write(''.join([str(self.sb_res),'\n']))
        f.write(''.join([str(self.c_co2),'\n']))
        f.write(''.join([str(self.c_bri),'\n']))
        f.write(''.join([str(self.c_roc),'\n']))
        if self.cap_rp_id == 0:
            f.write(''.join([str(self.cap_rp_id),'\n']))
        elif self.cap_rp_id == 1:
            lamb = 3.
            p_entry = 3000000.
            f.write(''.join([str(self.cap_rp_id), ', ',\
                    str(lamb), ', ', \
                    str(p_entry), '\n']))
            self.plot_cap_rp_bc(lamb, p_entry)
        # number of cells
        f.write(''.join([str(self.nx * self.ny), '\n']))
        for cel in self.grid_cells:
            g.write(''.join([str((cel.top_b - cel.bottom_b)),', ']))
            cel.write_cell(f)
        f.close()
        return 0
    def bc_cap(self, pentry, lamb):
        return pcap
    def plot_cap_rp_bc(self, lamb, p_entry):
        sb = np.linspace(self.sb_res,1.)
        pc = np.zeros(len(sb))
        krb = np.zeros(len(sb))
        krc = np.zeros(len(sb))
        for i in range(len(sb)):
            seff = (sb[i] - self.sb_res) / (1 - self.sb_res)
            pc[i] = p_entry * pow(seff, -1/lamb)
            krb[i] = pow(seff, (2 + 3 * lamb)/lamb)
            krc[i] = (1 - seff)**2 * (1 - pow(seff, (2 + lamb) / lamb))
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.plot(sb, pc)
        ax1.set_xlabel('sb []')
        ax1.set_ylabel('pcap [Pa]')
        plt.savefig('pcap.png')
        plt.clf()
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(sb, krb, label = 'krb')
        ax2.plot(sb, krc, label = 'krc')
        ax2.legend(loc=1)
        ax2.set_xlabel('sb')
        plt.savefig('relperm.png')
        return 0

class GridCell(object):
    def __init__(self, top_b, bottom_b, x, y, dx, dy, phi, k,\
            west_bc, east_bc, south_bc, north_bc, pressure): 
        self.top_b = top_b
        self.bottom_b = bottom_b
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy
        self.phi = phi
        self.k = k
        self.west_bc = west_bc
        self.east_bc = east_bc
        self.south_bc = south_bc
        self.north_bc = north_bc
        self.pressure = pressure

    def write_cell(self, f):
        f.write(''.join([str(self.x),', ']))
        f.write(''.join([str(self.y),', ']))
        f.write(''.join([str(self.dx),', ']))
        f.write(''.join([str(self.dy),', ']))
        # porosity
        f.write(''.join(['%.4f' % self.phi ,', ']))
        # NOTE: Enters permeability for all dimensions since the formation is isotropic
        # in the planar directions and the z permeability is not used in VESA.
        xperm = self.k
        yperm = self.k 
        zperm = self.k
        f.write(''.join(['%.0f' % xperm, ', ' ,'%.0f' % yperm, ', ',\
                '%.0f' % zperm, ', ']))
        f.write(''.join(['%.4f'  % self.bottom_b, ', ']))
        f.write(''.join(['%.4f' % self.top_b, ', ']))
        # initial CO2 saturation
        f.write('0.0, ')
        # past CO2 saturation
        f.write('0.0, ')
        #initial pressure at bottom of formation [Pa]
        f.write(''.join(['%.0f' % self.pressure, ', ']))
        f.write(''.join([str(self.east_bc), ', ']))
        f.write(''.join([str(self.north_bc), ', ']))
        f.write(''.join([str(self.west_bc), ', ']))
        f.write(''.join([str(self.south_bc), ', ']))
        f.write('\n')
        return f

def stats(data):
    sum_s = 0.0
    for value in data:
        sum_s += value
    mean = sum_s/len(data)
    sum_s = 0.0
    for value in data:
        sum_s += (value - mean)**2
    variance = sum_s/(len(data)-1)
    return(mean,variance)
