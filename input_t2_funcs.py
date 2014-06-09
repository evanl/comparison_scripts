#Author - Evan Leister
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import string
def write_separator(f, keyword):
    """ writes for the following keywords
    START ROCKS
    MULTI PARAM
    SOLVR GENER
    TIMES FOFT 
    GOFT  COFT 
    ELEME INCON
    ENDCY CONNE
    MESHMAKER
    SELEC MOP
    ENDFI
    one space is required after four letter keywords
    none required after MOP
    """
    if keyword == 'ENDFI':
        f.write(''.join([keyword,'----']))
        f = writeDashes(f)
    elif keyword == 'MOP':
        f.write('----*----1 MOP: 123456789*123456789*1234')
        f.write(' ---*----5----*----6----*----7----*----8')
        f.write('\r\n')
    elif keyword == 'SELEC':
        f.write('SELEC....2....3....4....5....6....7....8')
        f.write('....9...10...11...12...13...14...15...16')
        f.write('\r\n')
    elif keyword =='MESHMAKER':
        f.write('MESHMAKER')
        writeDashes(f)
        f.write('\r\n')
    elif (  
            keyword == 'START' or 
            keyword == 'ROCKS' or 
            keyword == 'MULTI' or 
            keyword == 'PARAM' or 
            keyword == 'SOLVR' or 
            keyword == 'GENER' or 
            keyword == 'TIMES' or 
            keyword == 'FOFT ' or 
            keyword == 'GOFT ' or 
            keyword == 'COFT ' or
            keyword == 'ELEME' or 
            keyword == 'CONNE' or 
            keyword == 'INCON' or  
            keyword == 'ENDCY'    ):
        f.write(''.join([keyword,'----']))
        f = writeDashes(f)
        f.write('\r\n')
    else :
        print "keyword = " + keyword
        raise NameError('Invalid Keyword')
    return f

def writeDashes(f):
    """ f must be an open file in w mode: writes this: 
    1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
    """
    f.write('1----*----2----*----3----*----4----*----5----*----6----*----7----*----8')
    return f

def write_rocks(f, name, density, porosity, xperm, yperm, zperm, \
    cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
    cap = 'vanGenuchten', rel_perm = 'vanGenuchten',\
    end = False):
    f.write(name + '    2')
    d = '{: 10.0f}'.format(density)
    p = '{: 10.2f}'.format(porosity)
    xp = '{: 10.2e}'.format(xperm)
    yp = '{: 10.2e}'.format(yperm)
    zp = '{: 10.2e}'.format(zperm)
    l = '{: 10.2f}'.format(thermk)
    sph = '{: 10.2f}'.format(specheat)
    f.write(d + p + xp + yp + zp + l + sph + '\r\n')

    # compressibility = 0
    f.write('   0.0e-10\r\n')

    # Capillary pressure
    if cap == 'vanGenuchten':
        cp_type = 7
    elif cap == 'none':
        cp_type = 8
    else:
        cp_type = 1

    # relative permeability
    if rel_perm == 'vanGenuchten':
        rp_type = 7
    else:
        rp_type = 1

    # write either form
    rp_str = 4 * ' ' + str(rp_type) + 5 * ' '
    s = ""
    for el in rp_vals:
        s += format_float_rocks(el)
    f.write(rp_str + s + '\r\n')
    cp_str = 4 * ' ' + str(cp_type) + 5 * ' '
    sc = ""
    for el in cp_vals:
        sc += format_float_rocks(el)
    f.write(cp_str + sc + '\r\n')

    # break the rocks line
    if end == True:
        f.write('\r\n')
    return f

def plot_relperm_cap(rp_vals, cp_vals, fmt = 'png',\
        rp = 'linear', cp = 'linear'):

    print "PLOTTING RELATIVE PERMEABILITY AND CAPILLARY PRESSURE CURVES"
    nvals = 100
    if rp == 'linear':
        sbres = rp_vals[0]
    else:
        sbres = rp_vals[1]
    sat = np.linspace(sbres, 1., nvals)
    pcap = np.zeros(nvals)
    krl = np.zeros(nvals)
    krg = np.zeros(nvals)

    if cp == 'linear':
        for i in range(len(sat)):
            pcap[i] = cap_linear(sat[i], cp_vals)
    elif cp == 'none':
        for i in range(len(sat)):
            pcap[i] = 0.
    else:
        for i in range(len(sat)):
            pcap[i] = -cap_vangenuchten(sat[i], cp_vals)
    if rp == 'linear':
        for i in range(len(sat)):
            krl[i], krg[i] = rel_perms_linear(sat[i], rp_vals)
    else:
        for i in range(len(sat)):
            krl[i], krg[i] = rel_perms_vangenuchten(sat[i], rp_vals)

    font = { 'size'   : 14}
    matplotlib.rc('font', **font)
    f = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    #f.suptitle('Capillary Pressure Curve')
    ax = f.add_subplot(111)
    ax.set_xlabel('Sw []')
    ax.set_ylabel('Capillary Pressure [Pa]')
    p = plt.plot(sat, pcap)
    f.savefig('capillary_pressure' + '.' + fmt)
    plt.clf()
    plt.close()

    g = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    #g.suptitle('Relative Permeability Curves')
    ax = g.add_subplot(111)
    ax.set_xlabel('Sw []')
    ax.set_ylabel('Relative Permeability []')
    p = plt.plot(sat, krl, label = 'liquid')
    p = plt.plot(sat, krg, label = 'gas')
    plt.legend()
    g.savefig('rel_perm_curves' + '.' + fmt)
    plt.clf()
    plt.close()

    return 0

def cap_linear(s, cp_vals):
    if cp_vals[2] < cp_vals[1]:
        print 'LINEAR CAPILLARY PRESSURE BOUND DONT MAKE NO SENSE'
        return 1
    if s <= cp_vals[1]:
        pcap = -cp_vals[0]
    elif s >= cp_vals[2]:
        pcap = 0
    else: 
        pcap = -cp_vals[0] *(cp_vals[2] - s) / (cp_vals[2] - cp_vals[1])
    return pcap

def cap_vangenuchten(sl, cp_vals):
    lamb = cp_vals[0]
    slr = cp_vals[1]
    p_0 = 1. / cp_vals[2]
    pmax = cp_vals[3]
    sls = cp_vals[4]

    ss = (sl - slr) / (sls - slr)

    pcap = -p_0 * pow(pow(ss,-1./lamb) - 1, 1. - lamb)
    if pcap < -pmax:
        pcap = -pmax
    elif pcap > 0.:
        pcap = 0.

    return pcap

def rel_perms_linear(sl, rp_vals):
    sg = 1. - sl
    l_slope = 1 / (rp_vals[2] - rp_vals[0])
    g_slope = 1 / (rp_vals[3] - rp_vals[1])

    if sl <= rp_vals[0]:
        krl = 0.
    elif sl >= rp_vals[2]:
        krl = 1.
    else:
        krl = (sl - rp_vals[0]) * l_slope

    if sg <= rp_vals[1]:
        krg = 0.
    elif sg >= rp_vals[3]:
        krg = 1.
    else:
        krg = (sg - rp_vals[1]) * g_slope

    return krl, krg

def rel_perms_vangenuchten(sl, rp_vals):
    lamb = rp_vals[0]
    s_lr = rp_vals[1]
    s_ls = rp_vals[2]
    s_gr = rp_vals[3]

    ss = (sl - s_lr) / (s_ls - s_lr)
    sh = (sl - s_lr) / ( 1. - s_lr - s_gr)

    if sl < s_ls:
        krl = np.sqrt(ss) * pow(1 - pow(1 - pow(ss, 1./lamb), lamb), 2.)
    else:
        krl = 1.

    if s_gr == 0.:
        krg = 1. - krl
    elif s_gr > 0.:
        krg = pow(1. - sh, 2.) * (1 - pow(sh, 2.))
    else: 
        print "s_gr < 0 in rel_perms_vangenuchten"
        return 1

    return krl, krg

def write_multi(f, mode = 'isothermal'):
    if mode == 'isothermal':
        write_separator(f, 'MULTI')
        f.write('    3    3    3    6\r\n')
    else:
        write_separator(f,'MULTI')
        f.write('    3    4    3    6\r\n')
    return f

def write_selec(f):
    write_separator(f, 'SELEC')
    f.write('    1     ' \
        + 10 * ' ' + 10 * ' ' + 10 * ' ' + \
          '         0' + '    0    0' + \
                '    0    0' + '    0    0\r\n')
    f.write('        .8        .8\r\n')
    return f

def write_start(f):
    write_separator(f,'START')
    write_separator(f,'MOP')

    return f

def write_param(f, init_check = False, pres = 110.5e5, salt = 3.2e-2, co2 = 0.0, temp = 37., \
    tmax = 63.1152e6, tolexp = -7):
    write_separator(f, 'PARAM')

    # line 1
    maxtstep = '9500'
    
    f.write('    ' + maxtstep + 2 * ' '   + '  99991000')
    if init_check == True:
        f.write(' 09000000 ')
    else :
        f.write(' 00000000 ')
    f.write(' 4    3   \r\n')

    # line 2
    f.write('          ')
    
    b = '{: 10.5e}'.format(tmax)
    l = list(b)
    # l[2], l[1] = l[1], l[2]
    del l[-2]
    del l[-2]
    c = ''.join(l)
    f.write(c)
    f.write('     0.8e6')
    #f.write('       -1.')
    #                         gravity
    f.write(20 * ' ' + '      9.81\r\n')
    #f.write('      1.e2\r\n')
    tolstr = '{:03d}'.format(tolexp)
    f.write('    1.E' + tolstr + '    1.E-01\r\n')

    # default initial conditions
    s1 = '{: 10.3e}'.format(pres)
    s2 = '{: 10.3e}'.format(salt)
    s3 = '{: 10.3e}'.format(co2 )
    s4 = '{: 10.1f}'.format(temp)
    f.write(10 * ' ' + s1)
    f.write(10 * ' ' + s2)
    f.write(10 * ' ' + s3)
    f.write(10 * ' ' + s4)
    f.write('\r\n')
    return f

def write_solvr(f, linsolve_int, preprocess_int, tolexp = -7):
    write_separator(f,'SOLVR')
    if tolexp < -12 or tolexp > -6:
        print "tolexp must be between -12 and -6"
        return 1
    tolstr = '{:03d}'.format(tolexp)
    if linsolve_int < 2 or linsolve_int > 6:
        print "pick a correct linear solver"
        return 1
    f.write(str(linsolve_int) + \
            '  Z'+ str(preprocess_int) + \
            '   O0    1.0e-1   1.0e' + tolstr +'\r\n')
    return f

def write_gener(f, eleme, phase = 'brine', mass_rate = .1585, column_inj = False,\
        kg_inflow = [], times = []):
    """
            So far only writes one constant rate. Need to fix later
    """
    write_separator(f,'GENER')
    if phase =='co2':
        p = 'COM3 '
    else:
        p = 'COM1 '

    if kg_inflow == [] :
        if column_inj == False:
            mr = '{: 10.4f}'.format(mass_rate)
            f.write(eleme + 'inj 1'+ 19 * ' ' + '1' + 5 * ' ' + p  + mr  + '\r\n')
        else:
            print mass_rate / len(column_inj)
            mr = '{: 10.4f}'.format(mass_rate / len(column_inj))
            for i in range(len(column_inj)):
                ir = '{:2d}'.format(i+1)
                f.write(column_inj[i] + 'inj' + ir + 19 * ' '  +\
                        '1' + 5 * ' ' + p  + mr  + '\r\n')
    else: 
        if len(kg_inflow) != len(times):
            print "mass rate and time arrays are not the same length"
            return 1

        # store the formatted strings in lists
        massprint = []
        timeprint = []
        for i in range(len(kg_inflow)):
            b = '{: 15.3e}'.format(times[i] * 24. * 3600.)
            l = list(b)
            del l[-2]
            if l[-2] == '+':
                l[-2] = '0'
            c = ''.join(l)
            timeprint.append(' ' + c)
            
            # quick and dirty formatting of mass inflow. 
            # only works if the mass flow rate is between -1 and +1. 
            b = '{: 15.3e}'.format(kg_inflow[i])
            l = list(b)
            del l[-2]
            if l[-2] == '+':
                l[-2] = '0'
            c = ''.join(l)
            massprint.append(' ' + c)

        # LTAB parameter
        ltab = '{:5d}'.format(len(kg_inflow))

        f.write(eleme + 'inj 1'+ 10 * ' '  + 5 * ' ' \
                + ltab + 5 * ' ' + p   + '\r\n')

        # write F1: Generation times
        for i in range(1, len(timeprint)+1):
            f.write(timeprint[i-1])
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')

        #write F2: Generation rates
        for i in range(1, len(massprint)+1):
            f.write(massprint[i-1])
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')

        # write F3: specific enthalpy
        # specific enthalpy of c02 chosen to be 50.
        # shouldn't matter since it's isothermal.
        specific_enthalpy = '{: 16.1e}'.format(50.)
        se = list(specific_enthalpy)
        del se[-3]
        specific_enthalpy = ''.join(se)
        for i in range(1, len(timeprint)+1):
            f.write(specific_enthalpy)
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')

    #close
    f.write('\r\n')
    return f

def write_times(f, timelist):
    numtimes = len(timelist)

    write_separator(f,'TIMES')
    f.write('   '+ str(numtimes)+'\r\n')

    count = 0
    for i in range(numtimes):
        seconds = timelist[i] * 24. * 3600.
        b = '{: 10.4e}'.format(seconds)
        l = list(b)
        del l[-2]
        del l[-2]
        c = ''.join(l)
        f.write(' ' +  c)
        count +=1
        if count == 8: 
            f.write('\r\n')
            count = 0

    f.write('\r\n')

    return f

def write_foft(f):
    write_separator(f, 'FOFT ')
    # place element ID's here to get some output
    f.write('\r\n')
    return f

def write_coft(f, sleipner):
    write_separator(f, 'COFT ')
    if sleipner == True:
        #f.write('bH732cH732\r\n')
        dummy =1
    else:
        #f.write('aB212bB212\r\n')
        f.write('xB212yB212\r\n')
    f.write('\r\n')
    return f

def write_goft(f):
    write_separator(f, 'GOFT ')
    f.write('\r\n')
    return f

def write_meshmaker(f , rect = True, flat = True, nx = 65, ny = 119, nz = 1):
    write_separator(f, 'MESHMAKER')
    f.write('XYZ\r\n')
    # degrees 
    f.write('        0.\r\n')

    f.write('NX     '+ str(nx)+'\r\n')
    count = 0
    for i in range(nx):
        f.write('       50.')
        count += 1
        if count % 8 == 0:
            f.write('\r\n')
    f.write('\r\n')

    f.write('NY     '+ str(ny) + '\r\n')
    count = 0
    for i in range(ny):
        f.write('       50.')
        count += 1
        if count % 8 == 0:
            f.write('\r\n')
    f.write('\r\n')
    if flat == True:
        f.write('NZ       1      28.0\r\n') 
    else:
        f.write('NZ     '+ str(nz) + '\r\n')
        count = 0
        for i in range(nz):
            f.write('       50.')
            count +=1
            if count % 8 == 0:
                f.write('\r\n')

    f.write('\r\n')
    f.write('\r\n') 
    return f

def format_float_mesh(val):
    """ takes in a double, returns a string to be entered into the MESH file
            ( output s   , val )
            ('-.5000E-00', -0.5)
            ('0.1000E+01', 1.0)
            ('0.3500E+03', 350.0)
            ('-.7500E-03', -0.00075)
            ('0.4100E-01', 0.041)
    """
    a = '{: .3E}'.format(val)
    l = list(a)
    l[1], l[2] = l[2], l[1]
    if l[-3] == '+':
        l[-1] = str( int(l[-1]) + 1 )
    elif l[-3] == '-':
        l[-1] = str( int(l[-1]) - 1 )
    else : 
        print "What you input was not a float or somthin"
        print "writingFunctions.format_float_mesh(val)"
        print "val = " + str(val)
    if l[0] != '-':
        l[0] = '0'
    s = ''.join(l)
    return s 

def format_float_rocks(val):
    return '{: .3e}'.format(val)
    

def format_float_incon(val):
    """ formats to fit in the INCON, includes blank chars before number
    only works for positive floating point values (I think)
    """
    a = '{:> 20.12E}'.format(val)
    l = list(a)
    l[2], l[3] = l[3], l[2]
    if l[-3] == '+':
        l[-1] = str( int(l[-1]) + 1 )
    elif l[-3] == '-':
        l[-1] = str( int(l[-1]) - 1 )
    else : 
        print "What you input was not a float or somthin"
        print "writingFunctions.format_float_mesh(val)"
        print "val = " + str(val)
    s = ''.join(l)
    return s

class T2InputGrid(object):

    def __init__(self, nx, ny, nz):
        self.elements = []
        self.boundary = []
        self.el_array = []
        self.x = {}
        self.y = {}
        self.z = {}
        self.z_top = {}
        self.z_bot = {}
        self.corners = {}
        self.mat = {}
        self.vol = {}
        self.area = {}
        self.pres = {}
        self.na_cl = {}
        self.x_co2 = {}
        self.temp = {}
        self.nx = nx
        self.ny = ny
        self.nz = nz

    class Corner(object):
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

        def get_x(self):
            return self.x

        def get_y(self):
            return self.y

        def get_z(self):
            return self.z


    def get_x_char(self, i):
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

    def get_y_char(self, j):
        if j > 259 or j < 0:
            print "change y indexing system"
            return 1
        else : 
            letters = string.uppercase
            c1 = letters[ j/10]
            c2 = str( j - 10 * (j/10))
            char = c1 + c2
            return char

    def get_z_char(self, k):
        if k > 51 or k < 0:
            print 'fix your z indecisis'
            return 1
        else:
            letters =  string.lowercase + string.uppercase
            char = letters[k]
            return char

    def get_element_chars(self, i, j, k):
        c1 = self.get_x_char( i)
        c2 = self.get_y_char( j)
        c3 = self.get_z_char( k)
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

    def fill_uniform_grid(self, porosity, dx, dy, dz, density = 1000.,\
            solubility_limit = 0.454104e-3, altered_cell = 'none'):
        # populates grid 
        self.nx_start = 0
        self.ny_start = 0
        self.nz_start = 0
        g = open('MESH', 'w')
        print "MESH created with:"
        for i in range(self.nx):
            temparray = []
            for j in range(self.ny):
                temparrayz = []
                for k in range(self.nz):
                    eleme = self.get_element_chars(i, j, k)
                    if i == 12:
                        if j == 12:
                            if k == 12:
                                print "INJECTION cell"
                                print eleme
                    self.elements.append(eleme)

                    xlocal = (dx / 2 + dx * i )
                    ylocal = (dy / 2 + dy * j ) 
                    zlocal = -(dz / 2 + dz * k ) - 859.5
                    if altered_cell != 'none':
                        #if eleme == 'cA2 3' or eleme == 'cA2 1'\
                                #or eleme == 'cA3 2' or eleme == 'cA1 2':
                        if eleme == altered_cell:
                            zlocal = zlocal + 5.
                    self.x[eleme] = xlocal
                    self.y[eleme] = ylocal
                    self.z[eleme] = zlocal
                    self.vol[eleme] = (dz * dy * dx)
                    self.area[eleme] = dy * dx
                    self.pres[eleme] = -zlocal * density * 10.0
                    self.na_cl[eleme] = 3.2e-2
                    self.x_co2[eleme] = solubility_limit
                    self.temp[eleme] = 32.
                    temparrayz.append(eleme)

                    #if k ==1 :
                        #self.mat[eleme] = 'sands'
                    #else: 
                        #self.mat[eleme] = 'shale'
                    self.mat[eleme] = 'sands'

                    vw = format_float_mesh(self.vol[eleme])
                    aw = format_float_mesh(self.area[eleme])
                    xw = format_float_mesh(self.x[eleme])
                    yw = format_float_mesh(self.y[eleme])
                    zw = format_float_mesh(self.z[eleme])
                temparray.append(temparrayz)
            self.el_array.append(temparray)

        return 0


    def e_cell_index(self,i,j,k):
        """ returns the index used to call the cell list"""
        nx = 65
        ny = 119
        ind = i + nx * j + nx * ny * k
        return ind 

    def fill_3d_grid(self, e_cells , temperature = 37., density = 1000.,\
            two_d = False,  gradient = 10 , solubility = 0.474e-1, \
            five_section = [], shale = True):
        """Fills the 3d grid with x, y, and z
        'gradient' specifies the hydrostatic gradient. [MPa/km]

        'solubility' specifies the initial dissolved CO2 in aqueous phase.

        'two_d', 
        if true, ensures that the grid being created is two-dimensional.

        'thickness', if specified, returns a constant thickness, also requires 
        'depth' to be specified in order to get the correct initial pressure.
        'gradient' is the pressure gradient, given in units of MPa/km

        if thickness and depth are not specified, 
        the sleipner vertical geometry is considered.
        """
        print "FILLING 3D grid based on Sleipner input data........"

        if five_section != []:
            self.nx_start = five_section[0]
            self.ny_start = five_section[1]
            self.nz_start = five_section[2]
            self.nx = self.nx_start + 5
            self.ny = self.ny_start + 5
            self.nz = self.nz_start + 5
        else:
            self.nx_start = 0
            self.ny_start = 0
            self.nz_start = 0

        if two_d == True:
            nzi = 1
        else:
            nzi = self.nz

        count = 0
        for i in range(self.nx_start, self.nx):
            temparray = []
            for j in range(self.ny-1, self.ny_start-1, -1):
                temparrayz = []
                sand_count = 0
                for k in range(self.nz_start, nzi):
                    ind = self.e_cell_index(i, j, k)
                    eleme = self.get_element_chars(i, j, k)
                    if shale == True or e_cells[ind].getXPermeability() > 1.:
                        if e_cells[ind].getXPermeability() > 1.:
                            sand_count +=1
                        count +=1
                        self.elements.append(eleme)
                        self.write_t2_cell(e_cells, i, j, k,\
                                temperature, density,\
                                two_d, gradient, solubility)
                        #el_array options
                        temparrayz.append(eleme)
                temparray.append(temparrayz)
            self.el_array.append(temparray)
        print "GRID FILLING COMPLETE"
        print str(count) + " ELEMENTS CREATED"
        print len(self.el_array), len(self.el_array[0]), len(self.el_array[0][0])
        return 0

    def write_t2_cell(self, e_cells, i, j, k,\
            temperature = 37., density = 1000.,\
            two_d = False,  gradient = 10 , solubility = 0.474e-1):
        ind = self.e_cell_index(i, j, k)
        eleme = self.get_element_chars(i, j, k)
        # sets material id
        if e_cells[ind].getXPermeability() > 1 or two_d == True:
            self.mat[eleme] = 'sands'
        else:
            self.mat[eleme] = 'shale'

        oc = e_cells[ind].getCorners()
        corners = []
        for c in oc:
            x, y = c.getXY()
            # FLIPPING ALL ZS IN THIS. 
            z = - c.getZ()
            nc = self.Corner(x, y, z)
            corners.append(nc)
        self.corners[eleme] = corners

        # getting more precise center and volume
        x = self.get_x_centroid(corners)
        y = self.get_y_centroid(corners)
        z = self.get_z_centroid(corners)
        volume = self.get_volume(x, y, z, corners)
        if two_d == True:
            pc_count = []
            for k in range(self.nz):
                pc_ind = self.e_cell_index(i,j,k)
                pc_ind = self.e_cell_index(i,j,k)
                if e_cells[pc_ind].getXPermeability() > 1.:
                    pc_ind_1 = self.e_cell_index(i,j,k+1)
                    pc_count.append(k)
                    ztop = -e_cells[pc_ind].getTopZ()
                    zbot = -e_cells[pc_ind_1].getTopZ()
                    dz =  ztop - zbot
                    if i == 32 or i == 25:
                        if j == 77:
                            print "i, j, k, dz, ", i, j, k, dz
            bot_ind = self.e_cell_index(i, j, pc_count[-1])
            top_ind = self.e_cell_index(i, j, pc_count[0])
            #volume = dx * dy * dz
            #new idea
            tc = e_cells[top_ind].getCorners()
            bc = e_cells[bot_ind].getCorners()
            corners = []
            for c_ind in range(4):
                x, y = tc[c_ind].getXY()
                z = -tc[c_ind].getZ()
                nc = self.Corner(x, y, z)
                corners.append(nc)
            for c_ind in range(4,8):
                x,y = bc[c_ind].getXY()
                z = -bc[c_ind].getZ()
                nc = self.Corner(x,y,z)
                corners.append(nc)
            self.corners[eleme] = corners
            x = self.get_x_centroid(corners)
            y = self.get_y_centroid(corners)
            z = self.get_z_centroid(corners)
            volume = self.get_volume(x, y, z, corners)
        self.x[eleme] = x
        self.y[eleme] = y
        self.z[eleme] = z
        pressure = -z * gradient * density
        self.vol[eleme] = volume
        self.area[eleme] = 0.0
        self.pres[eleme] = pressure
        self.na_cl[eleme] = 3.2e-2
        self.x_co2[eleme] = solubility
        self.temp[eleme] = temperature

    def get_x_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_x()
        return sum_c / count

    def get_y_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_y()
        return sum_c / count

    def get_z_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_z()
        return sum_c / count
  
    def get_dx(self, eleme, direc):
        """ returns the length of a grid cell in a particular direction.
        dir is either 1, 2 or 3 for x, y and z directions.
        i, j and k are the indices
        """
        if direc == 1 :
            corners = self.corners[eleme]
            dx = corners[0].get_x() - corners[1].get_x()
            return dx
        elif direc == 2 :
            corners = self.corners[eleme]
            dy = corners[0].get_y() - corners[2].get_y()
            return dy
        elif direc == 3 :
            z1 = abs(e_cells[self.e_cell_index(i,j,k)].getTopZ() - \
                    e_cells[self.e_cell_index(i,j,k)].getBottomZ())
            return z1
        else:
            raise Exception("Invalid direction, \n" + \
                    " Please specify 1, 2 or 3.\n")

    def get_volume(self, x, y, z, corners):
        """ uses the equation for volume of an orientable polyhedron
            V = 1/3 \sum_i x_i \dot n^hat_i A_i
        """
        face_map = ['west', 'south', 'east', 'north', 'bot', 'top']

        v_sum = 0.0
        for face in face_map:
            a = self.get_area(corners, face)
            centroid = self.get_face_center(x, y, z, corners, face)
            cent = np.asarray(centroid)
            vec = self.get_normal_vector(x, y, z, corners, face)
            v_sum += np.dot(cent, vec) * a

        vol = 1./3. * v_sum
        return vol

    def get_area(self, corners, face):
        """ returns the area of a cell face, east, west, etc
        """ 
        if face == 'west':
            x1 = corners[2].get_y()
            x2 = corners[0].get_y()
            y1 = corners[2].get_z()
            y2 = corners[0].get_z()
            y3 = corners[6].get_z()
            y4 = corners[4].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'south':
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_z()
            y2 = corners[3].get_z()
            y3 = corners[6].get_z()
            y4 = corners[7].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'east':
            x1 = corners[3].get_y()
            x2 = corners[1].get_y()
            y1 = corners[3].get_z()
            y2 = corners[1].get_z()
            y3 = corners[7].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'north':
            x1 = corners[0].get_x()
            x2 = corners[1].get_x()
            y1 = corners[0].get_z()
            y2 = corners[1].get_z()
            y3 = corners[4].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_y()
            y2 = corners[0].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[6].get_x()
            x2 = corners[7].get_x()
            y1 = corners[6].get_y()
            y2 = corners[4].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")

        return area

    def get_face_center(self, xc, yc, zc, corners, face):
        """ center vector location relative to polyhedron center
        """
        if face == 'west':
            nc = [corners[0], corners[2], corners[4], corners[6]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'south':
            nc = [corners[2], corners[3], corners[6], corners[7]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
            a = 2
        elif face == 'east':
            nc = [corners[3], corners[1], corners[7], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'north':
            nc = [corners[0], corners[1], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")

        vec = [xf - xc, yf - yc, zf - zc]
        return vec

    def get_normal_vector(self, x, y, z, corners, face):
        """ gets normal vector of face
        """
        if face == 'west':
            vec = [-1., 0., 0.]
        elif face == 'south':
            vec = [0., -1., 0.]
        elif face == 'east':
            vec = [1., 0., 0.]
        elif face == 'north':
            vec = [0., 1., 0.]
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [c[0]/mag, c[1]/mag, -1./mag]
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [-c[0]/mag, -c[1]/mag, 1./mag]
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")
        return vec

    def fit_plane(self, corners):
        """ takes four corner points and fits a plane least squares to them
            returns in form z = c[0] x + c[1] y + c[2]
        """
        x = []
        y = []
        z = []
        for c in corners:
            x.append(c.get_x())
            y.append(c.get_y())
            z.append(c.get_z())
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        A = np.column_stack((x, y, np.ones(x.size)))
        c, resid, rank, sigma = np.linalg.lstsq(A, z)
        return c, resid, rank, sigma

    def get_area_side(self, x1, x2, y1, y2, y3, y4):
        h = x2 - x1
        b1 = y4 - y2
        b2 = y3 - y1
        return 0.5 * h * (b1 + b2)
        
    def write_mesh(self, e_cel, two_d = False, uniform = False,\
            boundary_type = 1, shale = True,\
            type1_source_cell = 'none'):
        """ populates grid and writes ELEME block of MESH
        """
        g = open('MESH', 'w')
        print "MESH created with:"
        print "writing ELEMENT block of input data"
        g.write('ELEME\r\n')
        if two_d == True:
            nzi = 1
        else:
            if shale == False:
                nzi = 34
            else:
                nzi = self.nz
        count = 0 
        for i in range(0, self.nx - self.nx_start ):
            for j in range(0, self.ny - self.ny_start ):
                for k in range(0, nzi - self.nz_start):
                    eleme = self.el_array[i][j][k]
                    # if the cell is interior, writes it.
                    # interior check is below
                    # if not, adds it to the boundary group to be 
                    # written at the end
                    if boundary_type == 1 and \
                            (i == 0 or i == (self.nx - self.nx_start - 1 ) or \
                            j == 0 or j == (self.ny - self.ny_start - 1 )):
                        self.boundary.append(eleme)
                    elif boundary_type == 1 and type1_source_cell == eleme:
                        print "YAY"
                        print "element ", eleme, " is on the boundary!~"
                        print "YAY"
                        self.boundary.append(eleme)
                    else: 
                        vw = format_float_mesh(self.vol[eleme])
                        aw = format_float_mesh(self.area[eleme])
                        xw = format_float_mesh(self.x[eleme])
                        yw = format_float_mesh(self.y[eleme])
                        zw = format_float_mesh(self.z[eleme])
                        mat = self.mat[eleme]
                        g.write(eleme + 5 * ' ' + 5 * ' ' + mat + vw + aw + \
                            10 * ' ' + xw + yw + zw + '\r\n')
                        count +=1

        if boundary_type == 1:
            g.write('ina\r\n')

        for el in self.boundary:
            vw = format_float_mesh(self.vol[el])
            aw = format_float_mesh(self.area[el])
            xw = format_float_mesh(self.x[el])
            yw = format_float_mesh(self.y[el])
            zw = format_float_mesh(self.z[el])
            mat = self.mat[el]
            g.write(el + 5 * ' ' + 5 * ' ' + mat + vw + aw + \
                10 * ' ' + xw + yw + zw + '\r\n')
            count += 1 
        g.write('\r\n')

        print "ELEME: " + str(count) + " elements"
        print "ELEMENTS COMPLETE --------------- "
        print "Writing Connections.........."
        if two_d == True:
            nzi = 1
        else:
            if shale == False:
                nzi = 34
            else:
                nzi = self.nz
        g.write('CONNE\r\n')
        count = 0 
        for i in range(0, self.nx - self.nx_start ):
            for j in range(0, self.ny - self.ny_start ):
                for k in range(0, nzi - self.nz_start):
                    # z connection
                    if k != nzi - self.nz_start -1:
                        g, count = self.write_z_connection(g, count, \
                                i, j, k, uniform)
                    # y connection
                    if j != self.ny - self.ny_start -1:
                        g, count = self.write_y_connection(g, count, \
                                i, j, k, two_d, uniform, e_cel )
                    # x connection
                    if i != self.nx - self.nx_start -1: 
                        g, count = self.write_x_connection(g, count, \
                                i, j, k, two_d, uniform, e_cel )

        g.write('\r\n')
        g.write('\r\n')
        
        g.close()
        print "CONNE: " + str(count) +  " connections"
        print "CONNECTIONS COMPLETE -----------------"
        return 0

    def write_z_connection(self, g, count, i, j, k, uniform = False):
        """ should be legit now. 
        """
        if uniform == False:
            direc = 3  
            el1 = self.el_array[i][j][k]
            el2 = self.el_array[i][j][k+1]
            z1 = self.z[el1]
            z2 = self.z[el2]
            area = self.get_area(self.corners[el1], 'bot')
            mid = self.get_face_center(self.x[el1], self.y[el1], self.z[el1],\
                    self.corners[el1], 'bot')
            zmid = z1 + mid[2]
            dz1 = z1 - zmid
            dz2 = zmid - z2
            dz1w = format_float_mesh(dz1)
            dz2w = format_float_mesh(dz2)
            aw = format_float_mesh(area)
            beta = 1.0
            betaw = format_float_mesh(beta) 
            z_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dz1w + dz2w + aw + betaw + '\r\n'
        else:
            direc = 3  
            el1 =  self.el_array[i][j][k]
            el2 =  self.el_array[i][j][k+1]
            z1 = self.z[el1]
            z2 = self.z[el2]
            dzn = (z1 - z2)/2.
            dz1 = format_float_mesh((z1 - 0.5 * ( z1 + z2)))
            dz2 = format_float_mesh(-(z2 - 0.5 * (z1 + z2)))
            a1 = self.vol[el1] / dzn
            a2 = self.vol[el2] / dzn
            a = min(a1, a2)
            dzw = format_float_mesh(dzn) 
            aw = format_float_mesh(a)
            beta = 1.0
            betaw = format_float_mesh(beta) 
            z_string = el1 + el2 + 19 * ' ' + str(direc) + \
                            dz1 + dz2 + aw + betaw + '\r\n'
        g.write(z_string)
        count +=1
        return g, count

    def write_y_connection(self, g, count, i, j, k, \
            two_d = False, uniform = False, e_cells = 'uniform' ):
        """ need to incorporate z, inside this loop
        """ 
        if two_d == True:
            nzi = 1
        else:
            nzi = self.nz
        direc = 2  
        el1 =  self.el_array[i][j][k]
        el2 = self.el_array[i][j+1][k]
        if uniform == True or two_d == True:
            y1 = self.y[el1]
            y2 = self.y[el2]
            dyn = (y2 - y1)/2.
            a1 = self.vol[el1] / dyn
            a2 = self.vol[el2] / dyn
            a = min(a1, a2)
            dyw = format_float_mesh(dyn) 
            aw = format_float_mesh(a)
            z1 = self.z[el1]
            z2 = self.z[el2]
            dz = z1 - z2
            beta = dz / np.sqrt(pow(dz,2) + pow(dyn,2))
            betaw = format_float_mesh(beta) 
            y_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dyw + dyw + aw + betaw + '\r\n'
            g.write(y_string)
            count +=1
        else:
            area = self.get_area(self.corners[el1], 'south')
            dy1 = self.y[el1] - self.corners[el1][3].get_y()
            dy2 = self.corners[el2][1].get_y() - self.y[el2] 
            dy1w = format_float_mesh(dy1)
            dy2w = format_float_mesh(dy2)
            aw = format_float_mesh(area)
            dz = self.z[el1] - self.z[el2]
            beta = dz / np.sqrt(pow(dz,2) + pow(dy1 + dy2,2))
            betaw = format_float_mesh(beta)
            y_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dy1w + dy2w + aw + betaw + '\r\n'
            g.write(y_string)
            count +=1
            
        return g, count

    def write_x_connection(self, g, count, i, j, k, \
            two_d = True, uniform = False, e_cells = 'uniform' ):
        """ need to incorporate z, note, inside this loop
        """
        if two_d == True:
            nzi = 1
        else:
            nzi = self.nz
        direc = 1
        el1 =  self.el_array[i][j][k]
        el2 =  self.el_array[i+1][j][k]
        if uniform == True or two_d == True:
            x1 = self.x[el1]
            z1 = self.z[el1]
            x2 = self.x[el2]
            dxn = (x2 - x1)/2.
            a1 = self.vol[el1] / dxn
            a2 = self.vol[el2] / dxn
            a = min(a1, a2)
            dxw = format_float_mesh(dxn) 
            aw = format_float_mesh(a)
            z2 = self.z[el2]
            dz = z1 - z2
            beta = dz / np.sqrt(pow(dz,2) + pow(dxn,2))
            betaw = format_float_mesh(beta) 
            x_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dxw + dxw + aw + betaw + '\r\n'
            g.write(x_string)
            count +=1
        else:
            area = self.get_area(self.corners[el1], 'east')
            dx1 = self.corners[el1][3].get_x() - self.x[el1]
            dx2 = self.x[el2] - self.corners[el2][2].get_x()
            dx1w = format_float_mesh(dx1)
            dx2w = format_float_mesh(dx2)
            aw = format_float_mesh(area)
            dz = self.z[el1] - self.z[el2]
            beta = dz / np.sqrt(pow(dz,2) + pow(dx1 + dx2,2))
            betaw = format_float_mesh(beta)
            x_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dx1w + dx2w + aw + betaw + '\r\n'
            g.write(x_string)
            count +=1
             
        return g, count

    def write_incon(self, porosity):
        print "Writing INCON FILE"
        h = open('INCON','w')
        h.write('INCON -- INITIAL CONDITIONS FOR ' + str(len(self.elements)) + \
            ' ELEMENTS AT TIME  .000000E+00\r\n')
        for el in self.elements:
            poro = '{: 15.8E}'.format(porosity)
            h.write(el + 10 * ' '+ poro + '\r\n')
            x1 =  format_float_incon(self.pres[el])
            x2 =  format_float_incon(self.na_cl[el])
            x3 =  format_float_incon(self.x_co2[el])
            x4 =  format_float_incon(self.temp[el])
            h.write( x1 + x2 + x3 + x4 + '\r\n')
        h.write('\r\n')
        h.write('\r\n')
        h.close()
        print "INCON COMPLETE"
        return 0 

    def use_old_incon(self, hydro_directory, type1_source_cell = 'none',\
            saturation_fraction= 0.8):
        """ note that the hydro_directory does not require the '_dir/' suffix
        """
        print "Writing INCON with SAVE file from: " + hydro_directory
        with open(hydro_directory + '_dir/' + 'SAVE','r') as f:
            content = f.readlines()
        content[-2] = '\n'
        content[-1] = '\n'
        del content[0]
        f = open('INCON','w')
        f.write('INCON -- INITIAL CONDITIONS FOR ' + str(len(self.elements)) + \
            ' ELEMENTS AT TIME  .000000E+00\n')
        saturate = False
        for line in content:
            s = line.split()
            if saturate == True:
                print "saturating"
                print line
                gas_sat = 10 +  saturation_fraction
                brine_sat = 0.1 *(1. - saturation_fraction)
                brine_res = format_float_incon(brine_sat)
                gas_sat = format_float_incon(gas_sat)
                pres = format_float_incon(float(s[0]))
                temp = format_float_incon(float(s[3]))
                line = pres + brine_res + gas_sat + temp + '\r\n'
                print "changed line"
                print line
                saturate = False
            if s != [] and len(s[0]) == 3:
                s[0] = s[0] + '_' + s[1]
            if s != [] and type1_source_cell == s[0]:
                print "found one!"
                print line
                print type1_source_cell, s[0]
                saturate = True
            f.write(line)
        f.close()
        print "INCON from SAVE complete"
        return 0
    
    def plot_cells(self, show = True, two_d = False):
        print " Creating scatter plot of all cells."
        x = []
        y = []
        z = []
        pres = []
        for el in self.elements:
            x.append(self.x[el])
            y.append(self.y[el])
            z.append(self.z[el])
            pres.append(self.pres[el])

        xs = np.asarray(x)
        ys = np.asarray(y)
        zs = np.asarray(z)
        ps = np.asarray(pres)
        if two_d == True:
            xlist = []
            ylist = []
            zlist = []
            tempzlist = []
            plist = []
            tempplist = []
            count = 0
            for n in range(len(z)):

                tempzlist.append(z[n])
                tempplist.append(pres[n])
                if n < 119:
                    ylist.append(y[n])
                if n%119 == 0:
                    xlist.append(x[n])
                count +=1
                if count%119 == 0:
                    zlist.append(tempzlist)
                    plist.append(tempplist)
                    tempzlist = []
                    tempplist = []

            xg, yg = np.meshgrid(xlist, ylist)
            zg = np.asarray(zlist)
            zg = np.transpose(zg)
            pg = np.asarray(plist)
            pg = np.transpose(pg)


        if two_d == True:
            fig = plt.figure(figsize=(7.5,10), dpi = 480, facecolor = 'w',\
                edgecolor ='k')
            ax = fig.add_subplot(111)
            CS = ax.contourf(xg, yg, pg)
            CB = plt.colorbar(CS, shrink = 0.8, extend = 'both')
            title = 'contour_cells'
        else:
            fig = plt.figure(figsize=(11.,8.5), dpi = 960, facecolor = 'w',\
                edgecolor ='k')
            ax = fig.add_subplot(111, projection = '3d')
            p1 = ax.scatter(xs,ys,zs, s=10, facecolor = (.7,.7,.7), c=zs)
            plt.gray()
            ax.set_zlabel('Z')
            title = 'scatter_cells'

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        if show == True:
            plt.show()
        plt.savefig(title + '.pdf', format = 'pdf')
        plt.clf()
        plt.close()
        print "scatter plot complete"
        return 0

    def plot_cells_slice(self, direc = 1, ind = 0, show = True):
        """direc is the direction, ind is the index of the other direction
        for example, if direc == 1 and ind == 4,
        this routine will plot the cells in the x-z plane with y index of 4
        """
        print "plotting slice"
        x_mod = []
        z_mod = []
        for el in self.elements:
            if direc == 1:
                chars = self.get_y_char(ind)
                if el[1:3] == chars:
                    x_mod.append(self.x[el])
                    z_mod.append(self.z[el])
            elif direc == 2:
                chars = self.get_x_char(ind)
                if el[3:5] == chars:
                    x_mod.append(self.y[el])
                    z_mod.append(self.z[el])
            else:
                print "get the direction right, either 1 or two"

        xs = np.asarray(x_mod)
        zs = np.asarray(z_mod)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_aspect('equal')
        p1 = ax.scatter(xs,zs, s=20)
        if direc == 1:
            ax.set_xlabel('x')
        else:
            ax.set_xlabel('y')

        ax.set_ylabel('z')
        if show == True:
            plt.show()
        plt.savefig( str(direc) + '_' + str(ind) + 'scatter_cells_slice.png')
        plt.clf()
        plt.close()
        print "Scatter cells slice complete"
        return 0

if __name__ == '__main__':
    nx = 65
    ny = 119
    nz = 43
    grid = T2InputGrid(nx, ny, nz)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                char = grid.get_element_chars(i, j, k)
                if char == 'JH732':
                    print char, 'JH732'
                    print i, j, k
