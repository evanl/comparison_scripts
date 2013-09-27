import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import string
import output_t2_funcs as ot2

def read_error_elements(filename):
    f = open(filename, 'r')
    error_elements = []
    for i in range(6):
        line = f.readline()
        s = line.split()
    while line:
        if s[0] != 'step' and s[0] != 'number':
            if len(s[1]) != 5:
                s[1] = s[1] + " " + s[2]
            error_elements.append(s[1])
        line = f.readline()
        s = line.split()

    f.close()
    return error_elements

def index_scatter(x,y,z, show = False):
    xs = np.asarray(x)
    ys = np.asarray(y)
    zs = np.asarray(z)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    p1 = ax.scatter(xs,ys,zs, s=20)
    ax.set_xlabel('I')
    ax.set_ylabel('J')
    ax.set_zlabel('K')
    if show == True:
        plt.show()
    plt.savefig('ijkerrors.png')
    plt.clf()
    plt.clf()
    plt.close()
    return 0

def index_scatter_x_slice(x, y, z, ind, show = False):
    # get list of yz coordinates for slice i
    x_mod = []
    y_mod = []
    z_mod = []
    for i in range(len(y)):
        if y[i] == ind:
            x_mod.append(x[i])
            y_mod.append(y[i])
            z_mod.append(z[i])
    xs = np.asarray(x_mod)
    ys = np.asarray(y_mod)
    zs = np.asarray(z_mod)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p1 = ax.scatter(xs,zs, s=20)
    ax.set_xlabel('i')
    ax.set_ylabel('k')
    if show == True:
        plt.show()
    plt.savefig('ijkerror_slice.png')
    plt.clf()
    plt.close()
    return 0

def check_input_material(error_elements):
    # get list of mesh elements and materials
    mesh_elements = []
    mesh_material = {}
    f = open('MESH','r')
    wholefile = f.read()
    lines = wholefile.split('\r\n')
    for line in lines:
        if line != 'CONNE' and line != 'ELEME' and line != 'ina' and line != '':
            s = line.split()
            
            if len(s[0]) == 3:
                s[0] = s[0] + " "  + s[1]
                material = s[2][0:5]
            else: 
                material = s[1][0:5]
            mesh_elements.append(s[0])
            mesh_material[s[0]] = material
    #get list of sand errors
    error_sands = []
    for el in error_elements:
        if mesh_material[el] == 'sands':
            error_sands.append(el)

    return mesh_material, error_sands

def get_error_locs(error_elements, i_map, j_map, k_map):
    x = []
    y = []
    z = []
    for el in error_elements:
        l = list(el)
        il = l[3] + l[4]
        jl = l[1] + l[2]
        kl = l[0]
        
        x.append(i_map[il])
        y.append(j_map[jl])
        z.append(k_map[kl])

    return x, y, z

if __name__ == '__main__':
    # To use this, enable the check_3d_hydro_pressure function in 
    # process_t2_output main.
    # use the command:
    # $ process_t2_output.py > checkijk
    # copy the MESH file from the output directory back to the main directory
    nx = 65
    ny = 119
    nz = 43
    filename = 'checkijk'
    error_elements = read_error_elements(filename)
    mesh_material, error_sands = check_input_material(error_elements)

    i_map, j_map, k_map = ot2.get_index_maps()
    x, y, z = get_error_locs(error_sands, i_map, j_map, k_map)
    index_scatter(x, y, z, show = False)
    index = 2
    index_scatter_x_slice(x,y,z,index, show = False)
