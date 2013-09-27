#Authors - Karl Bandilla and Evan Leister
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import eclipse_cells as ec

def read_eclipse():
    print "reading Eclipse Files for Sleipner"
    # Sets up the list of cell objects. 
    inp = open('M9X1.grdecl', 'r')
    newline = inp.readline()
    while newline != 'SPECGRID\n':
        newline = inp.readline()
    newline = inp.readline()
    parameters = newline.split(' ')
    Nx = int(parameters[0])
    Ny = int(parameters[1])
    Nz = int(parameters[2])
    while newline != 'COORD\n':
        newline = inp.readline()
    counter = 0
    rows = []
    temprow = []
    while newline:
        newline = inp.readline()
        if newline == '\n':
            rows.append(temprow)
            temprow = []
        elif newline == '/\n':
            break
        else:
            parameters = newline.split(' ')
            temprow.append([float(parameters[0]), float(parameters[1])])

    while newline != 'ZCORN\n':
        newline = inp.readline()

    tempplanes = []
    templine = []
    temprows = []
    while newline:
        newline = inp.readline()
        if newline == 'ACTNUM\n':
            break
        parameters = newline.split(' ')
        for par in parameters[:-1]:
            templine.append(float(par))
        if len(templine) == 2 * Nx:
            counter = 0
            templist = []
            for i in range(0, Nx):
                templist.append(templine[counter])
                counter += 2
            templist.append(templine[-1])
            templine = []
            temprows.append(templist)
            if len(temprows) == Ny*2:
                counter = 0
                templist = []

                for i in range(0, Ny):
                    templist.append(temprows[counter])
                    counter += 2
                templist.append(temprows[-1])

                tempplanes.append(templist)
                temprows = []
    inp.close()
    planes = []
    counter = 0
    for i in range(0, Nz):
        planes.append(tempplanes[counter])
        counter += 2
    planes.append(tempplanes[-1])
    locations = []
    for row in rows:
        for val in row:
            locations.append(val)
    newplanes = []
    for plane in planes:
        tempplane = []
        for row in plane:
            for val in row:
                tempplane.append(val)
        newplanes.append(tempplane)
    Nodes = []
    for plane in newplanes:
        counter = 0
        for val in plane:
            Nodes.append(ec.Node(locations[counter][0], locations[counter][1], \
                    val))
            counter += 1

    counter = 0
    Cells = []
    for k in range(0, Nz):
        for j in range(0, Ny):
            for i in range(0, Nx):
                nearCornerID = i + j * (Nx+1) + k * (Nx+1) * (Ny+1)
                farCornerID = nearCornerID + Nx +1
                topNearCornerID = nearCornerID + (Nx+1) * (Ny+1)
                topFarCornerID = farCornerID + (Nx+1) * (Ny+1)
                Cells.append(ec.Cell(Nodes[nearCornerID], Nodes[nearCornerID+1], \
                        Nodes[farCornerID], Nodes[farCornerID+1], \
                        Nodes[topNearCornerID], Nodes[topNearCornerID+1], \
                        Nodes[topFarCornerID], Nodes[topFarCornerID+1], \
                        i+1, j+1, k+1))
                Cells[-1].setGeometry()


    #for cel in Cells:
    #    cel.checkConsistency()

    inp = open('M9X1_perm_X_mD_.inc', 'r')
    newline = inp.readline()
    newline = inp.readline()
    newline = inp.readline()
    parameters = newline.split(' ')
    tempPermX = []
    while newline:
        if newline == '/\n':
            break
        parameters = newline.split(' ')
        for par in parameters:
            if par != '\n':
                tempPermX.append(float(par))
        newline = inp.readline()
    inp.close()

    inp = open('M9X1_perm_Z_mD_.inc', 'r')
    newline = inp.readline()
    newline = inp.readline()
    newline = inp.readline()
    parameters = newline.split(' ')
    tempPermZ = []
    while newline:
        if newline == '/\n':
            break
        parameters = newline.split(' ')
        for par in parameters:
            if par != '\n':
                tempPermZ.append(float(par))
        newline = inp.readline()
    inp.close()

    inp = open('M9X1_poro___.inc', 'r')
    newline = inp.readline()
    newline = inp.readline()
    newline = inp.readline()
    parameters = newline.split(' ')
    tempPoro = []
    while newline:
        if newline == '/\n':
            break
        parameters = newline.split(' ')
        for par in parameters:
            if par != '\n':
                tempPoro.append(float(par))
        newline = inp.readline()
    inp.close()

    counter = 0
    for cel in Cells:
        cel.setXPermeability(tempPermX[counter])
        cel.setZPermeability(tempPermZ[counter])
        cel.setPorosity(tempPoro[counter])
        counter += 1


    # EVERYTHING ABOVE THIS IS FOR CELL setup. 

    #For plots of numerical layers

    #this set which layer gets displayed
    mylayer = 30
    singleLayer = []
    for cel in Cells:
        tempi, tempj, tempk = cel.getIJK()
        if tempk == mylayer:
            singleLayer.append(cel)


    vals = []
    locX = []
    locY = []
    tempvals = []
    counter = 0
    vals = []
    for cel in singleLayer:
        #change this to pick what value to contour or surface
        tempvals.append(-cel.getTopZ())
        counter += 1
        tempi, tempj, tempk = cel.getIJK()
        if counter == Nx:
            vals.append(tempvals)
            tempvals = []
            counter = 0
        if tempj == 1:
            locX.append(cel.getCenterX())
        if tempi == 1:
            locY.append(cel.getCenterY())
    xs = np.array(locX)
    ys = np.array(locY)
    X, Y = np.meshgrid(xs, ys)
    Z = np.asarray(vals)

    return Cells, Nx, Ny, Nz
