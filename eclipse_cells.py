#Authors - Karl Bandilla and Evan Leister
class Cell:
    def __init__(self, node1, node2, node3, node4, node5, \
                 node6, node7, node8, i, j, k):
        self.NearLeftBottom = node1
        self.NearRightBottom = node2
        self.FarLeftBottom = node3
        self.FarRightBottom = node4
        self.NearLeftTop = node5
        self.NearRightTop = node6
        self.FarLeftTop = node7
        self.FarRightTop = node8
        self.I = i
        self.J = j
        self.K = k

    def setXPermeability(self, val):
        self.XPermeability = val

    def setZPermeability(self, val):
        self.ZPermeability = val

    def setPorosity(self, val):
        self.Porosity = val

    def getXPermeability(self):
        return self.XPermeability

    def getZPermeability(self):
        return self.ZPermeability

    def getPorosity(self):
        return self.Porosity

    def getCenterX(self):
        return self.centerX

    def getCenterY(self):
        return self.centerY

    def getTopZ(self):
        return self.centerTop

    def getBottomZ(self):
        return self.centerBottom

    def getCorners(self):
        val = (self.NearLeftBottom, self.NearRightBottom,\
                self.FarLeftBottom , self.FarRightBottom,\
                self.NearLeftTop, self.NearRightTop,\
                self.FarLeftTop, self.FarRightTop)
        return val

    def setGeometry(self):
        nlx, nly = self.NearLeftBottom.getXY()
        nrx, nry = self.NearRightBottom.getXY()
        flx, fly = self.FarLeftBottom.getXY()
        frx, fry = self.FarRightBottom.getXY()
        self.centerX = (nlx + nrx + flx + frx) / 4.
        self.centerY = (nly + nry + fly + fry) / 4.
        nlb = self.NearLeftBottom.getZ()
        nrb = self.NearRightBottom.getZ()
        nlt = self.NearLeftTop.getZ()
        nrt = self.NearRightTop.getZ()
        flb = self.FarLeftBottom.getZ()
        frb = self.FarRightBottom.getZ()
        flt = self.FarLeftTop.getZ()
        frt = self.FarRightTop.getZ()
        self.centerBottom = (nlb + nrb +flb + frb) / 4.
        self.centerTop = (nlt + nrt +flt + frt) / 4.

    def checkConsistency(self):
        if self.NearLeftBottom.getXY() != self.NearLeftTop.getXY():
            print "Problem in cell ", self.I, self.J, self.K
        if self.NearRightBottom.getXY() != self.NearRightTop.getXY():
            print "Problem in cell ", self.I, self.J, self.K
        if self.FarLeftBottom.getXY() != self.FarLeftTop.getXY():
            print "Problem in cell ", self.I, self.J, self.K
        if self.FarRightBottom.getXY() != self.FarRightTop.getXY():
            print "Problem in cell ", self.I, self.J, self.K
        if self.NearLeftBottom.getZ() >= self.NearLeftTop.getZ():
            print "Problem in cell ", self.I, self.J, self.K
        if self.NearRightBottom.getZ() >= self.NearRightTop.getZ():
            print "Problem in cell ", self.I, self.J, self.K
        if self.FarLeftBottom.getZ() >= self.FarLeftTop.getZ():
            print "Problem in cell ", self.I, self.J, self.K
        if self.FarRightBottom.getZ() >= self.FarRightTop.getZ():
            print "Problem in cell ", self.I, self.J, self.K

    def getIJK(self):
        return self.I, self.J, self.K

    def writeCellData(self):
        tempstring = str(self.I) + ' ,' + str(self.J) + ' ,' + str(self.K) \
                     + ' ,'
        #near left
        tempx, tempy = self.NearLeftBottom.getXY()
        tempstring += str(tempx) + ', ' + str(tempy) + ', '
        tempstring += str(self.NearLeftBottom.getZ())+ ', '
        tempstring += str(self.NearLeftTop.getZ())+ ', '
        #near right
        tempx, tempy = self.NearRightBottom.getXY()
        tempstring += str(tempx) + ', ' + str(tempy) + ', '
        tempstring += str(self.NearRightBottom.getZ())+ ', '
        tempstring += str(self.NearRightTop.getZ())+ ', '
        #far right
        tempx, tempy = self.FarRightBottom.getXY()
        tempstring += str(tempx) + ', ' + str(tempy) + ', '
        tempstring += str(self.FarRightBottom.getZ())+ ', '
        tempstring += str(self.FarRightTop.getZ())+ ', '
        #far left
        tempx, tempy = self.FarLeftBottom.getXY()
        tempstring += str(tempx) + ', ' + str(tempy) + ', '
        tempstring += str(self.FarLeftBottom.getZ())+ ', '
        tempstring += str(self.FarLeftTop.getZ())+ ', '

        tempstring += str(self.XPermeability) + ', '        
        tempstring += str(self.ZPermeability) + ', '        
        tempstring += str(self.Porosity) + ' \n'

        return tempstring

class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def getXY(self):
        return self.x, self.y

    def getZ(self):
        return self.z
