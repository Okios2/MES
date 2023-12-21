from csv import *
from ElementUniwersalny import *

class GlobalData:
    def __init__(self):
        self.tab = []

    def addGlobalData(self, data):
        self.tab.append(data.replace('\n', ''))

    def showGlobalData(self):
        for i in range(len(self.tab)):
            print(self.tab[i])

    def showSpecificData(self, index):
        print(self.tab[index-1])

    def getValue(self, index):
        value = self.tab[index-1].split()
        value = int(value[1])
        return value

class Node:
    def __init__(self, x:float, y:float, id:int):
        self.x = x
        self.y = y
        self.id = id+1
        self.BC = None
        self.coords = (self.x, self.y)

    def setCoords(self, x:float, y:float):
        self.x = x
        self.y = y
        self.coords = (self.x, self.y)

    def setBC(self):
        self.BC = 1

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getID(self):
        return self.id

    def showCoords(self):
        print(self.coords)

    def showNodeCoords(self):
        print(str(self.id) + " " + str(self.coords))


class Element:
    def __init__(self, a, b, c, d):
        self.ID = [a, b, c, d]
        self.matH = None
        self.matHbc = []
        self.P = []
        self.nodelist = []
        self.nodeBClist = []
        self.matC = []
        for i in range(4):
            self.matHbc.append([])
            self.matC.append([])

    def editElement(self, newElement, index):
        self.ID[index] = newElement

    def showElement(self):
        print(self.ID)
        print("Matrix H")
        print(self.matH)
        print("Matrix Hbc")
        print(self.matHbc)
        print("Vector P")
        print(self.P)
        print("Matrix C")
        print(self.matC)

    def showOneElement(self, index):
        print(self.ID[index])

    def setMatH(self, matH):
        self.matH = matH

    def calcLenL(self, n1:Node, n2:Node):
        return sqrt((n1.x - n2.x)**2 + (n1.y - n2.y)**2)
    #Do zrobienia
    def calcC(self, n, rho, cw, grid):
        eu = ElementUniwersalny(n)
        weights = eu.surface.weights
        result = []

        x = []
        y = []
        for j in range(4):
            x.append(grid.Nodes[self.ID[j] - 1].x)
            y.append(grid.Nodes[self.ID[j] - 1].y)

        det = []
        for i in range(0, n):
            tdet = np.linalg.det(eu.matJ(x, y, i))
            det.append(tdet)

        pointMatC = []
        for i in range(0, n):
            pointMatC.append([])

        for i in range(0, n):
            pointMatC[i].append(rho * cw * np.outer(np.array(eu.N_matrix[i]), np.array(eu.N_matrix[i])) * det[i])

        pointMatC = np.array(pointMatC)
        if (n == 4):
            w1 = 1
            w2 = 1
            self.matC = sum(pointMatC[0] * w1 * w1 + pointMatC[1] * w2 * w1 + pointMatC[2] * w1 * w2 + pointMatC[
                3] * w2 * w2)
        if (n == 9):
            w1 = 5.0 / 9.0
            w2 = 8.0 / 9.0
            w3 = 5.0 / 9.0
            self.matC = sum(pointMatC[0] * w1 * w1 + pointMatC[1] * w1 * w2 + pointMatC[2] * w1 * w3 + pointMatC[3] * w2 * w1 + pointMatC[4] * w2 * w2 + pointMatC[5] * w2 * w3 + pointMatC[6] * w3 * w1 + pointMatC[7] * w3 * w2 + pointMatC[8] * w3 * w3)
        if (n == 16):
            w1 = (18.0 - sqrt(30.0)) / 36.0
            w2 = (18.0 + sqrt(30.0)) / 36.0
            w3 = (18.0 + sqrt(30.0)) / 36.0
            w4 = (18.0 - sqrt(30.0)) / 36.0
            self.matC = sum(pointMatC[0] * w1 * w1 + pointMatC[1] * w1 * w2 + pointMatC[2] * w1 * w3 + pointMatC[
                3] * w1 * w4 + pointMatC[4] * w2 * w1 + pointMatC[5] * w2 * w2 + pointMatC[6] * w2 * w3 + \
                        pointMatC[7] * w2 * w4 + pointMatC[8] * w3 * w1 + pointMatC[9] * w3 * w2 + pointMatC[
                            10] * w3 * w3 + pointMatC[11] * w3 * w4 + pointMatC[12] * w4 * w1 + pointMatC[
                            13] * w4 * w2 + pointMatC[14] * w4 * w3 + pointMatC[15] * w3 * w4)


    def calcHbc(self, n, alfa, T):
        L = 0
        eu = ElementUniwersalny(n)
        points = eu.surface.points
        weights = eu.surface.weights
        result = []
        resultP = []

        for i in self.nodeBClist:
            n1 = None
            n2 = None
            if(i == 0):
                n1 = self.nodelist[0]
                n2 = self.nodelist[1]
            if(i == 1):
                n1 = self.nodelist[2]
                n2 = self.nodelist[1]
            if(i == 2):
                n1 = self.nodelist[2]
                n2 = self.nodelist[3]
            if(i == 3):
                n1 = self.nodelist[3]
                n2 = self.nodelist[0]
            L = self.calcLenL(n1, n2)


            temp = alfa * weights[0]*np.outer(np.array(eu.surface.surfaceTab[i][0]), np.array(eu.surface.surfaceTab[i][0]))
            tempP = alfa * weights[0]*(np.array(eu.surface.surfaceTab[i][0])*T)

            for j in range(1, int(sqrt(n))):
                temp += alfa * weights[j]*np.outer(np.array(eu.surface.surfaceTab[i][j]), np.array(eu.surface.surfaceTab[i][j]))
                tempP += alfa * weights[j] * (np.array(eu.surface.surfaceTab[i][j])*T)

            result.append(temp)
            resultP.append(tempP)
        self.matHbc = sum(result)*L/2.
        self.P = sum(resultP) * L / 2.



class Grid(Node, Element):
    def __init__(self):
        self.Nodes = []
        self.Elements = []

    def addNode(self, node: Node):
        self.Nodes.append(node)

    def addElement(self, element: Element):
        self.Elements.append(element)


    def showNodes(self):
        for i in range(len(self.Nodes)):
            self.Nodes[i].showNodeCoords()

    def showElements(self):
        for i in range(len(self.Elements)):
            print("Element " + str(i+1))
            self.Elements[i].showElement()

    def setMatrixH(self, n, k):

        for i in range(len(self.Elements)):
            eu = ElementUniwersalny(n)
            x = []
            y = []
            for j in range(4):
                x.append(self.Nodes[self.Elements[i].ID[j]-1].x)
                y.append(self.Nodes[self.Elements[i].ID[j] - 1].y)
            self.Elements[i].setMatH(eu.matH(x, y, k))

    def setNodeList(self):
        for i in range(len(self.Elements)):
            for j in range(4):
                self.Elements[i].nodelist.append(self.Nodes[self.Elements[i].ID[j]-1])

        for i in range(len(self.Elements)):
            bc_list = [0, 0, 0, 0]

            for j in range(4):
                if self.Elements[i].nodelist[j].BC == 1:
                    bc_list[j] = 1

            if bc_list[0] == 1 and bc_list[1] == 1:
                self.Elements[i].nodeBClist.append(0)
            if bc_list[1] == 1 and bc_list[2] == 1:
                self.Elements[i].nodeBClist.append(1)
            if bc_list[2] == 1 and bc_list[3] == 1:
                self.Elements[i].nodeBClist.append(2)
            if bc_list[0] == 1 and bc_list[3] == 1:
                self.Elements[i].nodeBClist.append(3)

    def calcAllHbc(self, n, alfa, T):
        for i in self.Elements:
            i.calcHbc(n, alfa, T)

    def calcAllC(self, n, rho, cw):
        for i in self.Elements:
            i.calcC(n, rho, cw, self)


    def setFlagBC(self, BC):
        i = 0
        for node in self.Nodes:
            i+=1
            if i in BC:
                node.BC = 1


class SOE:

    def __init__(self, n, g: Grid):
        self.HG = np.zeros((n,n))
        self.p = np.zeros(n)

        for i in g.Elements:
            FinalH = i.matH + i.matHbc
            FinalH = FinalH[0]
            if(type(i.P) == float):
                i.P = [0,0,0,0]

            for j in range(4):
                self.p[i.ID[j]-1] += i.P[j]
                for k in range(4):
                    self.HG[i.ID[j]-1][i.ID[k]-1] += FinalH[j][k]

        self.result = np.linalg.solve(self.HG, self.p)

    def show(self):
        print(self.HG)
        print()
        print(self.p)
        print()
        print(self.result)

    def fT(self, n, g: Grid, data: GlobalData, T):
        self.HG = np.zeros((n,n))
        self.p = np.zeros(n)
        for i in g.Elements:
            FinalH = (i.matH + i.matHbc) + (i.matC / data.getValue(2))
            for j in range(4):
                element_t0 = []
                for id in i.ID:
                    element_t0.append(T[id - 1])
                self.p[i.ID[j] - 1] += i.P[j] + np.dot(i.matC / data.getValue(2), element_t0)[j]
                for k in range(4):
                    self.HG[i.ID[j] - 1][i.ID[k] - 1] += FinalH[0][j][k]
        self.result = np.linalg.solve(self.HG, self.p)

    def toVTK(self, data: GlobalData, g: Grid, name):
        nodes = data.getValue(9)
        elements = data.getValue(10)
        t = data.getValue(1)
        temp = np.full(nodes, data.getValue(6))
        iter = int((t/data.getValue(2))+1)

        for i in range(iter):
            f = open(f'vtk/{nameSchema}{iter}.vtk', 'w')
            content = f'''# vtk DataFile Version 2.0
            Unstructured Grid Example
            ASCII
            DATASET UNSTRUCTURED_GRID

            POINTS {nodes} float\n'''

            for i in range(nodes):
                content += f'{grid.Nodes[i].x} {grid.Nodes[i].y} 0.0\n'

            content += f'\nCELLS {nElements} {nElements * 5}\n'

            for i in range(elements):
                content += f'4 {g.Elements[i].ID[0]-1} {g.Elements[i].ID[1]-1} {g.Elements[i].ID[2]-1} {g.Elements[i].ID[3]-1}\n'

            content += f'\nCELL_TYPES {nElements}\n'

            for i in range(elements):
                content += f'9\n'

            content += f'''\nPOINT_DATA {nodes}
            SCALARS Temp float 1
            LOOKUP_TABLE default\n'''

            for i in range(nodes):
                content += f'{temp[i]}\n'

            temp = self.tempSystem(g, data.getValue(2), temp)

            f.write(content)













def read(filePath: str):
    file = open(filePath, "r")


#otwracie pliku i inicjalizacja obiektów GlobalData oraz Siatki(Grid)
# Test2_4_4_MixGrid.txt
# Test1_4_4.txt
file = open("Test2_4_4_MixGrid.txt", "r")
Set1 = GlobalData()
Grid1 = Grid()

line = file.readline()
while(line!="*Node\n"):
    Set1.addGlobalData(line)
    line = file.readline()

line = file.readline()
nodeList = []
while(line!="*Element, type=DC2D4\n"):
    coords = line.split(",")
    nodeList.append(coords)
    line = file.readline()

for i in range(len(nodeList)):
    nodei = Node(float(nodeList[i][1]), float(nodeList[i][2]), i)
    Grid1.addNode(nodei)


elementList = []
line = file.readline()
while(line!="*BC\n"):
    elementy = line.split(",")
    elementList.append(elementy)
    line = file.readline()

for i in range(len(elementList)):
    elementi = Element(int(elementList[i][1]), int(elementList[i][2]), int(elementList[i][3]), int(elementList[i][4]))
    Grid1.addElement(elementi)

line = file.readline()
BC = line.split(",")
for i in range(len(BC)):
    BC[i] = int(BC[i])


Set1.showGlobalData()
print()
Grid1.setFlagBC(BC)
Grid1.setNodeList()
Grid1.setMatrixH(9, 25)
Grid1.calcAllHbc(9, Set1.getValue(4), Set1.getValue(5))
Grid1.calcAllC(9, Set1.getValue(7), Set1.getValue(8))
Grid1.showNodes()
print()
Grid1.showElements()
print()
print("BC List")
print(BC)
s = SOE(16, Grid1)
#s.show()
s.fT(16,Grid1, Set1, np.full(16, Set1.getValue(6)))
s.show()
