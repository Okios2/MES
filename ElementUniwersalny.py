from math import *
import numpy as np

class Surface:

    def __init__(self, n):
        self.surfaceTab = []
        self.n = int(sqrt(n))
        for i in range(4):
            self.surfaceTab.append([])

        self.points = []
        self.weights = []
        if (self.n == 2):
            self.weights += [1.0, 1.0]
            self.points += [-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)]
        if (self.n == 3):
            self.weights += [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
            self.points += [-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)]
        if (self.n == 4):  # brak funkcji
            self.weights += [(18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0,
                             (18.0 - sqrt(30.0)) / 36.0]
            self.points += [-sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                            sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                            -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                            sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0)))]

        for i in range(4):
            if (i == 0):
                for j in range(self.n):
                    self.surfaceTab[0].append([])
                    self.surfaceTab[0][j].append(self.N1(self.points[j], -1.0))
                    self.surfaceTab[0][j].append(self.N2(self.points[j], -1.0))
                    self.surfaceTab[0][j].append(0.0)
                    self.surfaceTab[0][j].append(0.0)
            elif (i == 1):
                for j in range(self.n):
                    self.surfaceTab[1].append([])
                    self.surfaceTab[1][j].append(0.0)
                    self.surfaceTab[1][j].append(self.N2(1.0, self.points[j]))
                    self.surfaceTab[1][j].append(self.N3(1.0, self.points[j]))
                    self.surfaceTab[1][j].append(0.0)
            elif (i == 2):
                for j in range(self.n):
                    self.surfaceTab[2].append([])
                    self.surfaceTab[2][j].append(0.0)
                    self.surfaceTab[2][j].append(0.0)
                    self.surfaceTab[2][j].append(self.N3(self.points[j], 1.0))
                    self.surfaceTab[2][j].append(self.N4(self.points[j], 1.0))
            elif (i == 3):
                for j in range(self.n):
                    self.surfaceTab[3].append([])
                    self.surfaceTab[3][j].append(self.N1(-1.0, self.points[j]))
                    self.surfaceTab[3][j].append(0.0)
                    self.surfaceTab[3][j].append(0.0)
                    self.surfaceTab[3][j].append(self.N4(-1.0, self.points[j]))


    def N1(self, ksi, eta):
        return 1. / 4. * (1 - ksi)*(1 - eta)

    def N2(self, ksi, eta):
        return 1. / 4. * (1 + ksi)*(1 - eta)

    def N3(self, ksi, eta):
        return 1. / 4. * (1 + ksi)*(1 + eta)

    def N4(self, ksi, eta):
        return 1. / 4. * (1 - ksi)*(1 + eta)



class ElementUniwersalny:

    def __init__(self, l_matrix):
        self.matrixKsi = []
        self.matrixEta = []
        self.matrixX = []
        self.matrixY = []
        self.matrixj = []
        self.matrixH = []
        self.N_matrix = np.zeros((l_matrix, 4))
        self.l_matrix = int(l_matrix)
        self.surface = Surface(self.l_matrix)


        self.elementStruckt = []
        for i in range(4):
            self.elementStruckt.append([])
            self.elementStruckt.append([])
            self.matrixH.append([])



        for i in range(0, l_matrix):
            self.matrixKsi.append([])
            self.matrixEta.append([])
            self.matrixX.append([])
            self.matrixY.append([])
            self.matrixj.append([])



    #4x4
        if(l_matrix == 4):
            iter = 0
            for i in [-1. / sqrt(3), -1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3)]:
                self.matrixKsi[iter].append(self.fun1(i))
                self.matrixKsi[iter].append(self.fun2(i))
                self.matrixKsi[iter].append(self.fun3(i))
                self.matrixKsi[iter].append(self.fun4(i))
                iter+=1
            iter = 0
            for i in [-1. / sqrt(3), 1. / sqrt(3), -1. / sqrt(3), 1. / sqrt(3)]:
                self.matrixEta[iter].append(self.fun1(i))
                self.matrixEta[iter].append(self.fun4(i))
                self.matrixEta[iter].append(self.fun3(i))
                self.matrixEta[iter].append(self.fun2(i))
                iter+=1

            points = [-1. / sqrt(3),  1. / sqrt(3)]
            Ni = [self.N1, self.N2, self.N3, self.N4]
            for i in range(4):
                ksi = i%2
                eta = int(i/2)
                for j in range(4):
                    self.N_matrix[i][j] = Ni[j](points[ksi],points[eta])


    #4x9
        if(l_matrix == 9):
            iter = 0
            for i in [-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)]:
                self.matrixKsi[iter].append(self.fun1(i))
                self.matrixKsi[iter].append(self.fun2(i))
                self.matrixKsi[iter].append(self.fun3(i))
                self.matrixKsi[iter].append(self.fun4(i))
                iter+=1
            iter = 0
            for i in [-sqrt(3.0/5.0), -sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, 0.0, 0.0, sqrt(3.0/5.0), sqrt(3.0/5.0), sqrt(3.0/5.0)]:
                self.matrixEta[iter].append(self.fun1(i))
                self.matrixEta[iter].append(self.fun4(i))
                self.matrixEta[iter].append(self.fun3(i))
                self.matrixEta[iter].append(self.fun2(i))
                iter+=1

            points = [-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)]
            Ni = [self.N1, self.N2, self.N3, self.N4]
            for i in range(9):
                ksi = i%9
                eta = int(i/9)
                for j in range(4):
                    self.N_matrix[i][j] = Ni[j](points[ksi],points[eta])
    #4x16
        if (l_matrix == 16):
            iter = 0
            for i in [-sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0)))]:
                self.matrixKsi[iter].append(self.fun1(i))
                self.matrixKsi[iter].append(self.fun2(i))
                self.matrixKsi[iter].append(self.fun3(i))
                self.matrixKsi[iter].append(self.fun4(i))
                iter += 1
            iter = 0
            for i in [-sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),

                      sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0)))]:
                self.matrixEta[iter].append(self.fun1(i))
                self.matrixEta[iter].append(self.fun4(i))
                self.matrixEta[iter].append(self.fun3(i))
                self.matrixEta[iter].append(self.fun2(i))
                iter += 1

            points = [-sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), -sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      -sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))),

                      sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) - ((2.0 / 7.0) * sqrt(6.0 / 5.0))),
                      sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0))), sqrt((3.0 / 7.0) + ((2.0 / 7.0) * sqrt(6.0 / 5.0)))]
            Ni = [self.N1, self.N2, self.N3, self.N4]
            for i in range(16):
                ksi = i%16
                eta = int(i/16)
                for j in range(4):
                    self.N_matrix[i][j] = Ni[j](points[ksi],points[eta])



    #pochodne funkcji kształtu
    def fun1(self, x):
        return -1. / 4. * (1 - x)

    def fun2(self, x):
        return 1. / 4. * (1 - x)

    def fun3(self, x):
        return 1. / 4. * (1 + x)

    def fun4(self, x):
        return -1. / 4. * (x + 1)
    #funkcje kształtu
    def N1(self, x, y):
        return 1. / 4. * (1 - x) * (1 - y)

    def N2(self, x, y):
        return 1. / 4. * (1 + x) * (1 - y)

    def N3(self, x, y):
        return 1. / 4. * (1 + x) * (1 + y)

    def N4(self, x, y):
        return 1. / 4. * (1 - x) * (1 + y)

#Obliczanie macierzy H
    def dxdKsi(self, tabx, point):
        wynik = 0
        # print("iterKsi")
        for i in range(0, len(tabx)):
            wynik += tabx[i]*self.matrixKsi[point][i]
        #     print(tabx[i])
        #     print(wynik, "\n")
        # print("Suma = " , wynik)
        return wynik

    def dxdEta(self, tabx, point):
        wynik = 0
        # print("IterEta")
        for i in range(0, len(tabx)):
            wynik += tabx[i]*self.matrixEta[point][i]
        #     print(tabx[i])
        #     print(wynik, "\n")
        # print("Suma = " , wynik)
        return wynik

    def matJ(self, tabx, taby, point):
        matrixj = []
        for i in range(2):
            matrixj.append([])
        matrixj[0].append(self.dxdKsi(tabx, point))
        matrixj[0].append(self.dxdKsi(taby, point))
        matrixj[1].append(self.dxdEta(tabx, point))
        matrixj[1].append(self.dxdEta(taby, point))
        return matrixj

    def matH(self, tabx, taby, k):
        matJ = []
        #tablica macierzy jakobianów
        for i in range(0, self.l_matrix):
            matJ.append(self.matJ(tabx, taby, i))

        #tablica wyznaczników jakobianów
        tabDetJ = []
        # tablica 1/jakobian
        tabDetRevJ = []
        # macierz [[y/eta -x/ksi][-y/eta x/ksi]]
        tempMat = []
        for i in range(0, len(matJ)):
            tabDetJ.append(np.linalg.det(matJ[i]))
            tabDetRevJ.append(1./(tabDetJ[i]))
            tempMat = [[matJ[i][1][1], -matJ[i][0][1]], [-matJ[i][1][0], matJ[i][0][0]]]
            self.matrixj[i].append(np.array(tempMat)*tabDetRevJ[i])


        for i in range(0, self.l_matrix):
            for j in range(4):
                self.matrixX[i].append(self.matrixj[i][0][0][0]*self.matrixKsi[i][j]+self.matrixj[i][0][0][1]*self.matrixEta[i][j])
                self.matrixY[i].append(self.matrixj[i][0][1][0] * self.matrixKsi[i][j] + self.matrixj[i][0][1][1] * self.matrixEta[i][j])


        x = np.array(self.matrixX)
        y = np.array(self.matrixY)


        pointMatH = []
        for i in range(0, self.l_matrix):
            pointMatH.append([])

        for i in range(0, self.l_matrix):
            # transponowane wiersze macierzy X i Y
            tx = np.array([x[i]])
            ty = np.array([y[i]])
            pointMatH[i].append(k*((self.matrixX[i]*tx.T)+(self.matrixY[i]*ty.T))*tabDetJ[i])

        pointMatH = np.array(pointMatH)

        if(self.l_matrix == 4):
            w1 = 1
            w2 = 1
            self.matrixH = pointMatH[0]*w1*w1+pointMatH[1]*w2*w1+pointMatH[2]*w1*w2+pointMatH[3]*w2*w2
        if (self.l_matrix == 9):
            w1 = 5.0/9.0
            w2 = 8.0/9.0
            w3 = 5.0/9.0
            self.matrixH = pointMatH[0]*w1*w1+pointMatH[1]*w1*w2+pointMatH[2]*w1*w3+pointMatH[3]*w2*w1+pointMatH[4]*w2*w2+pointMatH[5]*w2*w3+pointMatH[6]*w3*w1+pointMatH[7]*w3*w2+pointMatH[8]*w3*w3
        if (self.l_matrix == 16):
            w1 = (18.0-sqrt(30.0))/36.0
            w2 = (18.0+sqrt(30.0))/36.0
            w3 = (18.0+sqrt(30.0))/36.0
            w4 = (18.0-sqrt(30.0))/36.0
            self.matrixH = pointMatH[0]*w1*w1+pointMatH[1]*w1*w2+pointMatH[2]*w1*w3+pointMatH[3]*w1*w4+pointMatH[4]*w2*w1+pointMatH[5]*w2*w2+pointMatH[6]*w2*w3+pointMatH[7]*w2*w4+pointMatH[8]*w3*w1+pointMatH[9]*w3*w2+pointMatH[10]*w3*w3+pointMatH[11]*w3*w4+pointMatH[12]*w4*w1+pointMatH[13]*w4*w2+pointMatH[14]*w4*w3+pointMatH[15]*w3*w4

        return self.matrixH


    def showKsi(self):
        print("Funkcje kształtu dN/dKsi")
        for j in range(4):
            print("dN" + str(j+1) + "/dKsi")
            for i in range(len(self.matrixKsi)):
                print(self.matrixKsi[i][j])
        print("")

    def showEta(self):
        print("Funkcje kształtu dN/dEta")
        for j in range(4):
            print("dN" + str(j+1) + "/dEta")
            for i in range(len(self.matrixEta)):
                print(self.matrixEta[i][j])
        print("")


matrix1 = ElementUniwersalny(4)
# print(matrix1.matrixKsi)
# print(matrix1.matrixEta)
matrix2 = ElementUniwersalny(9)
# print(matrix2.matrixKsi)
# print(matrix2.matrixEta)
matrix3 = ElementUniwersalny(16)

# matrix1.showKsi()
# matrix1.showEta()

#[x[0, 0.025, 0.025, 0]y[0, 0, 0.025, 0.025]]

# dla 2 pkt
# print("2 pkt:")
# k = matrix1.matH([0, 0.025, 0.025, 0], [0, 0, 0.025, 0.025], 30)
# print(matrix1.matrixH[0])
#
# # dla 3 pkt
# print("3 pkt:")
# k = matrix2.matH([0, 0.025, 0.025, 0], [0, 0, 0.025, 0.025], 30)
# print(matrix2.matrixH[0])

# print("4 pkt:")
# k = matrix3.matH([0, 0.025, 0.025, 0], [0, 0, 0.025, 0.025], 30)
# print(matrix3.matrixH[0])

#[x[0, 0.025, 0.025, 0]y[0, 0, 0.025, 0.025]]

#print(matrix1.surface.surfaceTab)


#TESTY

# node0 = [0,0.0376100652, -0.0573899336]
# node1 = [1,0.0, -0.0496918149]
# node2 = [2,0.0, -0.0949999988]
# node3 = [3,0.0453081839, -0.0949999988]
#
# list = [node0, node1, node2, node3]
#
#list2 = [0.0376100652, 0.0, 0.0, 0.0453081839]
#list3 = [-0.0573899336, -0.0496918149, -0.0949999988, -0.0949999988]


# for i in range(0,4):
#     print(i)
#     print(matrix1.matJ(list2, list3, i))

# print(list2)
# print(list3)
# print(matrix1.matJ(list2, list3, 0))

#print(matrix1.matH(list2, list3, 25))

print(matrix1.matrixj)


