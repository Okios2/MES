from math import *
from inspect import signature

class Gauss:
    def __init__(self, n):
        self.n = n #ilość punktów
        self.points = []
        self.weights = []
        if(n == 2):
            self.weights += [1.0,1.0]
            self.points += [-sqrt(1.0/3.0), sqrt(1.0/3.0)]
        if(n == 3):
            self.weights += [5.0/9.0, 8.0/9.0, 5.0/9.0]
            self.points += [-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)]
        if(n == 4): #brak funkcji
            self.weights += [(18.0+sqrt(30.0))/36.0, (18.0+sqrt(30.0))/36.0,  (18.0-sqrt(30.0))/36.0, (18.0-sqrt(30.0))/36.0]
            self.points += [-sqrt((3.0/7.0)-((2.0/7.0)*sqrt(6.0/5.0))), sqrt((3.0/7.0)-((2.0/7.0)*sqrt(6.0/5.0))), -sqrt((3.0/7.0)+((2.0/7.0)*sqrt(6.0/5.0))), sqrt((3.0/7.0)+((2.0/7.0)*sqrt(6.0/5.0)))]
        if(n == 5): #brak funkcji
            self.weights += [128.0/225.0, (322.0+13.0*sqrt(70.0))/900.0, (322.0+13.0*sqrt(70.0))/900.0, (322.0-13.0*sqrt(70.0))/900.0, (322.0-13.0*sqrt(70.0))/900.0]
            self.points += [0, -1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)), 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)), -1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)), 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0))]

    def solve(self, fun):
        self.function = fun
        self.result = 0
        if(len(signature(fun).parameters)==1):
            for i in range(self.n):
                self.result += self.weights[i]*fun(self.points[i])
        else:
            for i in range(self.n):
                for j in range(self.n):
                    self.result += self.weights[i]*self.weights[j]*fun(self.points[i], self.points[j])


    def getResult(self):
        return self.result



def fun(x):
    return 5*x**2+3*x+6

def fun2(x, y):
    return 5*(x**2)*(y**2)+3*x*y+6


test = Gauss(5)
test.solve(fun)
print(test.getResult())

test2 = Gauss(5)
test2.solve(fun2)
print(test2.getResult())





