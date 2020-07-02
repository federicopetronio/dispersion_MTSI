import numpy as np
import scipy as sc
from tqdm import tqdm_notebook as tqdm

import matplotlib.pyplot as plt
import multiprocessing as mp

class DirectSolver():

    def __init__(self, xmin, xmax, ymin, ymax, costFunc, method="Nelder-Mead"):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.costFunc = costFunc
        self.domain = (xmin, ymin, xmax, ymax)
        self.x0 = None
        self.method = method

        if self.method in ["Powell", "Nelder-Mead"]:
            self.jac=None
            self.optionslow = {"fatol": 1e-12, "maxiter":800}
            self.optionshigh = {"fatol": 1e-16, "maxiter":1000}
        else:
            self.jac='2-point'
            self.optionslow = {}
            self.optionshigh = {}


    def setks(self, kx, ky, kz):
        self.kx = kx
        self.ky = ky
        self.kz = kz

    def funct(self, xs):
        x = np.abs(self.costFunc(xs[0] + 1j*xs[1], self.kx, self.ky, self.kz))
        if x is np.nan:
            return 10^300
        return x



    def nelder_mead(self):
        x0 = self.get_x0()
        all_x_i = [x0[0]]
        all_y_i = [x0[1]]
        all_f_i = [self.funct(x0)]

        def store(X):
            x, y = X
            all_x_i.append(x)
            all_y_i.append(y)
            all_f_i.append(self.funct(X))

        sc.optimize.minimize(self.funct, x0, jac=self.jac, method=self.method , callback=store, options=self.optionslow)
        r = self.funct([all_x_i[-1], all_y_i[-1]])
        if r > 0.001:
            all_x_i = [x0[0]]
            all_y_i = [x0[1]]
            all_f_i = [self.funct(x0)]
            sc.optimize.minimize(self.funct, x0, jac=self.jac, method=self.method , callback=store, options=self.optionshigh)

        return all_x_i, all_y_i, all_f_i

    def set_x0(self, x0):
        self.x0 = x0

    def get_x0(self):
        if self.x0 is None:
            return np.array([(self.xmax+self.xmin)/2, (self.ymin+self.ymax)/2])
        else:
            return self.x0

    def solve(self):
        """We use nelder_mead because I am not sure of the Jacobian of the matrix"""
        all_x_i, all_y_i, all_f_i = self.nelder_mead()
        x = [all_x_i[-1], all_y_i[-1]]
        self.set_x0(x)
        return x




class SolveKy:
    def __init__(self, costFunc, kx, kz, method, x0=[0.1,0.0], constfuncRef=None, wrfunct=None):
        self.DS = DirectSolver(0.0001,1,0,1,costFunc = costFunc, method=method)
        self.constfuncRef = constfuncRef
        self.DSref = DirectSolver(0.0001,1,0,1,costFunc=constfuncRef, method=method)
        self.DS.set_x0(x0)
        self.DSref.set_x0(x0)
        self.kx = kx
        self.kz = kz
        self.wrfunct = wrfunct

    def __call__(self, ky):

        if self.constfuncRef is not None:
            self.DSref.setks(self.kx, ky, self.kz)
            x0 = self.DSref.solve()
            self.DS.set_x0(x0)

        else:
            if self.wrfunct is not None:
                self.DS.set_x0([self.wrfunct(ky), 0.01])
            else:
                if ky > 0.5:
                    self.DS.set_x0([ky/2, 0.1])

        self.DS.setks(self.kx, ky, self.kz)
        x = self.DS.solve()
        return x


def solvekys(costFunc, kx, kz, kymin=0.01, kymax = 2, Nkys= 100, method="Nelder-Mead", plot=False, constfuncRef=None, parall=True, wrfunct=None):
    DS = DirectSolver(0.0001,1,0,1,costFunc = costFunc, method=method)
    DS.set_x0([0,0.0])
    kys = np.linspace(kymin, kymax, Nkys)
    xs = np.zeros((len(kys), 2))

    if parall:
        tmpObj = SolveKy(costFunc, kx, kz, method)
        x0 = tmpObj((kymin))

        with mp.Pool() as pool:
            xs = pool.map(SolveKy(costFunc, kx, kz, method, x0=x0, constfuncRef=constfuncRef, wrfunct=wrfunct), kys)
        xs = np.array(xs)
    else:
        DS = DirectSolver(0.0001,1,0,1,costFunc = costFunc, method=method)

        if wrfunct is not None:
            DS.set_x0([1.5*wrfunct(kymin), 0.1])
        else:
            if kymax > 0.5:
                DS.set_x0([kymin/2, 0.1])

        for i,ky in enumerate(kys):
            DS.setks(kx, ky, kz)
            x0 = DS.solve()
            xs[i,:] = x0
            DS.set_x0(x0)


    if plot:
        plt.figure()
        plt.plot(kys, xs[:,0])
        plt.plot(kys, xs[:,1])

    else:
        return kys, xs
