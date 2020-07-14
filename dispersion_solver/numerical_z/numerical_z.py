# -*- coding: utf-8 -*-

"""Main module."""

import numpy as np
import numba
from scipy.integrate import  quad, simps

class NumericZ:
    """A class that can store the information from the numerical calculation of Z
    Then useing the __call__  methode to compute the Z function.
    """

    def __init__(self, F_function, N=64):
        self.F = F_function
        self.N = N
        self.delt = 1

        self.L = np.sqrt(self.N / np.sqrt(2))
        self.compute_a_real()
        self.compute_a_complex()

        self.imag_tol = 100000

    def setimagetol(self, value):
        self.imag_tol = value

    def compute_a_real(self):
        """Compute the a_n coefficients for the real axis"""
        N = self.N
        n = np.arange(-N, N)
        v = self.L * np.tan(np.pi * (n + 1 / 2) / (2 * N))
        FF = self.F(v)
        a = np.fft.fft(np.fft.fftshift(FF * (self.L - 1.j * v)))
        a = np.exp(- 1.j * n * np.pi / (2 * N)) * np.fft.fftshift(a)
        a = np.flipud(1.j * np.sign(n + 1 / 2) * a) / (2 * N)

        self.areal = a

    def compute_a_complex(self):
        """Compute the a_n coefficients for the complex domain"""
        N = self.N  # Two times more points for the complex domain

        n = np.arange(-N, N)
        v = self.L * np.tan(np.pi * (n) / (2 * N))
        FF = self.F(v) * (self.L ** 2 + v ** 2)

        a = np.fft.fft(np.fft.fftshift(FF)) / (2 * N)
        a0 = a[0]
        a = np.flipud(a[1:N])

        self.aimag0 = a0
        self.aimag = a

    def hilb1(self, zf):
        """Hilbert transfor on the real axis"""
        t = (self.L + 1.j * zf) / (self.L - 1.j * zf)
        #  Real axis
        h = polyval(self.areal, t) / (t ** self.N * (self.L - 1.j * zf))
        return h

    def hilb2(self, zf):

        t = np.array((self.L + 1.j * zf) / (self.L - 1.j * zf))
        if t.size > 0:
            p = polyval(self.aimag, t)
            h2 = 1.j * (2 * p / (self.L - 1.j * zf) ** 2 + (self.aimag0 / self.L) / (self.L - 1.j * zf))
        else:
            h2 = np.array([])
        return h2

    def hilb(self, w):
        """hilb"""
        zvalues = np.zeros(w.shape, dtype=np.complex)

        #  Real axis
        realindex = np.nonzero(np.imag(w) == 0)

        if len(realindex) > 0:
            wreal = w[realindex]
            h = self.hilb1(wreal)
            FFreal = self.F(wreal)
            zvalues[realindex] = h + 1.j*FFreal

        # Complex domain
        compindex = np.nonzero(np.imag(w) != 0)

        if len(compindex) > 0:
            wcomp = w[compindex]
            zf = w[compindex]

            negativimag = np.nonzero(zf.imag < 0)
            positivimag = np.nonzero(zf.imag != 0)

            #zf[negativimag] = np.conj(zf[negativimag])

            h2 = self.hilb2(zf)
            #FFcomp = self.F(wcomp[negativimag])

            try:
                #zvalues[negativimag] = np.conj(h2[negativimag]) #+ self.delt*2.j*FFcomp
                zvalues[positivimag] = h2[positivimag]
            except TypeError:
                print(zvalues)
                print(h2)
                print(negativimag)


        return np.pi * zvalues

    def __call__(self, w):

        if not isinstance(w, np.ndarray):
            w = np.array([w])
        try:
            z = np.zeros(len(w), dtype=complex)
        except TypeError:
            print(w, type(w))
            print(w.shape)

        #indexreal = np.nonzero(np.abs(w.imag) <= self.imag_tol)
        # indexreal = np.arange(w.shape[0])

        #if len(indexreal) > 0:
        z = self.hilb(w)

        indexcomp = np.argwhere(np.abs(w.imag) > self.imag_tol)

        # indexcomp = []
        def integrand(x, zi):
            return self.F(x) / (x - zi)

        if len(indexcomp) > 0:
            try:
                for i, zi in enumerate(w[indexcomp]):
                    tmp_v = complex_quadrature(integrand, -10, 10, zi=zi)
                    # tmp_v = complex_quadrature(integrand, -10,10, args=(zi,), epsabs=1.49e-04, epsrel=1.49e-04,
                    #                                    full_output = 1)
                    z[indexcomp[i]] = tmp_v

            except IndexError:
                print(i, indexcomp[i], zi)

        return z


@numba.jit(nopython=True)
def polyval(p, x):
    """Faster (?) than numpy.polyval
    """

    y = np.zeros_like(x)
    for i in range(len(p)):
        y = y * x + p[i]
    return y




def complex_quadrature_quad(func, a, b, **kwargs):
    def real_func(x, *args):
        return func(x, *args).real

    def imag_func(x, *args):
        return func(x, *args).imag

    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return real_integral[0] + 1j * imag_integral[0]


def complex_quadrature(func, a, b, **kwargs):
    v_axis = np.linspace(-10, 10, 512)
    f_points = func(v_axis, **kwargs)
    integral = simps(f_points, x=v_axis)

    return integral
