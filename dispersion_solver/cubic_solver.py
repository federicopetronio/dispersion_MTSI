import numpy as np

"""TO USE WITH THE ANALITIC APPROX DISPERSION"""
def cubic_solver(a,b,c,d):

    solution = np.zeros(3, dtype=complex)

    delta_0 = b**2 - 3*a*c
    delta_1 = 2*b**3 - 9*a*b*c +27*a**2*d
    C = ((delta_1-(delta_1**2-4*delta_0**3)**0.5)/2)**(1/3)

    print(delta_0,delta_1,C)

    epss = (-1 + (-3)**0.5)/2

    for ind in [0,1,2]:
        solution[ind] = -1/(3*a) * (b + epss**ind * C + delta_0/(epss**ind * C))
    return solution
