import numpy as np
from scipy import integrate


def integrate_f_over_reference(f,deg):
    #Integrate a function over reference triangle using
    #Gauss-Legendre quadrature of degree deg.
    xi, wxi = np.polynomial.legendre.leggauss(deg)
    eta, weta = np.polynomial.legendre.leggauss(deg)
    #Following Rathod et al. 

    x = (1+xi)/2.
    y = [[(1-xi[i])*(1+eta[j])/4. for j in range(len(eta))] for i in range(len(xi))]
    
    gauss_leg = 0.
    for i in range(len(xi)):
        for j in range(len(eta)):
            gauss_leg += wxi[i]*weta[j]*f(x[i],y[i][j])*(1-xi[i])/8.
            #np.sum(w * f(s)) * 0.5*(b - a)

    return gauss_leg




if __name__=="__main__":
    f =lambda x,y: (x+y)**(-0.5)
    deg=30
    gauss_le = integrate_f_over_reference(f,deg)
    print(gauss_le)


