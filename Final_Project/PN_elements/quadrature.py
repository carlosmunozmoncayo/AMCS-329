import numpy as np

def integrate_over_reference(f,deg):
    #Integrate a function f over reference triangle using
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
    return gauss_leg

#def integrate_contribution_any_triangle(triangle,DOF, f,i):
    #Integrate f(x)*phi_i(x) over a triangle
    #i denotes the DOF in triangle corresponding to 
    

    





if __name__=="__main__":
    f =lambda x,y: (x+y)**(-0.5)
    deg=3
    gauss_le = integrate_over_reference(f,deg)
    print(gauss_le)


