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
            gauss_leg += wxi[i]*weta[j]*f(np.array(x[i],y[i][j]))*(1-xi[i])/8.
    return gauss_leg

#def integrate_contribution_any_triangle(triangle,DOF, f,i):
    #Integrate f(x)*phi_i(x) over a triangle
    #i denotes the DOF in triangle corresponding to 
    
def integrate_f_phi(shape_function,affine_mapping,f,degree):
    def compose_f_phi_affine(x_ref):
            x = affine_mapping.A@x_ref + affine_mapping.b
        return f(x)*shape.phi(x)
    return affine_mapping.detA*integrate_over_reference(compose_f_phi_affine,degree)

def integrate_gradphiI_dot_gradphiJ(shapef_i,shapef_j,affine_mapping,degree):
    def compose_grad_phi_affine(x_ref):
        x = affine_mapping.A@x_ref + affine_mapping.b
        return np.dot(shapef_i.grad_phi(x),shapef_j.grad_phi(x))
    return affine_mapping.detA*integrate_over_reference(compose_grad_phi_affine,degree)


    





if __name__=="__main__":
    f =lambda x,y: (x+y)**(-0.5)
    deg=3
    gauss_le = integrate_over_reference(f,deg)
    print(gauss_le)


