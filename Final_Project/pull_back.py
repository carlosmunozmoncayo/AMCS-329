import numpy as np
from scipy.spatial import Delaunay


#Points of reference triangle:
#xr1=(0,1)
#xr2=(0,0)
#xr3=(1,0)

#Assuming an affine push forward function, i.e.,
#x = Ft(xri) = A xri+b

def p_back_forward(triangle,points):
    #Returns important quantities given a triangle and list of points
    #output: Jac(F_T), |det(Jac(F_T)|, Jac^{-1}(F_T), |det(Jac^{-1}(F_T)|
    [x1,x2,x3] = [points[triangle[i]] for i in range(3)]
    b = x2
    A = np.dstack((x2-x1,x3-x1))
    DFT_inv = np.linalg.inv(A)
    return A, np.abs(np.linalg.det(A)), DFT_inv, np.abs(np.linalg.det(DFT_inv))



#########################
#Computing Lagrange basis
#########################

def Vandermonde_matrix(DOF,poly_order):
    #Computes Vandermonde matrix given degrees of freedom
    #for single triangle
    
    #Check if we have enough degrees of freedom
    n_DOF = len(DOF)
    n_coeffs =  np.sum(range(1,poly_order+2))
    if n_DOF == n_coeffs:
        A = np.zeros((n_DOF,n_DOF))
        for i in range(n_DOF):
            A[i,:] = get_row_vandermonde(DOF[i],poly_order)
    else:
        print("Not enough or too much DOF for Lagrange basis of order: ",poly_order)
    return A

def get_row_vandermonde(point,poly_order):
    #Gives back a Vandermonde row, i.e. 1 x y x^2 xy y^2 x^3 x^2y xy^2 y^3 ...
    x,y = point[0],point[1]
    row_vandermonde = [1.]
    #row_vandermonde = [0 ]
    for sub_order in range(1,poly_order+1):
        #Adding terms of order=sub_order to Vandermonde row
        for k in range(sub_order+1):
            #row_vandermonde.append((sub_order-k,k))
            row_vandermonde.append(x**(sub_order-k)*y**(k))
    return np.array(row_vandermonde)

#binomial_prod_without_coeff(1,1,3)   
print(Vandermonde(DOF=[(0,0),(1,0),(0,1)],poly_order=1))
