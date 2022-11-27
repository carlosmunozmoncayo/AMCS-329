import numpy as np
from DOF import get_all_DOF_reference
#########################
# Computing Lagrange basis
#########################

#########################
#Main functions
#########################
def set_phi(coeffs,exps):
    #Input: coefficients (e.g., a row of get_coefficients_Lagrange_basis)
    #       and exponents (e.g., the first output of get_exponents) for a polynomial function
    #Output: A polynomial function that can be evaluated at a point (x,y) in R2
    def phi(p):
        x = p[0]**exps[:,0]
        y = p[1]**exps[:,1]
        s = np.sum(coeffs*x*y) 
        return(s)
    return phi

def set_grad_phi(coeffs,exps,coeffs_dx,exps_dx, coeffs_dy,exps_dy):
    def grad_phi(p):
        xdx = p[0]**exps_dx
        ydx = p[1]**exps[:,1]
        
        xdy = p[0]**exps[:,0]
        ydy = p[1]**exps_dy
        return np.array([np.sum(coeffs*coeffs_dx*xdx*ydx),np.sum(coeffs*coeffs_dy*xdy*ydy)])
    return grad_phi



#######################
#Auxiliary functions 1
#######################
def get_coefficients_Lagrange_basis(poly_degree=1, DOF=[(0,0),(1,0),(0,1)], reference_tri=True):
    if reference_tri: #and poly_degree>1:
        DOF = get_all_DOF_reference(poly_degree)[0]
    
    A = Vandermonde_matrix(DOF, poly_degree)
    #The ith row of alpha will contain the coefficients of the 
    #ith Lagrange function phi_i, where
    #phi_i(x,y) = a0 + a1 x + a2 y + a3 x^2  + a4 xy + a5 y^2 + a6 x^3 + ...
    alpha = []
    n_DOF = len(DOF)
    for i in range(n_DOF):
        #We'll have function with value 1 on each DOF, so we
        #we will need to solve n = len(DOF) systems with our Vandermonde
        #matrix on the LHS and a Dirac delta on the RHS
        b = np.zeros(n_DOF)
        b[i] = 1.
        alpha.append(np.linalg.solve(A,b))
    return np.array(alpha), DOF

def get_exponents(poly_degree):
    #Gives back a Vandermonde row, i.e. 1 x y x^2 xy y^2 x^3 x^2y xy^2 y^3 ...
    exponents = [(0,0)]
    exponentsdx = [(0,0)]
    exponentsdy = [(0,0)]
    for sub_degree in range(1,poly_degree+1):
        #Adding terms of degree=sub_degree to Vandermonde row
        for k in range(sub_degree+1):
            exponents.append((sub_degree-k,k))
    exponents = np.array(exponents)
    return np.array(exponents) 

def get_coefficients_grad_phi(poly_degree):
    #Automatic differentiation of a polynomial
    #Given a polynomial function phi given by
    #phi_i(x,y) = a0 + a1 x + a2 y + a3 x^2  + a4 xy + a5 y^2 + a6 x^3 + ...
    #This function returns the coefficients that will go in front of each term
    #when differentiating with respect to x and y, i.e., it gives
    #c0,c1,c2,..., and d0,d1,d2,... such that
    #\partial_x(phi(x,y))= a0 + c1a1 x + c2a2 y + c3a3 x^2  + ...
    #\partial_y(phi(x,y))= a0 + d1a1 x + d2a2 y + d3a3 x^2  + ...
    #Also returns the corresponding new exponents
    n_coeffs =  np.sum(range(1,poly_degree+2))

    #Obtaining exponents of polynomial
    #exponents = [(0,0)]
    coeffs_partial_x = np.zeros(n_coeffs,dtype=np.intc)
    coeffs_partial_y = np.zeros(n_coeffs,dtype=np.intc)
    
    i=1
    for sub_degree in range(1,poly_degree+1):
        for k in range(sub_degree+1):
            coeffs_partial_x[i] = sub_degree-k
            coeffs_partial_y[i] = k
            i+=1

    #The exponents of x in phi_x
    exponents_x = np.where(coeffs_partial_x>0,coeffs_partial_x-1,0)
    #The exponents of y in phi_y
    exponents_y = np.where(coeffs_partial_y>0,coeffs_partial_y-1,0)

    return coeffs_partial_x, exponents_x, coeffs_partial_y, exponents_y



#######################
#Auxiliary functions 2
#######################
def Vandermonde_matrix(DOF, poly_degree):
    #Computes Vandermonde matrix given degrees of freedom
    #for single triangle
    
    #Check if we have enough degrees of freedom
    n_DOF = len(DOF)
    n_coeffs =  np.sum(range(1,poly_degree+2))
    if n_DOF == n_coeffs:
        A = np.zeros((n_DOF,n_DOF))
        for i in range(n_DOF):
            A[i,:] = get_row_vandermonde_matrix(DOF[i],poly_degree)
    else:
        raise ValueError(f"Not enough or too much DOF ({n_DOF}) for Lagrange basis of degree: {poly_degree}. {n_coeffs} DOF needed")
    return A

def get_row_vandermonde_matrix(point, poly_degree):
    #Gives back a Vandermonde row, i.e. 1 x y x^2 xy y^2 x^3 x^2y xy^2 y^3 ...
    x,y = point[0],point[1]
    row_vandermonde = [1.]
    for sub_degree in range(1,poly_degree+1):
        #Adding terms of degree=sub_degree to Vandermonde row
        for k in range(sub_degree+1):
            row_vandermonde.append(x**(sub_degree-k)*y**(k))
    return np.array(row_vandermonde)


#A nice reference I am using
#http://femwiki.wikidot.com/elements:lagrange-elements
if __name__=="__main__":
    #print(Vandermonde_matrix(DOF=[(0,0),(1,0),(0,1),(.5,0),(.5,.5),(0,.5)],poly_degree=1))
    #print(get_coefficients_Lagrange_basis(poly_degree=2))
    print(get_coefficients_grad_phi(2))
