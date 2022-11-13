import numpy as np
#########################
# Computing Lagrange basis
#########################

def get_coefficients_Lagrange_basis(poly_order=1, DOF=[(0,0),(1,0),(0,1)],reference_tri=True):
    if reference_tri and poly_order>1:
        #Getting degrees of freedom for 
        DOF = []
        points_edge = poly_order+1 #Including vertices
        #Adding horizontal and vertical points
        delta = 1./poly_order
        x,y = -delta, -delta
        for iy in range(poly_order+1):
            y += delta
            for jx in range(poly_order+1-iy):
                x += delta
                DOF.append((x,y))
            x = -delta
    
    A = Vandermonde_matrix(DOF,poly_order)
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

def Vandermonde_matrix(DOF,poly_order):
    #Computes Vandermonde matrix given degrees of freedom
    #for single triangle
    
    #Check if we have enough degrees of freedom
    n_DOF = len(DOF)
    n_coeffs =  np.sum(range(1,poly_order+2))
    if n_DOF == n_coeffs:
        A = np.zeros((n_DOF,n_DOF))
        for i in range(n_DOF):
            A[i,:] = get_row_vandermonde_matrix(DOF[i],poly_order)
    else:
        raise ValueError(f"Not enough or too much DOF ({n_DOF}) for Lagrange basis of order: {poly_order}. {n_coeffs} DOF needed")
    return A

def get_row_vandermonde_matrix(point,poly_order):
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




if __name__=="__main__":
    #print(Vandermonde_matrix(DOF=[(0,0),(1,0),(0,1),(.5,0),(.5,.5),(0,.5)],poly_order=1))
    print(get_coefficients_Lagrange_basis(poly_order=2))
