import numpy as np
from scipy.sparse import dok_matrix

def assemble_stiffness_matrix_load_vector(tri,f):
    points = tri.points
    Np = len(points)
    #Define sparse stiffness matrix
    A_stiff_mat = dok_matrix((Np, Np), dtype=np.float32)
    #Define load vector
    b_load_vec = np.zeros(Np)

    def grad_phi(triangle):
        #returns grad phi (over K) for phi corresponding to each
        #point in the triangle K.
        #phi_i=ai+bi*x+ci*y (over K)
        #grad_phi_i=[bi,ci]^t (over K)
        [i,j,k] = triangle
        [xi,yi],[xj,yj],[xk,yk] = points[i],points[j],points[k]
        twice_area=np.abs(xi*(yj-yk)+xj*(yk-yi)+xk*(yi-yj))
        bi,ci=(yj-yk)/twice_area,(xk-xj)/twice_area
        bj,cj=(yk-yi)/twice_area,(xi-xk)/twice_area
        bk,ck=(yi-yj)/twice_area,(xj-xi)/twice_area
        return [[bi,ci],[bj,cj],[bk,ck]],twice_area/2.

    def int_prod_grad_phi(coeffs_i,coeffs_j,area):
        bi,bj,ci,cj = coeffs_i[0],coeffs_j[0],coeffs_i[1],coeffs_j[1]
        return (bi*bj+ci*cj)*area
    
    def int_prod_f_phi(triangle,area,f):
        #This returns the integral (over K) of f*phi for each point
        #in the triangle K.
        #Just using a node quadrature for the moment
        [i,j,k] = triangle
        [xi,yi],[xj,yj],[xk,yk] = points[i],points[j],points[k]
        return [(1./3.)*f(xi,yi)*area,(1./3.)*f(xj,yj)*area,(1./3.)*f(xk,yk)*area]
        
    #Iterating through all the triangles
    for triangle in tri.simplices:
        coeffs,area = grad_phi(triangle) #For stifness matrix
        loads = int_prod_f_phi(triangle,area,f)  #For load vector
        for a in range(3):
            i = triangle[a]
            b_load_vec[i] += loads[a]
            for b in range(3):
                j = triangle[b]
                A_stiff_mat[i,j] += int_prod_grad_phi(coeffs[a],coeffs[b],area)
                
    return A_stiff_mat, b_load_vec
