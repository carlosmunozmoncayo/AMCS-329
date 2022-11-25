#Defining degrees of freedom that are not at the vertices
import numpy as np
from itertools import combinations
from pull_back import pull_back_forward


def get_edges_DOF(triangulation, points, k):
    #Get degrees of freedom over all the domain. This function returns the DOF
    #that are not vertices or inner points for each triangle.
    #############
    #k is the number of interior poins at each edge, which is
    #1 less than the degree of the polynomials from our basis
    #############
    
    indexes_connected_points = [] #To not repeat the same edge from two triangles
    #list of sets where each set has the form {i,j}, this indicates that the points points[i] and points[j]
    #are connected through an edge

    DOF_on_edges = [] #To store DOF on edges, list of 1x2 arrays
    
    #########
    #Getting coefficiens for convex combinations (keeping just interior points)
    ########
    coeffs = np.linspace(0., 1., k+2)[1:-1] 

    for triangle in triangulation.simplices:
        sets_of_pairs = [set(i) for i in combinations(triangle,2)] #Combinations of 2 points (list of sets of 2 indices)
        for pair in sets_of_pairs:
            if pair not in indexes_connected_points:
                indexes_connected_points.append(pair)
                i,j = tuple(pair)
                x,y = points[i], points[j]
                for lamb in coeffs:
                    DOF_on_edges.append(lamb*x+(1-lamb)*y)
    return np.array(DOF_on_edges)

def get_inner_DOF(triangulation, points, DOF_ref_inner):
    ######################
    #Returns a list of all DOF that are interior points
    #and a list of lists of indexes that indicate to which triangle each
    #interior DOF belongs
    ######################
    DOF_inner = []
    interior_index = []
    i = 0
    for triangle in triangulation.simplices:
        triangle_indexes = []
        A,b,_,_,_ = pull_back_forward(triangle,points)
        for x in DOF_ref_inner:
            DOF_inner.append(A@x+b)
            triangle_indexes.append(i)
            i+=1
        interior_index.append(triangle_indexes)
    return np.array(DOF_inner), np.array(interior_index, dtype=int) 

def get_all_DOF_reference(poly_degree):
    #Get all the degrees of freedom in the reference triangle
    #Also returns just the inner points for practical purposes
    ##########################################################
    tol = 1.e-5
    #Adding horizontal and vertical points
    DOF_reference = []
    DOF_reference_inner = []
    DOF_reference_edges = []
    DOF_reference_vertices = []
    delta = 1./poly_degree
    x,y = -delta, -delta
    for iy in range(poly_degree+1):
        y += delta
        for jx in range(poly_degree+1-iy):
            x += delta
            DOF_reference.append((x,y))
            if x>tol and y >tol and y+x<(1-tol) :
                DOF_reference_inner.append((x,y))
            elif np.allclose([0.,0.],[x,y]) or np.allclose([0.,1.],[x,y]) or np.allclose([1.,0.],[x,y]) :
                DOF_reference_vertices.append((x,y))
            else:
                DOF_reference_edges.append((x,y))
        x = -delta
    return (np.array(DOF_reference),
            np.array(DOF_reference_inner),
            np.array(DOF_reference_edges),
            np.array(DOF_reference_vertices))
        

    

