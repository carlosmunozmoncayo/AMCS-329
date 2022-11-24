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
    DOF_inner = []
    for triangle in triangulation.simplices:
        A,b,_,_,_ = pull_back_forward(triangle,points)
        for x in DOF_ref_inner:
            DOF_inner.append(A@x+b)
    return np.array(DOF_inner) #Need to return also indexes as required by the class I created

def get_all_DOF_reference(poly_degree):
    #Get all the degrees of freedom in the reference triangle
    #Also returns just the inner points for practical purposes
    ##########################################################

    #Adding horizontal and vertical points
    DOF_reference = []
    DOF_reference_inner = []
    delta = 1./poly_degree
    x,y = -delta, -delta
    for iy in range(poly_degree+1):
        y += delta
        for jx in range(poly_degree+1-iy):
            x += delta
            DOF_reference.append((x,y))
            if x>0. and y >0. and y+x<.9999 :
                DOF_reference_inner.append((x,y))
        x = -delta
    return np.array(DOF_reference), np.array(DOF_reference_inner) #Need to return also edge points as required by the class I created

        

    

