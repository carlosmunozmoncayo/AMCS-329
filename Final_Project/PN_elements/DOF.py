#Defining degrees of freedom that are not at the vertices

import itertools

def get_edges_DOF(triangulation, points, k):

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
        sets_of_pairs = [set(i) for i in itertools.combinations(triangle,2)] #Combinations of 2 points (list of sets of 2 indices)
        for pair in sets_of_pairs:
            if pair not in indexes_connected_points:
                indexes_connected_points.append(pair)
                i,j = tuple(pair)
                x,y = points[i], points[j]
                for lamb in coeffs:
                    DOF_on_edges.append(lamb*x+(1-lamb)*y)
    return DOF_on_edges
