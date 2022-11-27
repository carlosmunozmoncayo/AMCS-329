import numpy as np
from scipy.spatial import Delaunay #To get triangulation
from quadrature import integrate_f_phi, integrate_gradphiI_dot_gradphiJ
from MyTri import MyTri


def Assemble_Ab(mytri, f, poly_degree):
    n_DOF = len(mytri.global_DOF)
    A = np.zeros((n_DOF,n_DOF))
    b = np.zeros(n_DOF)
    global_DOF = mytri.global_DOF
    for triangle in range(len(mytri.vertices)):
        #Iterating through all the triangles
        vertices = mytri.vertices[triangle]
        interior = mytri.interior[triangle]
        edges = mytri.edges[triangle]
        affine = mytri.affine_list[triangle]
        indexes = np.concatenate((vertices,interior,edges))
        shape_funcs = mytri.list_shape_functions[triangle]
        for index in range(len(indexes)):
            i = indexes[index]
            shapef_i = shape_funcs[index]
            b[i] += integrate_f_phi(shapef_i,affine,f,poly_degree)
            for index2 in range(len(indexes)):
                j = indexes[index2]
                shapef_j =shape_funcs[index2]
                A[i,j] += integrate_gradphiI_dot_gradphiJ(shapef_i,shapef_j,affine,poly_degree)
    return A,b
            


        
        
        
