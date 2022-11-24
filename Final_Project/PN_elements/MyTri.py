import numpy as np
from scipy.spatial import Delaunay #To get triangulation
from pull_back import fast_pull_back



class affine_mapping:
    def __init__(self, x1,x2,x3):
        self.A, self.b, self.A_inv, self.detA, self.detA_inv = fast_pull_back(x1,x2,x3)

class MyTri:
    def __init__(self, 
                 vertices_index, vertices_points,
                 interior_index, interior_points, 
                 edges_points):
        #To store the DOF
        self.global_DOF = list(vertices_points)+list(interior_points)+list(edges_points)

        self.nvertices = len(vertices_points) #Number of DOF in vertices
        self.ninterior = len(interior_points) #Number of interior DOF
        self.nedges = len(edges_points)    #Number of DOF on edges

        self.vertices = np.copy(vertices)
        self.interior = interior_index+self.nvertices
        self.edges = []

        self.affine_list = []

    def fill_affine_list():

    def fill_edges_DOF(DOF_edges_reference):
        #Requires affine_list to be already full
        for i in range(nvertices):
            edge_indexes = []
            affine = self.affine_list[i]
            for xt in DOF_edges_reference:
                x = affine.A@xt +affine.b
                j = THIS DOES NOT EXIST YET #get index of x in edges_points
                edge_indexes.append(j)
            self.edges.append(edge_indexes)
        self.edges = np.array(self.edges,dtype=int)+selft.nvertices+self.ninterior


        
