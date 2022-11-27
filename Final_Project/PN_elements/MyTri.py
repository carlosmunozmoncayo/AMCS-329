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
        #To store the DOF, this is just a list of np.arrays
        self.global_DOF = np.array(list(vertices_points)+list(interior_points)+list(edges_points))
        #self.global_DOF = [np.array(x) for x in self.global_DOF] Maybe not necessary

        self.vertices_points = vertices_points #Number of DOF in vertices
        self.interior_points = interior_points #Number of interior DOF
        self.edges_points    = edges_points    #Number of DOF on edges

        #List indexes of points in global DOF
        self.vertices = np.copy(vertices_index) 
        self.interior = interior_index+len(self.vertices_points)
        self.edges = [] #Don't forget to convert to np array
        
        #List of affine mappings
        self.affine_list = [] #self.fill_affine_list()

        #List of DOF in reference triangle
        self.DOF_ref_triangle = []

        #Lists to construct nodal basis at DOF in reference triangle
        self.coeffs_phi_ref_triangle = [] #[[a0,a1,a2,...],[a0,a1,a2,...],...]
        self.exponents_phi  = [] #[[(ex0, ey0), (ex1,ey1),...],[(ex0, ey0), (ex1,ey1),...],...]
        #each element of the nodal basis will have the form:
        #phi(x,y) = sum(i) ai*(x^exi)*(y^eyi)

        #List of coefficients gradients for nodal basis at DOF in reference triangle
        self.coeffs_grad_phi_ref_triangle = [] #[[a0,a1,a2,...],[a0,a1,a2,...],...]
        self.exponents_grad_phi  = [] #[[(ex0, ey0), (ex1,ey1),...],[(ex0, ey0), (ex1,ey1),...],...]

    def fill_affine_list(self):
        if (len(self.affine_list)) == 0:
            for i in range(len(self.vertices)):
                #Getting indexes
                #I cannot do x1,x2,x3 = self.global_DOF[self.vertices[i]]
                indexes = self.vertices[i]
                #Getting points
                x1,x2,x3 = self.global_DOF[indexes]#self.global_DOF[a], self.global_DOF[b], self.global_DOF[c]
                self.affine_list.append(affine_mapping(x1,x2,x3))
            
    def fill_edges_DOF(self, DOF_edges_reference):
        #Requires affine_list to be already full
        for i in range(len(self.vertices)):
            edge_indexes = []
            affine = self.affine_list[i]
            for xt in DOF_edges_reference:
                x = affine.A@xt +affine.b
                j = self.get_index_in_edge_points(x) #get index of x in edges_points
                edge_indexes.append(j)
            self.edges.append(edge_indexes)
        self.edges = np.array(self.edges,dtype=int)+len(self.vertices_points)+len(self.interior_points)

    def fill_DOF_ref_triangle(self, all_DOF_reference_triangle):
        self.DOF_ref_triangle = np.copy(all_DOF_reference_triangle)

    def get_index_in_edge_points(self, x):
        distances = np.array([np.linalg.norm(point-x) for point in self.edges_points])
        return np.argmin(distances)



        
