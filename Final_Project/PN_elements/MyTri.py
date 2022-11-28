import numpy as np
from scipy.spatial import Delaunay #To get triangulation
from pull_back import fast_pull_back
from lagrange_basis import set_phi, set_grad_phi,\
        get_coefficients_Lagrange_basis, get_coefficients_grad_phi, get_exponents



class affine_mapping:
    def __init__(self, x1,x2,x3):
        self.A, self.b, self.A_inv, self.detA, self.detA_inv = fast_pull_back(x1,x2,x3)

class shape_function:
    def __init__(self,coeffs,exps,coeffs_dx,exps_dx, coeffs_dy,exps_dy):
        self.phi = set_phi(coeffs,exps)
        self.grad_phi = set_grad_phi(coeffs,exps,coeffs_dx,exps_dx, coeffs_dy,exps_dy)

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


        #List of lists of shape_functions objects [[sf,sf,sf...],[sf,sf,sf,...],...]
        #Each list contains the shape_function object corresponding to its degrees of freedom
        #Using the indexing [vertices, interior, edges]
        self.list_shape_functions = []


        #List of boundary indexes wrt to global_DOF
        self.boundary_indexes = []


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

    def fill_shape_function_list(self, poly_degree):
        if len(self.list_shape_functions)>0:
            return 0
        for i in range(len(self.vertices)):
            #Getting all the DOF for this triangle
            DOF_triangle = np.concatenate((self.global_DOF[self.vertices[i]],self.global_DOF[self.interior[i]],self.global_DOF[self.edges[i]]))
            #This is done in a simple setting in one of the Jupyter notebooks (Checking_Lagrange_basis)
            coeffs,_ = get_coefficients_Lagrange_basis(poly_degree=poly_degree, DOF=DOF_triangle, reference_tri=False)
            coeffs_dx,exps_dx, coeffs_dy,exps_dy = get_coefficients_grad_phi(poly_degree=poly_degree)
            exps = get_exponents(poly_degree=poly_degree)
            
            #Local list of shape functions
            local_list_sp = []
            
            for coeff_row in coeffs:
                local_list_sp.append(shape_function(coeff_row,exps,coeffs_dx,exps_dx, coeffs_dy,exps_dy))
            self.list_shape_functions.append(local_list_sp)

    def get_boundary_indexes(self, scipy_tri, poly_degree):
        if len(self.boundary_indexes) >0:
            return 0
        # Find edges at the boundary
        boundary_indexes = set()
        points_edges_on_boundary = [] #points that are not vertices
        for i in range(len(scipy_tri.neighbors)):
            for k in range(3):
                if (scipy_tri.neighbors[i][k] == -1):
                    nk1,nk2 = (k+1)%3, (k+2)%3
                    index1=scipy_tri.simplices[i][nk1]
                    index2=scipy_tri.simplices[i][nk2]
                    boundary_indexes.add(index1)
                    boundary_indexes.add(index2)
                    x = self.vertices_points[index1]
                    y = self.vertices_points[index2]
                    coeffs = np.linspace(0., 1., poly_degree+1)[1:-1]
                    for lamb in coeffs:
                        points_edges_on_boundary.append(lamb*x+(1-lamb)*y)
        #We already have the indexes of the boundary vertices, now we have to get
        #the indexes (wrt to global_DOF) of the edge boundary points
        n = len(self.vertices_points)+len(self.interior_points)
        for point in points_edges_on_boundary:
            index_edge_points = self.get_index_in_edge_points(point)
            boundary_indexes.add(index_edge_points+n) 

        self.boundary_indexes = np.array(list(boundary_indexes))

                



        
