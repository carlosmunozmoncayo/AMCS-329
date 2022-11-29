#Use the finite element method to solve the 2d Poisson problem lambda(u)=-f in omega and u=g in \partial{Omega}
#For the moment considering P1 elements 

def main():
    #External libraries
    import numpy as np
    from scipy.spatial import Delaunay #Just to get triangulation
    import matplotlib.pyplot as plt

    #My libraries
    from pull_back import pull_back_forward
    from DOF import get_all_DOF_reference, get_edges_DOF, get_inner_DOF
    from MyTri import MyTri, affine_mapping
    from problems import get_Dirichlet_problem
    from assemble_Ab import Assemble_Ab

    #Defining number of nodes: Np=mx*my
    mx = 10
    my = 10
    #Creating square domain
    xmin = -1.5
    xmax = 1.5
    ymin = -1.5
    ymax = 1.5
    x = np.linspace(xmin,xmax,mx)
    y = np.linspace(ymin,ymax,my)
    #Create a grid (similar to cartesian product)
    X,Y = np.meshgrid(x,y)
    X,Y = X.flatten(), Y.flatten()
    vertices= np.dstack((X,Y)).reshape(-1, 2)
    polynomial_degree = 2
    tri = Delaunay(vertices)
    #Retrieving DOFs and available indexes
    DOF_edges = get_edges_DOF(tri, vertices, polynomial_degree-1)
    all_DOF_reference, DOF_ref_inner, DOF_ref_edges, DOF_ref_vertices = get_all_DOF_reference(polynomial_degree)
    DOF_inner, interior_index = get_inner_DOF(tri,vertices,DOF_ref_inner)

    #Defining object of class MyTri to store our data
    mytri = MyTri(vertices_index=tri.simplices, vertices_points=vertices,
                 interior_index=interior_index, interior_points=DOF_inner,
                 edges_points=DOF_edges)
    mytri.edges_points
    mytri.fill_affine_list()
    mytri.fill_edges_DOF(DOF_edges_reference=DOF_ref_edges)
    mytri.fill_DOF_ref_triangle(all_DOF_reference_triangle=all_DOF_reference)
    mytri.fill_shape_function_list(poly_degree=polynomial_degree)
    mytri.get_boundary_indexes(scipy_tri=tri,poly_degree=polynomial_degree)
    mytri.fill_shape_function_list(poly_degree=polynomial_degree)
    print(f"NDOF: {len(mytri.global_DOF)}")

    u_exact, f, g = get_Dirichlet_problem(num=1)

    A,b = Assemble_Ab(mytri=mytri, f=f, poly_degree=polynomial_degree)

    for i in mytri.boundary_indexes:
        A[i,:] = 0.
        A[i,i] = 1.
        x,y = mytri.global_DOF[i]
        b[i] = g(x,y)

    u_FE = np.linalg.solve(A,b)

    from plot import plot
    prob_type=1
    plot(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,prob_type=prob_type, polynomial_degree=polynomial_degree,
         X=X,Y=Y,mytri=mytri,tri=tri,u_exact=u_exact,u_FE=u_FE,g=g)

if __name__=="__main__":
   main() 
