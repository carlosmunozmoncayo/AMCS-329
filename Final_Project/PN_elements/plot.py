
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import matplotlib.tri as mtri
from problems import get_Dirichlet_problem

def plot(xmin,xmax,ymin,ymax,prob_type, polynomial_degree,
         X,Y,mytri,tri,u_exact,u_FE,g):
    original_poly_degree = polynomial_degree
    fig = plt.figure(dpi=150,figsize=(8,8))
    axs = []

    ###################################################
    #Plotting P^n solution over domain using triangulation
    ###################################################
    #Constructing an object to plot triangulation
    axs.append(fig.add_subplot(221,projection='3d'))
    
    #Plot FE solution over triangulation
    axs[0].plot_trisurf(mytri.global_DOF[:,0],mytri.global_DOF[:,1],u_FE, cmap=plt.cm.winter)
    axs[0].set_zlim(-1., 1.)
    axs[0].set_title(f"FE $P^{polynomial_degree}$ solution")



    ################
    #Retrieving some data
    ################
    #Number of nodes: Np=mx*my
    mx = int(np.sqrt(len(X)))
    my = int(np.sqrt(len(Y)))

    ###################################################
    #Plotting a exact solution over domain using triangulation
    ###################################################
    axs.append(fig.add_subplot(222,projection='3d'))
    #Defining number of nodes: Np=mx*my
    delta = min(30,3*min(mx,my))
    #Creating square domain
    x = np.linspace(xmin,xmax,delta)
    y = np.linspace(ymin,ymax,delta)
    #Create a grid (similar to cartesian product)
    X1,Y1 = np.meshgrid(x,y)
    X2,Y2 = X1.flatten(), Y1.flatten()
    axs[1].plot_trisurf(X2,Y2,u_exact(X2,Y2), cmap=plt.cm.winter)
    axs[1].set_title('Exact solution')
    



    #################
    #Getting  P1 solution
    ################
    #External libraries
    from scipy.spatial import Delaunay #Just to get triangulation

    #My libraries
    from pull_back import pull_back_forward
    from DOF import get_all_DOF_reference, get_edges_DOF, get_inner_DOF
    from MyTri import MyTri, affine_mapping
    from assemble_Ab import Assemble_Ab
    x = np.linspace(xmin,xmax,mx)
    y = np.linspace(ymin,ymax,my)
    #Create a grid (similar to cartesian product)
    X,Y = np.meshgrid(x,y)
    X,Y = X.flatten(), Y.flatten()
    vertices= np.dstack((X,Y)).reshape(-1, 2)
    tri = Delaunay(vertices)
    #Retrieving DOFs and available indexes
    polynomial_degree = 1
    DOF_edges = get_edges_DOF(tri, vertices, polynomial_degree-1)
    all_DOF_reference, DOF_ref_inner, DOF_ref_edges, DOF_ref_vertices = get_all_DOF_reference(polynomial_degree)
    DOF_inner, interior_index = get_inner_DOF(tri,vertices,DOF_ref_inner)

    #Defining object of class MyTri to store our data
    mytri1 = MyTri(vertices_index=tri.simplices, vertices_points=vertices,
                 interior_index=interior_index, interior_points=DOF_inner,
                 edges_points=DOF_edges)
    mytri1.edges_points
    mytri1.fill_affine_list()
    mytri1.fill_edges_DOF(DOF_edges_reference=DOF_ref_edges)
    mytri1.fill_DOF_ref_triangle(all_DOF_reference_triangle=all_DOF_reference)
    mytri1.fill_shape_function_list(poly_degree=polynomial_degree)
    mytri1.get_boundary_indexes(scipy_tri=tri,poly_degree=polynomial_degree)
    mytri1.fill_shape_function_list(poly_degree=polynomial_degree)

    u_exact, f, g = get_Dirichlet_problem(prob_type)

    A,b = Assemble_Ab(mytri=mytri1, f=f, poly_degree=1)
    for i in mytri1.boundary_indexes:
        A[i,:] = 0.
        A[i,i] = 1.
        x,y = mytri1.global_DOF[i]
        b[i] = g(x,y)
    u_FE_P1 = np.linalg.solve(A,b) 

    #Constructing an object to plot triangulation
    plt_triang=mtri.Triangulation(x=X,y=Y,triangles=tri.simplices)

    ###################################################
    #Plotting a surface over domain using triangulation
    ###################################################
    #Constructing an object to plot triangulation
    axs.append(fig.add_subplot(223,projection='3d'))
    #Plot FE solution over triangulation
    axs[2].plot_trisurf(X,Y,u_FE_P1,triangles=plt_triang.triangles, cmap=plt.cm.winter)
    axs[2].set_title('FE $P^1$ solution')


    

    ######################
    #Plotting triangulation
    ######################
    axs.append(fig.add_subplot(224))
    axs[3].triplot(vertices[:,0],vertices[:,1],tri.simplices, color='green', lw=0.8)
    if original_poly_degree>2:
        axs[3].plot(mytri.interior_points[:,0],mytri.interior_points[:,1], 'b.', markersize=2)
    if original_poly_degree>1:
        axs[3].plot(mytri.edges_points[:,0],mytri.edges_points[:,1], 'r.' , markersize = 2)
    axs[3].plot(mytri.vertices_points[:,0],mytri.vertices_points[:,1], 'k.', markersize = 2)
    #Plotting boundary points
    axs[3].plot(mytri.global_DOF[mytri.boundary_indexes][:,0],
            mytri.global_DOF[mytri.boundary_indexes][:,1],'yo', markersize = 2)
    axs[3].axis('scaled')
    plt.tight_layout() 
    plt.show()
    plt.close()










