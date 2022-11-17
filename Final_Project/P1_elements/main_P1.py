#Use the finite element method to solve the 2d Poisson problem lambda(u)=-f in omega and u=g in \partial{Omega}
#For the moment considering P1 elements 


##################
# External libraries
##################
import numpy as np
from scipy.spatial import Delaunay #To get triangulation

##################
# My libraries
##################
from problems import get_Dirichlet_problem
from assemble import assemble_stiffness_matrix_load_vector

def main(problem='Dirichlet', prob_type=1, plot_exact='no_tri', plot_line=True):
    from plot import plot
    #Setting parameters
    print("Defining domain...")
    
    #Defining problem    
    if problem == 'Dirichlet':
        u_exact,f,g = get_Dirichlet_problem(prob_type)

    #Defining number of nodes: Np=mx*my
    mx = 12
    my = 13
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
    points= np.dstack((X,Y)).reshape(-1, 2)
    print("Defining triangulation...")
    #Creating a Delaunay triangulation with the points
    tri = Delaunay(points)
    print("Constructing stiffness matrix and load vector")
    A_stiff_mat, b_load_vec = assemble_stiffness_matrix_load_vector(tri,f)
    A_stiff_mat = A_stiff_mat.toarray()

    # Find vertices at the boundary
    boundary = set()
    for i in range(len(tri.neighbors)):
        for k in range(3):
            if (tri.neighbors[i][k] == -1):
                nk1,nk2 = (k+1)%3, (k+2)%3 
                boundary.add(tri.simplices[i][nk1])
                boundary.add(tri.simplices[i][nk2])
    
    #Imposing Dirichlet BCs
    for i in boundary:
        A_stiff_mat[i,:] = 0.
        A_stiff_mat[i,i] = 1.
        b_load_vec[i] = g(points[i,0],points[i,1])  
    
    u_FE = np.linalg.solve(A_stiff_mat,b_load_vec)
    print("Error at nodes:", 0.5*((xmax-xmin)/(mx+1))*((ymax-ymin)/(my+1))*np.linalg.norm(u_FE-u_exact(X,Y),2))
    print("Plotting results")
    #plot(X,Y,tri,Z=u_exact(X,Y))
    plot(X,Y,tri,u_exact,u_FE,g, plot_exact,plot_line)

    return 0




if __name__=="__main__":
   main() 
