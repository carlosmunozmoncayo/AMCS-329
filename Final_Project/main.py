import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import dok_matrix
#Use the finite element method to solve the Poisson problem lambda(u)=-f in 2D
#For the moment considering P1 elements and homogeneous boundary conditions

def main():
    #Setting parameters
    print("Defining domain...")
    #RHS function
    #f = lambda x,y : -2*(x**2+y**2-2) #RHS of problem with sol. (-x**2+1)*(-y**2+1)
    #u_exact = lambda x,y : (-x**2+1)*(-y**2+1)
    u_exact = lambda x,y, : np.sin(x*np.pi)*np.sin(y*np.pi)
    f = lambda x,y : 2*np.pi**2*u_exact(x,y)
    #Defining number of nodes: Np=mx*my
    mx = 13
    my = 12
#Creating square domain
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
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

    # Find edges at the boundary
    boundary = set()
    for i in range(len(tri.neighbors)):
        for k in range(3):
            if (tri.neighbors[i][k] == -1):
                nk1,nk2 = (k+1)%3, (k+2)%3 
                boundary.add(tri.simplices[i][nk1])
                boundary.add(tri.simplices[i][nk2])
    #Imposing homogeneous Dirichlet BCs
    for i in boundary:
        A_stiff_mat[i,:] = 0.
        A_stiff_mat[:,i] = 0.
        A_stiff_mat[i,i] = 1.
        b_load_vec[i] = 0.    
    u_FE = np.linalg.solve(A_stiff_mat,b_load_vec)
    print("Plotting results")
    #plot(X,Y,tri,Z=u_exact(X,Y))
    plot(X,Y,tri,u_exact(X,Y),u_FE)

    return 0

def assemble_stiffness_matrix_load_vector(tri,f):
    points = tri.points
    Np = len(points)
    #Define sparse stiffness matrix
    A_stiff_mat = dok_matrix((Np, Np), dtype=np.float32)
    #Define load vector
    b_load_vec = np.zeros(Np)

    def grad_phi(triangle):
        #returns grad phi (over K) for phi corresponding to each
        #point in the triangle K.
        #phi_i=ai+bi*x+ci*y (over K)
        #grad_phi_i=[bi,ci]^t (over K)
        [i,j,k] = triangle
        [xi,yi],[xj,yj],[xk,yk] = points[i],points[j],points[k]
        twice_area=np.abs(xi*(yj-yk)+xj*(yk-yi)+xk*(yi-yj))
        bi,ci=(yj-yk)/twice_area,(xk-xj)/twice_area
        bj,cj=(yk-yi)/twice_area,(xi-xk)/twice_area
        bk,ck=(yi-yj)/twice_area,(xj-xi)/twice_area
        return [[bi,ci],[bj,cj],[bk,ck]],twice_area/2.

    def int_prod_grad_phi(coeffs_i,coeffs_j,area):
        bi,bj,ci,cj = coeffs_i[0],coeffs_j[0],coeffs_i[1],coeffs_j[1]
        return (bi*bj+ci*cj)*area
    
    def int_prod_f_phi(triangle,area,f):
        #This returns the integral (over K) of f*phi for each point
        #in the triangle K.
        #Just using a node quadrature for the moment
        [i,j,k] = triangle
        [xi,yi],[xj,yj],[xk,yk] = points[i],points[j],points[k]
        return [(1./3.)*f(xi,yi)*area,(1./3.)*f(xj,yj)*area,(1./3.)*f(xk,yk)*area]
        
    #Iterating through all the triangles
    for triangle in tri.simplices:
        coeffs,area = grad_phi(triangle) #For stifness matrix
        loads = int_prod_f_phi(triangle,area,f)  #For load vector
        for a in range(3):
            i = triangle[a]
            b_load_vec[i] += loads[a]
            for b in range(3):
                j = triangle[b]
                A_stiff_mat[i,j] += int_prod_grad_phi(coeffs[a],coeffs[b],area)
                
    return A_stiff_mat, b_load_vec

def plot(X,Y,tri,u_exact,u_FE):
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    #fig, axs =plt.subplots(nrows=1,ncols=3)#,figsize=(12,6))
    fig = plt.figure()
    axs = []
    #gs = gridspec.GridSpec(nrows=1,ncols=3)  #width_ratios=[1, 1.5,]

    import matplotlib.tri as mtri
    #Constructing an object to plot triangulation
    plt_triang=mtri.Triangulation(x=X,y=Y,triangles=tri.simplices)
        
    ###################################################
    #Plotting a surface over domain using triangulation
    ###################################################
    #Constructing an object to plot triangulation
    axs.append(fig.add_subplot(221,projection='3d'))
    #ax=plt.axes(projection='3d')
    axs[0].plot_trisurf(X,Y,u_FE,triangles=plt_triang.triangles, cmap=plt.cm.winter)
    axs[0].set_title('FE solution')

    ###################################################
    #Plotting a exact solution over domain using triangulation
    ###################################################
    axs.append(fig.add_subplot(222,projection='3d'))
    #ax=plt.axes(projection='3d')
    axs[1].plot_trisurf(X,Y,u_exact,triangles=plt_triang.triangles, cmap=plt.cm.winter)
    axs[1].set_title('Exact solution')


    ######################
    #Plotting triangulation
    ######################
    
    axs.append(fig.add_subplot(223))
    axs[2].triplot(X, Y, tri.simplices)
    axs[2].plot(X, Y, 'r.')
    axs[2].axis('scaled')
    axs[2].set_title('Triangulation and DOF')

    plt.tight_layout()
    plt.show()




if __name__=="__main__":
   main() 
