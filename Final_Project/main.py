import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import dok_matrix


def main():
    #Setting parameters
    print("Defining domain...")
    #Defining number of nodes: Np=mx*my
    mx = 12
    my = 13
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
    plot(X,Y,tri)
    A=assemble_stiffnes_matrix(tri).toarray()
    b=assemble_load_vector(tri,1)
    return 0

def assemble_stiffness_matrix(tri):
    points = tri.points
    Np = len(points)
    #Define sparse matrix
    A = dok_matrix((Np, Np), dtype=np.float32)

    def grad_phi(triangle):
        #returns grad phi (over K) for phi corresponding to each
        #point in the triangle K.
        #phi_i=ai+bi*x+ci*y
        #grad_phi_i=[bi,ci]^t
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
        
    #Iterating through all the triangles
    for triangle in tri.simplices:
        #Getting indexes of points in triangle
        #[i,j,k] = triangle
        coeffs,area = grad_phi(triangle)
        for a in range(3):
            for b in range(3):
                i,j = triangle[a],triangle[b]
                A[i,j] += int_prod_grad_phi(coeffs[i],coeffs[j],area) 
    return A

def assemble_load_vector(tri,f):
    return f

def plot(X,Y,tri,Z='test'):
    if Z=='test':
        #Creating a dummy function to plot over triangles
        f=lambda x,y: np.sin(x)**2*np.cos(y)**3
        #Evaluating dummy function
        Z=f(X,Y)
    print("Plotting...")
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    fig =plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.5]) 
    ######################
    #Plotting triangulation
    ######################
    ax2d = fig.add_subplot(gs[0])
    ax2d.triplot(X, Y, tri.simplices)
    ax2d.plot(X, Y, 'o')
    
    ###################################################
    #Plotting a surface over domain using triangulation
    ###################################################
    import matplotlib.tri as mtri
    #Constructing an object to plot triangulation
    plt_triang=mtri.Triangulation(x=X,y=Y,triangles=tri.simplices)
    ax3d = fig.add_subplot(gs[1],projection='3d')
    #ax=plt.axes(projection='3d')
    ax3d.plot_trisurf(X,Y,Z,triangles=plt_triang.triangles, cmap=plt.cm.winter)
    plt.tight_layout()
    plt.show()







if __name__=="__main__":
   main() 
