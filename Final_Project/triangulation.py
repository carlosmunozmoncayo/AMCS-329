import numpy as np
from scipy.spatial import Delaunay

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



if __name__=="__main__":
    print("Plotting...")
    import matplotlib.pyplot as plt
    fig =plt.figure(figsize=(10,5))
    ######################
    #Plotting triangulation
    ######################
    ax2d = fig.add_subplot(121)
    ax2d.triplot(points[:,0], points[:,1], tri.simplices)
    ax2d.plot(points[:,0], points[:,1], 'o')
    
    ###################################################
    #Plotting a surface over domain using triangulation
    ###################################################
    import matplotlib.tri as mtri
    #Creating a dummy function to plot over triangles
    f=lambda x,y: np.sin(x)**2*np.cos(y)**3
    #Evaluating dummy function
    Z=f(X,Y)
    #Constructing an object to plot triangulation
    plt_triang=mtri.Triangulation(x=X,y=Y,triangles=tri.simplices)
    ax3d = fig.add_subplot(122,projection='3d')
    #ax=plt.axes(projection='3d')
    ax3d.plot_trisurf(X,Y,Z,triangles=plt_triang.triangles, cmap=plt.cm.winter)
    plt.show()
