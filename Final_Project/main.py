import numpy as np
from scipy.spatial import Delaunay

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
    return 0

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
