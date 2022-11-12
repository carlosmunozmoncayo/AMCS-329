
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import matplotlib.tri as mtri

def plot(X,Y,tri,u_exact,u_FE,g,plot_exact,plot_line):
    
    #fig, axs =plt.subplots(nrows=1,ncols=3)#,figsize=(12,6))
    fig = plt.figure()
    axs = []
    #gs = gridspec.GridSpec(nrows=1,ncols=3)  #width_ratios=[1, 1.5,]

    
    #Constructing an object to plot triangulation
    plt_triang=mtri.Triangulation(x=X,y=Y,triangles=tri.simplices)

    ################
    #Retrieving some data
    ################
    #Number of nodes: Np=mx*my
    mx = int(np.sqrt(len(X)))
    my = int(np.sqrt(len(Y)))
    #Creating square domain
    xmin = X[0]
    xmax = X[-1]
    ymin = Y[0]
    ymax = Y[-1]



        
    ###################################################
    #Plotting a surface over domain using triangulation
    ###################################################
    #Constructing an object to plot triangulation
    axs.append(fig.add_subplot(221,projection='3d'))
    #Plot FE solution over triangulation
    axs[0].plot_trisurf(X,Y,u_FE,triangles=plt_triang.triangles, cmap=plt.cm.winter)

    if plot_line:
        #Plotting a line on boundary
        x = np.linspace(xmin,xmax,mx)
        y = np.linspace(ymin,ymax,my)
        X1,Y1 = np.meshgrid(x,y)
        for i in [0,-1]:
            axs[0].plot(X1[i,:],Y1[i,:],g(X1[i,:],Y1[i,:]),color='k',lw=4)
            axs[0].plot(X1[:,i],Y1[:,i],g(X1[:,i],Y1[:,i]),color='k',lw=4) 
    axs[0].set_title('FE solution')



    ###################################################
    #Plotting a exact solution over domain using triangulation
    ###################################################
    axs.append(fig.add_subplot(222,projection='3d'))
    #ax=plt.axes(projection='3d')
    if plot_exact == 'tri':
        axs[1].plot_trisurf(X,Y,u_exact(X,Y),triangles=plt_triang.triangles, cmap=plt.cm.winter)
    else:
        #Defining number of nodes: Np=mx*my
        delta = min(30,3*min(mx,my))
        #Creating square domain
        x = np.linspace(xmin,xmax,delta)
        y = np.linspace(ymin,ymax,delta)
        #Create a grid (similar to cartesian product)
        X1,Y1 = np.meshgrid(x,y)
        X2,Y2 = X1.flatten(), Y1.flatten()
        axs[1].plot_trisurf(X2,Y2,u_exact(X2,Y2), cmap=plt.cm.winter)
        if plot_line:
            #Plotting line in boundary
            for i in [0,-1]:
                axs[1].plot(X1[i,:],Y1[i,:],u_exact(X1[i,:],Y1[i,:]),color='k',lw=4)
                axs[1].plot(X1[:,i],Y1[:,i],u_exact(X1[:,i],Y1[:,i]),color='k',lw=4)
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


