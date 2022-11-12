import numpy as np
from scipy.spatial import Delaunay


#Points of reference triangle:
#xr1=(0,1)
#xr2=(0,0)
#xr3=(1,0)

#Assuming an affine push forward function, i.e.,
#x = Ft(xri) = A xri+b

def pull_back(triangle,points):
    #Returns important quantities given a triangle and list of points
    [x1,x2,x3] = [points[triangle[i]] for i in range(3)]
    b = x2
    A = np.dstack((x3-x2,x1-x2))
    JFT_inv = np.linalg.inv(A)
    return JFT_inv, np.abs(np.linalg.det(JFT_inv))

    


