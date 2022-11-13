import numpy as np

############################
# Pull back and push forward functions 
############################

#Points of reference triangle:
#xr1=(0,0)
#xr2=(1,0)
#xr3=(0,1)

#Assuming an affine push forward function, i.e.,
#x = Ft(xr) = A xr+b

def p_back_forward(triangle,points):
    #Returns important quantities given a triangle and list of points
    #output: Jac(F_T), |det(Jac(F_T)|, Jac^{-1}(F_T), |det(Jac^{-1}(F_T)|

    #Extremely inefficient btw, it would be better to receive the three points
    #of the triangle directly isntead of passing the list of points each time.
    #It suffices for our simple applications anyway.
    [x1,x2,x3] = [points[triangle[i]] for i in range(3)]
    b = x2
    A = np.dstack((x2-x1,x3-x1))
    DFT_inv = np.linalg.inv(A)
    return A, np.abs(np.linalg.det(A)), DFT_inv, np.abs(np.linalg.det(DFT_inv))



