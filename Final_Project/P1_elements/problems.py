import numpy as np


def get_Dirichlet_problem(num):
    if num == 1 :
        a=1
        f = lambda x,y : 2*(a*np.pi)**2*np.sin(x*a*np.pi)*np.sin(y*a*np.pi)
        u_exact = lambda x,y, : np.sin(x*a*np.pi)*np.sin(y*a*np.pi)
        g = lambda x,y : u_exact(x,y)
    elif num == 2 :
        u_exact = lambda x,y, : np.sin(x+y)*np.cos(x*y)
        f = lambda x,y : (x**2+y**2+2)*np.sin(x+y)*np.cos(x*y)+2*(x+y)*np.sin(x+y)*np.cos(x*y)
        g = lambda x,y : u_exact(x,y)
    return u_exact,f,g 


