import numpy as np 

class EulerVec:
    def __init__(x, y, z):
        self.euler_vec = np.array([x, y, z])
    
    def to_so3():
        return np.array([[0, -self.euler_vec[2], self.euler_vec[1]], 
                        [self.euler_vec[2], 0, -self.euler_vec[0]], 
                        [-self.euler_vec[1], self.euler_vec[0], 0]])
