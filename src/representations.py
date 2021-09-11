import numpy as np
from numpy.core.records import array 

class Rotation:
    def __init__(self, x, y, z):
        self.euler_vec = np.array([x, y, z])
    
    @property
    def so3(self):
        return np.array([[0, -self.euler_vec[2], self.euler_vec[1]], 
                        [self.euler_vec[2], 0, -self.euler_vec[0]], 
                        [-self.euler_vec[1], self.euler_vec[0], 0]])
    @property
    def SO3(self):
        theta = np.abs(self.euler_vec)
        norm_vector = self.euler_vec / theta
        return np.eye(3) - np.sin(theta) * self.so3 + (1 - np.cos(theta)) * (self.so3 @ self.so3)

class Transformation:

    def __init__(self):
        pass

    @property
    def se3(self):
        pass
    
    @property
    def SE3(self):
        pass
