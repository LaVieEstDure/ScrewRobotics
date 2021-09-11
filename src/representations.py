import numpy as np
from numpy.core.records import array 

class EulerVec:
    def __init__(self, x, y, z):
        self.euler_vec = np.array([x, y, z])
    
    def to_so3(self):
        return np.array([[0, -self.euler_vec[2], self.euler_vec[1]], 
                        [self.euler_vec[2], 0, -self.euler_vec[0]], 
                        [-self.euler_vec[1], self.euler_vec[0], 0]])

class ScrewVec:
    def __init__(self, wx, wy, wz, vx, vy, vz):
        self.screw_vec = np.array([wx, wy, wz, vx, vy, vz])
    
    def to_se3(self):
        pass

class Transformation:

    def __init__(self):
        pass

    @property
    def se3(self):
        pass
    
    @property
    def SE3(self):
        pass