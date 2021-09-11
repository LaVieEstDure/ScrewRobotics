import numpy as np
from numpy.core.records import array 

class w:
    def __init__(self, x, y, z):
        self.euler_vec = np.array([x, y, z])

    @classmethod
    def identity():
        self.euler_vec = np.array([0, 0, 0])

    @property
    def so3(self):
        return np.array([[0, -self.euler_vec[2], self.euler_vec[1]], 
                        [self.euler_vec[2], 0, -self.euler_vec[0]], 
                        [-self.euler_vec[1], self.euler_vec[0], 0]])
    @property
    def SO3(self):
        theta = np.abs(self.euler_vec)
        return np.eye(3) - np.sin(theta)*self.so3/theta + (1-np.cos(theta))*(self.so3/theta @ self.so3/theta)

class Transformation:
    def __init__(self, w1, w2, w3, x, y, z):
        self.w = w(w1, w2, w3)
        self.v = np.array([x, y, z])

    @classmethod
    def from_w_xyz(w, x, y, z):
        self.w = w 
        self.v = np.array([x, y, z]).T

    @property
    def se3(self):
        return np.vstack(np.hstack(self.w.so3, self.v), np.zeros(1, 4))
    
    @property
    def SE3(self):
        pass
