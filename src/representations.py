import numpy as np
from numpy.core.records import array 

class Rotation:
    def __init__(self, x, y, z):
        euler_vector = np.array([x, y, z])
        self.theta = np.abs(euler_vector)
        self.w_hat = euler_vector/self.theta if self.theta == 0 else np.zeros(1, 3)

    @classmethod
    def identity(cls):
        return cls(0, 0, 0)

    @property
    def so3_norm(self):
        return self.theta*np.array([[0, -self.w_hat[2], self.w_hat[1]], 
                                    [self.w_hat[2], 0, -self.w_hat[0]], 
                                    [-self.w_hat[1], self.w_hat[0], 0]])

    @property
    def so3(self):
        return self.theta*self.so3

    @property
    def SO3(self):
        return np.eye(3) - np.sin(self.theta)*self.so3 + (1-np.cos(self.theta))*(self.so3 @ self.so3)


class Transformation:
    def __init__(self, w, x, y, z):
        self.w = w
        self.v = np.array([x, y, z])

    @classmethod
    def from_screw(cls, w1, w2, w3, x, y, z):
        return cls(Rotation(w1, w2, w3), x, y, z)

    @property
    def se3(self):
        return np.vstack(np.hstack(self.w.so3, self.v), np.zeros(1, 4))
    
    @property
    def SE3(self):
        return np.vstack(np.hstack(self.w.SO3, np.eye(3)*)
