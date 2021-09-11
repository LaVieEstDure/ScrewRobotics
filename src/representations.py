import numpy as np
from numpy.core.records import array 
from math import *
from functools import lru_cache
tol = 1e-9

class Rotation:
    def __init__(self, wx, wy, wz):
        w = np.array([[wx, wy, wz]]).T
        self.theta = np.linalg.norm(w)
        self.w_hat = w / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T

    @classmethod
    def identity(cls):
        return cls(0, 0, 0)

    @property
    @lru_cache(maxsize=1)
    def so3_norm(self):
        return np.array([[             0, -self.w_hat[2],  self.w_hat[1]], 
                         [ self.w_hat[2],              0, -self.w_hat[0]], 
                         [-self.w_hat[1],  self.w_hat[0],              0]])

    '''
    Returns the transformation in so3 (Lie Algebra)
    '''
    @property
    @lru_cache(maxsize=1)
    def so3(self):
        return self.theta * self.so3_norm

    '''
    Returns the transformation in SO3 (Lie Group)
    '''
    @property
    @lru_cache(maxsize=1)
    def SO3(self):
        R = np.eye(3) + sin(self.theta) * self.so3_norm + (1 - cos(self.theta)) * (self.so3_norm @ self.so3_norm)
        R[np.abs(R) < tol] = 0.0
        return R

class Transformation:
    def __init__(self, wx, wy, wz, vx, vy, vz):
        w = np.array([[wx, wy, wz]]).T
        v = np.array([[vx, vy, vz]]).T
        mag_w = np.linalg.norm(w)
        mag_v = np.linalg.norm(v)
        '''
        Set theta to be the magnitude of w
        If w = 0 (no rotation) set theta to be the magnitude of v
        If v = 0 (no translation) set theta to be 0 ¯\_(ツ)_/¯
        self.theta = mag_w if mag_w > tol else mag_v if mag_v > tol else 0
        '''
        self.theta = mag_w if mag_w > tol else mag_v if mag_v > tol else 0
        self.w_hat = w / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T
        self.v_hat = v / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T

        self.rotation = Rotation(wx, wy, wz)

    @property
    @lru_cache(maxsize=1)
    def se3_norm(self):
        return np.vstack((np.hstack((self.rotation.so3_norm, self.v_hat)), np.array([[0, 0, 0, 0]])))
    

    '''
    Returns the transformation in so3 (Lie Algebra)
    '''
    @property
    @lru_cache(maxsize=1)
    def se3(self):
        return self.theta * self.se3_norm

    @property
    @lru_cache(maxsize=1)
    def SE3(self):
        R = np.eye(3) + sin(self.theta) * self.rotation.so3_norm + (1 - cos(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm)
        p = (np.eye(3) * self.theta + (sin(self.theta) * self.rotation.so3_norm + (1 - cos(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm))) @ self.v_hat
        T = np.vstack((np.hstack((R, p)), np.array([[0, 0, 0, 1]])))
        T[np.abs(T) < tol] = 0.0
        return T