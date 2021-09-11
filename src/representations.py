import numpy as np
from numpy.core.records import array 
from math import *
from functools import lru_cache
<<<<<<< HEAD
from copy import deepcopy
tol = 1e-9
=======
# Tolerance
tol = 1e-10
>>>>>>> 814ba73a4aaf8a0cc012e6d779457e402a8fef30

class Rotation:
    def __init__(self, wx, wy, wz):
        """
        Defines a rotation from its Euler Vector
        w = [wx wy wz]^T

        The rotation is given by an axis of rotation w_hat and angle of rotation theta
        
        The euler vector is decomposed as w = w_hat * theta,
        where theta is the angle of rotation
        and w_hat is the axis of rotation

        If w = 0, the axis is undefined and the angle is zero
        This corresponds to the identity rotation matrix.
        """
        w = np.array([[wx, wy, wz]]).T
        self.theta = np.linalg.norm(w)
        self.w_hat = w / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T

    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0, 0, 0)

    @property
    @lru_cache(maxsize=1)
    def so3_norm(self):
        """Returns the transformation in so3 (element of Lie Algebra), when theta = 1"""
        return np.array([[             0, -self.w_hat[2],  self.w_hat[1]], 
                         [ self.w_hat[2],              0, -self.w_hat[0]], 
                         [-self.w_hat[1],  self.w_hat[0],              0]])

    @property
    @lru_cache(maxsize=1)
    def so3(self):
        """Returns the transformation in so3 (element of Lie Algebra)"""
        return self.theta * self.so3_norm

    @property
    @lru_cache(maxsize=1)
    def SO3(self):
        """Returns the transformation in SO3 (element of Lie Group)"""
        R = np.eye(3) + sin(self.theta) * self.so3_norm + (1 - cos(self.theta)) * (self.so3_norm @ self.so3_norm)
        R[np.abs(R) < tol] = 0.0
        return R

class Transformation:
    def __init__(self, wx, wy, wz, vx, vy, vz):
        """
        Defines a transformation (rotation + translation) from its Screw Vector
        w = [wx wy wz vx vy vz]^T

        The rotation is given by an axis of rotation w_hat and angle of rotation theta
        The translation has two terms dueto the linear motion at the origin induced by
        rotation about the axis, and the other due to the translation along the screw axis
        
        The screw vector is decomposed as V = S_hat * theta,
        where theta is the angle of rotation for revolute joints
        or the distance of translation for prismatic joints

        If w == 0 and v != 0, this corresponds to a translation.
        If w != 0 and v == -w x q, this corresponds to a rotation about the axis centred at q
        If w != 0 and v == 0, this corresponds to a rotation about the origin
        if w == 0 and v == 0, this correponds to the identity element of the group (no transformation)
        if w != 0 and v != 0, this correponds to an arbitrary screw transformation
        """
        w = np.array([[wx, wy, wz]]).T
        v = np.array([[vx, vy, vz]]).T
        mag_w = np.linalg.norm(w)
        mag_v = np.linalg.norm(v)
        self.theta = mag_w if mag_w > tol else mag_v if mag_v > tol else 0
        self.w_hat = w / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T
        self.v_hat = v / self.theta if self.theta > tol else np.array([[0, 0, 0]]).T
        self.rotation = Rotation(wx, wy, wz)

    def __mul__(self, b):
        if isinstance(b, (int, float, complex)) and not isinstance(b, bool):
            new = deepcopy(self)
            new.theta *= b
            return new
    
    __rmul__ = __mul__

    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0, 0, 0, 0, 0, 0)

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
        R = np.eye(3) + sin(self.theta) * self.rotation.so3_norm + \
            (1 - cos(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm)
        p = (np.eye(3) * self.theta + (sin(self.theta) * self.rotation.so3_norm + \
            (1 - cos(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm))) @ self.v_hat
        T = np.vstack((np.hstack((R, p)), np.array([[0, 0, 0, 1]])))
        T[np.abs(T) < tol] = 0.0
        return T

