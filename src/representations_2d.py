import numpy as np
from numpy.core.records import array 
from math import *
from functools import lru_cache
from copy import deepcopy

import matplotlib.pyplot as plt

import modern_robotics as mr

# Tolerance
tol = 1e-10

class Rotation2D:
    def __init__(self, w):
        """
        Defines a rotation from its axis

        The rotation is given by an axis of rotation w_hat and angle of rotation theta
        
        The euler vector is decomposed as w = w_hat * theta,
        where theta is the angle of rotation
        and w_hat is the axis of rotation

        If w = 0, the axis is undefined and the angle is zero
        This corresponds to the identity rotation matrix.
        """
        self.theta = abs(w)
        self.w_hat = w / self.theta if self.theta > tol else 0

    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0)

    @property
    @lru_cache(maxsize=1)
    def ad(self):
        """
        Returns the lie bracket (skew matrix)
        """
        w = self.w_hat * self.theta
        Ad_w = Rotation2D(w).so2
        return Ad_w

    @property
    @lru_cache(maxsize=1)
    def Ad(self):
        """
        Returns the adjoint map associated with SO2
        """
        Ad_R = self.SO2
        return Ad_R

    @classmethod
    def from_SO2(cls, R):
        return cls(atan2(R[1, 0], R[0, 0]))

    @property
    @lru_cache(maxsize=1)
    def so2_norm(self):
        """Returns the transformation in so2 (element of Lie Algebra), when theta = 1"""
        return np.array([
            [0, -self.w_hat], 
            [self.w_hat, 0]
        ])

    @property
    @lru_cache(maxsize=1)
    def so2(self):
        """Returns the transformation in so2 (element of Lie Algebra)"""
        return self.theta * self.so2_norm

    @property
    @lru_cache(maxsize=1)
    def SO2(self):
        """Returns the transformation in SO2 (element of Lie Group)"""
        R = cos(self.theta) * np.eye(2) + sin(self.theta) * self.so2_norm
        R[np.abs(R) < tol] = 0.0
        return R

class Transformation2D:
    def __init__(self, w, vx, vy):
        """
        Defines a transformation (rotation + translation) from its Screw Vector
        w = [w vx vy]^T

        The rotation is given by an axis of rotation w_hat and angle of rotation theta
        The translation has two terms dueto the linear motion at the origin induced by
        rotation about the axis, and the other due to the translation along the screw axis
        
        The screw vector is decomposed as V = S_hat * theta,
        where theta is the angle of rotation for revolute joints
        or the distance of translation for prismatic joints

        If w == 0 and v != 0, this corresponds to a translation.
        If w != 0 and v == -[w] q, this corresponds to a rotation about the axis centred at q
        If w != 0 and v == 0, this corresponds to a rotation about the origin
        if w == 0 and v == 0, this correponds to the identity element of the group (no transformation)
        if w != 0 and v != 0, this correponds to an arbitrary screw transformation
        """
        v = np.array([[vx, vy]]).T
        mag_w = abs(w)
        mag_v = np.linalg.norm(v)
        self.theta = mag_w if mag_w > tol else mag_v if mag_v > tol else 0
        self.w_hat = w / self.theta if self.theta > tol else 0
        self.v_hat = v / self.theta if self.theta > tol else np.array([[0, 0]]).T
        self.rotation = Rotation2D(w)

    def __mul__(self, b):
        if isinstance(b, (int, float, complex)) and not isinstance(b, bool):
            new = deepcopy(self)
            new.theta *= b
            return new
    
    __rmul__ = __mul__

    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0, 0, 0)

    @property
    @lru_cache(maxsize=1)
    def se2_norm(self):
        return np.vstack((np.hstack((self.rotation.so2_norm, self.v_hat)), np.array([[0, 0, 0]])))
    
    '''
    Returns the transformation in so2 (Lie Algebra)
    '''
    @property
    @lru_cache(maxsize=1)
    def se2(self):
        return self.theta * self.se2_norm

    @property
    @lru_cache(maxsize=1)
    def SE2(self):
        if abs(self.rotation.w_hat) < tol:
            R = np.eye(2)
            G = self.theta * np.eye(2)
        else:
            R = cos(self.theta) * np.eye(2) + sin(self.theta) * self.rotation.so2_norm
            G = np.eye(2) * sin(self.theta) + (1 - cos(self.theta)) * self.rotation.so2_norm
        p = G @ self.v_hat
        T = np.vstack((np.hstack((R, p)), np.array([[0, 0, 1]])))
        T[np.abs(T) < tol] = 0.0
        return T



    @property
    @lru_cache(maxsize=1)
    def ad(self):
        """
        Returns the lie bracket
        """
        w = self.w_hat * self.theta
        v = self.v_hat * self.theta
        w_bracket = Rotation2D(w[0, 0], w[1, 0], w[2, 0]).so2
        v_bracket = Rotation2D(v[0, 0], v[1, 0], v[2, 0]).so2
        ad_V = np.vstack((np.hstack((w_bracket, np.zeros((3, 3)))), np.hstack((v_bracket, w_bracket))))
        return ad_V

    @property
    @lru_cache(maxsize=1)
    def Ad(self):
        """
        Returns the adjoint map associated with SE(2)
        """
        T = self.SE2
        R = T[:2, :2]
        p = np.array([T[:2, 2]]).T
        skew = np.array([[0, -1], [1, 0]])
        Ad_T = np.vstack((np.array([[1, 0, 0]]), np.hstack((-skew @ p, R))))
        return Ad_T

    @classmethod
    def from_SE2(cls, T):
        """
        Finds SE(2) representation from homogeneous transformation matrix
        """
        R = T[:2, :2]
        p = np.array([T[:2, 2]]).T
        rotation = Rotation2D.from_SO2(R)
        w_hat = rotation.w_hat
        if abs(w_hat) < tol:
            return cls(0, p[0, 0], p[1, 0])
        else:
            theta = rotation.theta
            w = w_hat * theta
            w_mat = rotation.so2_norm
            G_inv = sin(theta) / (2.0 * (1.0 - cos(theta))) - rotation.so2_norm / 2.0
            v_hat = G_inv @ p
            v = v_hat * theta
            tf = cls(w, v[0, 0], v[1, 0])
            return tf

if __name__ == "__main__":
    T = np.array([
        [1, 0, 3],
        [0, 1, 0],
        [0, 0, 1]
    ])
    tf = Transformation2D.from_SE2(T)
    print(tf.SE2)