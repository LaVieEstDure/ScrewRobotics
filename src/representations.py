import numpy as np
from numpy.core.records import array 
from math import *
from functools import lru_cache
from copy import deepcopy

import matplotlib.pyplot as plt

import modern_robotics as mr

# Tolerance
tol = 1e-10

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
    def from_numpy(cls, arr):
        assert arr.shape == (3, 3)
        
        return cls(0, 0, 0)


    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0, 0, 0)

    @property
    @lru_cache(maxsize=1)
    def ad(self):
        """
        Returns the lie bracket (skew matrix)
        """
        w = self.w_hat * self.theta
        Ad_w = Rotation(w[0, 0], w[1, 0], w[2, 0]).so3
        return Ad_w

    @property
    @lru_cache(maxsize=1)
    def Ad(self):
        """
        Returns the adjoint map associated with SO3
        """
        Ad_R = self.SO3
        return Ad_R

    @classmethod
    def from_SO3(cls, R):
        if abs(np.linalg.norm(R - np.ones(3))) < tol:
            return cls(0, 0, 0)
        if abs(np.trace(R) + 1.0) < tol:
            wx = (pi / sqrt(2 * (1 + R[2][2]))) * R[0,2]
            wy = (pi / sqrt(2 * (1 + R[2][2]))) * R[1,2]
            wz = (pi / sqrt(2 * (1 + R[2][2]))) * (1.0 + R[2,2])
            return cls(wx, wy, wz)
        if False:
            cls(0,0,0)
        else:
            theta = acos((np.trace(R) - 1.0) / 2.0)
            w_hat_so3 = (R - R.T) / (2 * sin(theta))
            w_hat = np.array([w_hat_so3[2,1], w_hat_so3[0, 2], w_hat_so3[1, 0]])
            return cls(*(w_hat*theta))

    @property
    @lru_cache(maxsize=1)
    def so3_norm(self):
        """Returns the transformation in so3 (element of Lie Algebra), when theta = 1"""
        return np.array([[             0, -self.w_hat[2,0],  self.w_hat[1,0]], 
                         [ self.w_hat[2,0],              0, -self.w_hat[0,0]], 
                         [-self.w_hat[1,0],  self.w_hat[0,0],              0]])

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

    @classmethod
    def from_numpy(cls, arr):
        assert arr.shape == (4, 4)
        return cls(0, 0, 0, 0, 0, 0)

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
        p = (np.eye(3) * self.theta + (1 - cos(self.theta)) * self.rotation.so3_norm + \
            (self.theta - sin(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm)) @ self.v_hat
        T = np.vstack((np.hstack((R, p)), np.array([[0, 0, 0, 1]])))
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
        w_bracket = Rotation(w[0, 0], w[1, 0], w[2, 0]).so3
        v_bracket = Rotation(v[0, 0], v[1, 0], v[2, 0]).so3
        ad_V = np.vstack((np.hstack((w_bracket, np.zeros((3, 3)))), np.hstack((v_bracket, w_bracket))))
        return ad_V

    @property
    @lru_cache(maxsize=1)
    def Ad(self):
        """
        Returns the adjoint map associated with SE3
        """
        T = self.SE3
        R = T[:3, :3]
        p = np.array([T[:3, 3]]).T
        p_rot = Rotation(p[0, 0], p[1, 0], p[2, 0])
        Ad_T = np.vstack((np.hstack((R, np.zeros((3, 3)))), np.hstack((p_rot.so3 @ R, R))))
        return Ad_T

    @classmethod
    def from_SE3(cls, T):
        """
        Finds SE(3) representation from homogeneous transformation matrix
        """
        R = T[:3, :3]
        p = np.array([T[:3, 3]]).T
        rotation = Rotation.from_SO3(R)
        w_hat = rotation.w_hat
        if np.linalg.norm(w_hat) < tol:
            return cls(0, 0, 0, p[0, 0], p[1, 0], p[2, 0])
        else:
            theta = rotation.theta
            w = w_hat * theta
            w_mat = rotation.so3_norm
            v_hat = (np.eye(3) / theta - w_mat / 2.0 + (1.0 / theta - (1.0 / np.tan(theta / 2.0)) / 2.0) * (w_mat @ w_mat)) @ p
            v = v_hat * theta
            tf = cls(w[0, 0], w[1, 0], w[2, 0], v[0, 0], v[1, 0], v[2, 0])
            return tf


class Transformation2D:
    def __init__(self, w, wz, vx, vy):
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
        w = np.array([[wx, wy, wz]]).T
        v = np.array([[vx, vy, vz]]).T
        mag_w = np.linalg.norm(w)
        mag_v = np.linalg.norm(v)
        self.theta = mag_w if mag_w > tol else mag_v if mag_v > tol else 0
        self.w_hat = w / self.theta if self.theta > tol else 0
        self.v_hat = v / self.theta if self.theta > tol else np.array([[0, 0]]).T
        self.rotation = Rotation(wx, wy, wz)

    @classmethod
    def from_numpy(cls, arr):
        assert arr.shape == (4, 4)
        return cls(0, 0, 0, 0, 0, 0)

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
        p = (np.eye(3) * self.theta + (1 - cos(self.theta)) * self.rotation.so3_norm + \
            (self.theta - sin(self.theta)) * (self.rotation.so3_norm @ self.rotation.so3_norm)) @ self.v_hat
        T = np.vstack((np.hstack((R, p)), np.array([[0, 0, 0, 1]])))
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
        w_bracket = Rotation(w[0, 0], w[1, 0], w[2, 0]).so3
        v_bracket = Rotation(v[0, 0], v[1, 0], v[2, 0]).so3
        ad_V = np.vstack((np.hstack((w_bracket, np.zeros((3, 3)))), np.hstack((v_bracket, w_bracket))))
        return ad_V

    @property
    @lru_cache(maxsize=1)
    def Ad(self):
        """
        Returns the adjoint map associated with SE3
        """
        T = self.SE3
        R = T[:3, :3]
        p = np.array([T[:3, 3]]).T
        p_rot = Rotation(p[0, 0], p[1, 0], p[2, 0])
        Ad_T = np.vstack((np.hstack((R, np.zeros((3, 3)))), np.hstack((p_rot.so3 @ R, R))))
        return Ad_T

    @classmethod
    def from_SE3(cls, T):
        """
        Finds SE(3) representation from homogeneous transformation matrix
        """
        R = T[:3, :3]
        p = np.array([T[:3, 3]]).T
        rotation = Rotation.from_SO3(R)
        w_hat = rotation.w_hat
        if np.linalg.norm(w_hat) < tol:
            return cls(0, 0, 0, p[0, 0], p[1, 0], p[2, 0])
        else:
            theta = rotation.theta
            w = w_hat * theta
            w_mat = rotation.so3_norm
            v_hat = (np.eye(3) / theta - w_mat / 2.0 + (1.0 / theta - (1.0 / np.tan(theta / 2.0)) / 2.0) * (w_mat @ w_mat)) @ p
            v = v_hat * theta
            tf = cls(w[0, 0], w[1, 0], w[2, 0], v[0, 0], v[1, 0], v[2, 0])
            return tf

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

def screw_demo():
    ax = plt.axes(projection='3d')
    tf = Transformation(0, 3/sqrt(5), 4/sqrt(5), 0.0, -1.0, 2.0)
    x = []
    y = []
    z = []
    ix_vec = []
    jx_vec = []
    kx_vec = []
    iy_vec = []
    jy_vec = []
    ky_vec = []
    iz_vec = []
    jz_vec = []
    kz_vec = []
    for i in range(100):
        T = ((2 * pi * i / 50.0) * tf).SE3 @ Transformation(0, 0, 0, 2, 0, 0).SE3
        x.append(T[0, 3])
        y.append(T[1, 3])
        z.append(T[2, 3])
        ix_vec.append(T[0, 0])
        iy_vec.append(T[1, 0])
        iz_vec.append(T[2, 0])
        jx_vec.append(T[0, 1])
        jy_vec.append(T[1, 1])
        jz_vec.append(T[2, 1])
        kx_vec.append(T[0, 2])
        ky_vec.append(T[1, 2])
        kz_vec.append(T[2, 2])
    ax.scatter3D(x, y, z, color='k')
    ax.quiver3D(x, y, z, ix_vec, iy_vec, iz_vec, color='r', length=0.3, arrow_length_ratio=0.1)
    ax.quiver3D(x, y, z, jx_vec, jy_vec, jz_vec, color='g', length=0.3, arrow_length_ratio=0.1)
    ax.quiver3D(x, y, z, kx_vec, ky_vec, kz_vec, color='b', length=0.3, arrow_length_ratio=0.1)
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([0, 10])
    axisEqual3D(ax)
    plt.show()
    

if __name__ == "__main__":

    V = Transformation(0, 0, 1, 0, -1, 0)
    V_mr = np.array([0, 0, 1, 0, -1, 0])
    print(V.ad)
    print(mr.ad(V_mr))