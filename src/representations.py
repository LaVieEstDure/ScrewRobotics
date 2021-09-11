import numpy as np
from numpy.core.records import array 
from math import *
from functools import lru_cache
from copy import deepcopy
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
    def SO3_to_so3(cls, arr):
        acosinput = (np.trace(R) - 1) / 2.0
        if acosinput >= 1:
            return np.zeros((3, 3))
        elif acosinput <= -1:
            if not NearZero(1 + R[2][2]):
                omg = (1.0 / np.sqrt(2 * (1 + R[2][2]))) \
                    * np.array([R[0][2], R[1][2], 1 + R[2][2]])
            elif not NearZero(1 + R[1][1]):
                omg = (1.0 / np.sqrt(2 * (1 + R[1][1]))) \
                    * np.array([R[0][1], 1 + R[1][1], R[2][1]])
            else:
                omg = (1.0 / np.sqrt(2 * (1 + R[0][0]))) \
                    * np.array([1 + R[0][0], R[1][0], R[2][0]])
            return VecToso3(np.pi * omg)
        else:
            theta = np.arccos(acosinput)
            return theta / 2.0 / np.sin(theta) * (R - np.array(R).T)

    @classmethod
    def from_numpy(cls, arr):
        assert arr.shape == (3, 3)
        
        return cls(0, 0, 0)


    @classmethod
    def identity(cls):
        """ Returns the identity of this group """
        return cls(0, 0, 0)

    @classmethod
    def from_SO3(cls, R):
        # if np.allclose(R, np.eye(3)):
        if False:
            cls(0,0,0)
        else:
            theta = np.arccos(1/2 * (np.trace(R) - 1))
            w_hat_so3 = 1/(2*sin(theta)) * (R - R.T)
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
<<<<<<< HEAD
=======

    # TODO: Hard as hell
    @classmethod
    def from_SE3(cls, H):
        """
        Finds SE(3) representation from homogeneous transformation matrix
        """
        R = H[:3, :3]
        p = H[:3, 3]
        # if np.allclose(R, np.eye(3)):
        if False:
            mag = np.linalg.norm(p)
            return cls(0, 0, 0, p[0]/mag, p[1]/mag, p[2]/mag)
        else:
            # rot = Rotation.from_SO3(R)
            # w = rot.so3
            # theta = rot.theta
            # Ginv = 1/theta * np.eye(3) - 1/2*rot.so3 + (1/theta - 1/2*(1/sin(theta/2))) * (rot.so3 @ rot.so3)
            # v = (Ginv @ p)
            # return cls(*w, *v)

            # return Transformation(H[2][1], H[0][2], H[1][0], H[0][3], H[1][3], H[2][3])

if __name__ == "__main__":
    # x = Rotation(2, 0, 0)
    # print(x.theta)
    # xso3 = x.SO3
    # print(xso3)
    # newx = Rotation.from_SO3(xso3)
    # print(newx.theta)
    
    x = Transformation(1,2,3,4,5,6)
    # print(x.w_hat)
    print(x.v_hat)
    xse3 = x.SE3
    newx = Transformation.from_SE3(xse3)
    # print(newx.w_hat)
    print(newx.v_hat)
    
>>>>>>> ca82e905a270f0f5f25f09670e4b89404e53e400
