from enum import Enum
from functools import reduce
from operator import matmul
from representations import Transformation
import numpy as np
from math import *

class Frame(Enum):
    SPACE_FRAME = 0
    BODY_FRAME = 1

class Robot:
    def __init__(self, M_list, screw_list, frame):
        assert len(M_list) == len(screw_list) + 2
        self.n = len(screw_list)
        self.M_list = M_list
        self.screw_list = [*map(lambda x: Transformation(*x), screw_list)]
        self.frame = frame

    def forward_kinematics(self, theta_list):
        non_norm = [theta * self.screw_list[i] for (i, theta) in enumerate(theta_list)]
        T = [None] * (self.n + 2)
        T[0] = np.eye(4)
        for i in range(1, self.n + 2):
            index = i if i != self.n + 1 else i - 1
            SE3_mats = map(lambda x: x.SE3, non_norm[:index])
            transformation = reduce(matmul, SE3_mats)
            if self.frame == Frame.SPACE_FRAME:
                T[i] = transformation @ self.M_list[i]
            elif self.frame == Frame.BODY_FRAME:
                T[i] = self.M_list[i] @ transformation
            else:
                raise ValueError("Wrong frame specified!")
        return T
                

if __name__ == "__main__":
    M0 = np.array([[1, 0,  0, 0],
                 [ 0, 1,  0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0,  0, 1]])

    M1 = np.array([[1, 0,  0, 0.5],
                 [ 0, 1,  0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0,  0, 1]])
    M2 = np.array([[1, 0,  0, 1.5],
                 [ 0, 1,  0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0,  0, 1]])
    M3 = np.array([[1, 0,  0, 2.5],
                 [ 0, 1,  0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0,  0, 1]])

    M4 = np.array([[1, 0,  0, 3.0],
                 [ 0, 1,  0, 0],
                 [ 0, 0, 1, 0],
                 [ 0, 0,  0, 1]])

    M_list = [M0, M1, M2, M3, M4]
    
    S1 = np.array([0, 0, 1.0,  0, 0.0, 0])
    S2 = np.array([0, 0, 1.0,  0, -1.0, 0])
    S3 = np.array([0, 0, 1.0,  0, -2.0, 0])

    screw_list = [S1, S2, S3]
    theta_list = [0, pi / 2, pi / 2]

    r = Robot(M_list, screw_list, Frame.SPACE_FRAME)

    T = r.forward_kinematics(theta_list)
