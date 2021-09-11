from enum import Enum
from functools import reduce
from operator import matmul
from representations import Transformation
import numpy as np

class Frame(Enum):
    SPACE_FRAME = 0
    BODY_FRAME = 1

class Robot:
    def __init__(self, M, screw_list, frame):
        self.M = M
        self.screw_list = [*map(lambda x: Transformation(*x), screw_list)]
        self.frame = frame

    def forward_kinematics(self, theta_list):
        non_norm = [theta*self.screw_list[i] for i, theta in enumerate(theta_list)]
        se3mats = map(lambda x: x.SE3, non_norm)
        transformation = reduce(matmul, se3mats)
        if self.frame == Frame.SPACE_FRAME:
            return transformation @ self.M
        elif self == Frame.BODY_FRAME:
            return self.M @ transformation
        else:
            raise ValueError("Wrong frame specified!")

if __name__ == "__main__":
    M = np.array([[-1, 0,  0, 0],
                 [ 0, 1,  0, 6],
                 [ 0, 0, -1, 2],
                 [ 0, 0,  0, 1]])

    screw_vecs = np.array([[0, 0,  1,  4, 0,    0],
                           [0, 0,  0,  0, 1,    0],
                           [0, 0, -1, -6, 0, -0.1]])

    joint_states = np.array([np.pi / 2.0, 3, np.pi])

    r = Robot(M, screw_vecs, Frame.SPACE_FRAME)

    end_effector_loc = r.forward_kinematics(joint_states)

