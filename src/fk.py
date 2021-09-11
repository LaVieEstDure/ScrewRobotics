from enum import Enum
from functools import reduce
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
        non_normalised = [theta*self.screw_list[i] for i, theta in enumerate(theta_list)]
        transformation = reduce(lambda t1, t2: t1.SE3 @ t2.SE3, non_normalised)
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
    Slist = np.array([[0, 0,  1,  4, 0,    0],
                      [0, 0,  0,  0, 1,    0],
                      [0, 0, -1, -6, 0, -0.1]])
    thetalist = np.array([np.pi / 2.0, 3, np.pi])

    r = Robot(M, Slist, Frame.SPACE_FRAME)
    print(r.forward_kinematics(thetalist))