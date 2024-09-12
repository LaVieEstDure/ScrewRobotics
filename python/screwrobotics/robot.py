from copy import error
from enum import Enum
from functools import reduce
import operator
import numpy as np
from math import pi
from typing import List, Callable, Optional
from representations import Transformation
import modern_robotics as mr


class Frame(Enum):
    SPACE_FRAME = 0
    BODY_FRAME = 1


class Robot:
    def __init__(self, M_list: List[np.array], screw_list: List[np.array], frame: Frame, 
                ik: Optional[Callable[[np.array], List[float]]] = None, num_ik_iterations=5000):
        assert len(M_list) == len(screw_list) + 2
        self.n = len(screw_list)
        self.M_list = M_list
        if frame == Frame.SPACE_FRAME:
            self.S_list = [*map(lambda x: Transformation(*x), screw_list)]
            self.B_list = np.array([(-1.0 * Transformation.from_SE3(M_list[-1])).Ad @ S.V_hat for S in self.S_list])
            print(self.B_list)
        if frame == Frame.BODY_FRAME:
            self.B_list = [*map(lambda x: Transformation(*x), screw_list)]
        
        self.frame = frame
        self.ik = ik
        self.num_iterations = num_ik_iterations
        
    def forward_kinematics_all_joints(self, theta_list: List[float]) -> List[np.array]:
        S_mats = [theta * self.S_list[i] for (i, theta) in enumerate(theta_list)]
        T = [None] * (self.n + 2)
        T[0] = np.eye(4)
        for i in range(1, self.n + 2):
            index = i if i != self.n + 1 else i - 1
            SE3_mats = map(lambda x: x.SE3, S_mats[:index])
            transformation = reduce(operator.matmul, SE3_mats)
            if self.frame == Frame.SPACE_FRAME:
                T[i] = transformation @ self.M_list[i]
            elif self.frame == Frame.BODY_FRAME:
                T[i] = self.M_list[i] @ transformation
            else:
                raise ValueError("Wrong frame specified!")
        return T

    def jacobian(self, theta_list: List[float]) -> np.array:
        J = np.zeros((6, self.n))
        screw_mats = [theta * self.S_list[i] for (i, theta) in enumerate(theta_list)]
        for i in range(self.n):
            SE3_mats = map(lambda x: x.SE3, screw_mats[:i])
            E = reduce(operator.matmul, SE3_mats, np.eye(4))
            tf = Transformation.from_SE3(E)
            J[:, i:i+1] = tf.Ad @ self.S_list[i].V_hat
        return J
            
    def forward_kinematics(self, theta_list: List[float]) -> np.array:
        return self.forward_kinematics_all_joints(theta_list)[-1]

    def inverse_kinematics(self, Tsd: np.array, theta_list=None) -> List[float]:
        tol = 1e-5
        if self.ik:
            return self.ik(Tsd)
        else:
            if theta_list is None:
                theta_list = np.random.rand(self.n)
            for _ in range(self.num_iterations):
                Tsb = self.forward_kinematics(theta_list)
                trans_bd = Transformation.from_SE3(np.linalg.inv(Tsb) @ Tsd)
                Tbd = trans_bd.SE3
                Ad_Tbd = trans_bd.Ad
                Vb = trans_bd.V_hat * trans_bd.theta
                Vs = Ad_Tbd @ Vb
                error_w = np.linalg.norm([Vs[0,0], Vs[1,0], Vs[2,0]]) 
                error_v = np.linalg.norm([Vs[3,0], Vs[4,0], Vs[5,0]])
                if error_w < tol and error_v < tol:
                    return theta_list
                else:
                    Js = self.jacobian(theta_list)
                    delta_theta_list = np.reshape(np.linalg.pinv(Js) @ Vs, (self.n,))
                    theta_list += delta_theta_list
            return theta_list


if __name__ == "__main__":
    M0 = np.array([[1, 0,  0, 0],
                   [0, 1,  0, 0],
                   [0, 0, 1, 0],
                   [0, 0,  0, 1]])

    M1 = np.array([[1, 0,  0, 0.5],
                   [0, 1,  0, 0],
                   [0, 0, 1, 0],
                   [0, 0,  0, 1]])
    M2 = np.array([[1, 0,  0, 1.5],
                   [0, 1,  0, 0],
                   [0, 0, 1, 0],
                   [0, 0,  0, 1]])
    M3 = np.array([[1, 0,  0, 2.5],
                   [0, 1,  0, 0],
                   [0, 0, 1, 0],
                   [0, 0,  0, 1]])

    M4 = np.array([[1, 0,  0, 3.0],
                   [0, 1,  0, 0],
                   [0, 0, 1, 0],
                   [0, 0,  0, 1]])

    M_list = [M0, M1, M2, M3, M4]

    S1 = np.array([0, 0, 1.0,  0, 0.0, 0])
    S2 = np.array([0, 0, 1.0,  0, -1.0, 0])
    S3 = np.array([0, 0, 1.0,  0, -1.0, 0])

    screw_list = [S1, S2, S3]

    rob = Robot(M_list, screw_list, Frame.SPACE_FRAME)
    theta_list = [np.random.random(), np.random.random(), np.random.random()]
    T = rob.forward_kinematics(theta_list)
    theta_sol = rob.inverse_kinematics(T)
    T_sol = rob.forward_kinematics(theta_sol)
    print(np.round(T, 3))
    print(np.round(T_sol, 3))