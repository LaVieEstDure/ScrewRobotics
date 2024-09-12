#ifndef ROBOT_HPP
#define ROBOT_HPP

#include <iostream>
#include <Eigen/Dense>
#include <representations/transformations.hpp>
#include <representations/rotations.hpp>
#include <representations/aliases.hpp>
#include <vector>

class Robot {
private:
    std::vector<Mat4x4> M_list;
    std::vector<Mat4x4> screw_list;
public:
    Robot(std::vector<Mat4x4> _M_list, std::vector<Vec6> _screw_list) : M_list(_screw_list), screw_list(_screw_list) {
        n = _screw_list.size();
    }
}

/*

class Robot:
    def __init__(self, M_list: List[np.array], screw_list: List[np.array], frame: Frame, 
                ik: Optional[Callable[[np.array], List[float]]] = None, num_ik_iterations=20):
        assert len(M_list) == len(screw_list) + 2
        self.n = len(screw_list)
        self.M_list = M_list
        self.screw_list = [*map(lambda x: Transformation(*x), screw_list)]
        self.num_joints = len(self.screw_list)
        self.frame = frame
        self.ik = ik
        self.num_iterations = num_ik_iterations

    def forward_kinematics_all_joints(self, theta_list: List[float]) -> List[np.array]:
        non_norm = [theta * self.screw_list[i]
                    for (i, theta) in enumerate(theta_list)]
        T = [None] * (self.n + 2)
        T[0] = np.eye(4)
        for i in range(1, self.n + 2):
            index = i if i != self.n + 1 else i - 1
            SE3_mats = map(lambda x: x.SE3, non_norm[:index])
            transformation = reduce(operator.matmul, SE3_mats)
            if self.frame == Frame.SPACE_FRAME:
                T[i] = transformation @ self.M_list[i]
            elif self.frame == Frame.BODY_FRAME:
                T[i] = self.M_list[i] @ transformation
            else:
                raise ValueError("Wrong frame specified!")
        return T

    def forward_kinematics(self, theta_list: List[float]) -> np.array:
        return self.forward_kinematics_all_joints(theta_list)[-1]

    def inverse_kinematics(self, pose: np.array) -> List[float]:
        tol = 1e-8
        if self.ik:
            return self.ik(pose)
        else:
            joint_values = np.random.rand(self.num_joints)
            end_effector_frame = self.forward_kinematics(joint_values)
            Tsbinv_Tsd = Transformation.from_SE3((-Transformation.from_SE3(pose)).SE3 @ end_effector_frame)
            Vs = Transformation.from_SE3(end_effector_frame).Ad @ Tsbinv_Tsd.s
            error_w = np.linalg.norm([Vs[0,0], Vs[1,0], Vs[2,0]]) 
            error_v = np.linalg.norm([Vs[3,0], Vs[4,0], Vs[5,0]])
            for i in range(self.num_iterations):
                if error_w < tol and error_v < tol:
                    return joint_values
                break


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
    S3 = np.array([0, 0, 1.0,  0, -2.0, 0])

    screw_list = [S1, S2, S3]
    theta_list = [0, pi / 2, pi / 2]

    r = Robot(M_list, screw_list, Frame.SPACE_FRAME)
    fk = r.forward_kinematics([0.1, 0.2, 0.3])
    r.inverse_kinematics(fk)

*/

#endif