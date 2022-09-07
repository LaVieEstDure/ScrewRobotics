#include <iostream>
#include <Eigen/Dense>
#include <rotations/rotations.hpp>
#include <transformations/transformations.hpp>

int main() {
    Eigen::Matrix<float, 4, 4> T;
    T << 1 / sqrt(2), 1 / sqrt(2), 0, 1,
         -1 / sqrt(2), 1 / sqrt(2), 0, 2,
         0, 0, 1, -0.5,
         0, 0, 0, 1;
    auto tf = trans::from_SE3(T);

}