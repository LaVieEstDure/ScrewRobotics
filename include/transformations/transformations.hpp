#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include <iostream>
#include <Eigen/Dense>
#include <transformations/transformations.hpp>
#define TOL 1e-6

template <class Type>
class Transformation{
public:
    Transformation(Type wx, Type wy, Type wz, Type vx, Type vy, Type vz);
    Type theta;
    Eigen::Matrix<Type, 3, 1> w_hat;
    Eigen::Matrix<Type, 3, 1> v_hat;
    Rotation<Type> rotation;
    //Transformation<Type> identity();
    //Transformation<Type> from_array(Type* array);
    //Eigen::Matrix<Type, 6, 6> ad(Type* array);
    //Eigen::Matrix<Type, 6, 6> Ad(Type* array);
    //Transformation<Type>& operator=(const Transformation<Type>& other);
    Eigen::Matrix<Type, 4, 4> se3_norm();
    Eigen::Matrix<Type, 4, 4> se3();
    //Eigen::Matrix<Type, 4, 4> SE3();
};


namespace trans {

template <class Type>
Transformation<Type> from_SE3(const Eigen::Matrix<Type, 4, 4> &T);

template <class Type>
Transformation<Type> from_SE3(const Eigen::Matrix<Type, 4, 4> &T) {
    Eigen::Matrix<Type, 3, 3> R = T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Eigen::Matrix<Type, 3, 1> p = T(Eigen::seq(0, 2), 3);
    Rotation<Type> rotation = rot::from_SO3(R);
    Eigen::Matrix<Type, 3, 1> w_hat = rotation.w_hat;
    if (w_hat.norm() < TOL) {
        return Transformation<Type>(0, 0, 0, p(0, 0), p(1, 0), p(2, 0));
    } else {
        Type theta = rotation.theta;
        Eigen::Matrix<Type, 3, 1> w = w_hat * theta;
        Eigen::Matrix<Type, 3, 3> w_mat = rotation.so3_norm();
        Eigen::Matrix<Type, 3, 3> G_inv = Eigen::Matrix<Type, 3, 3>::Identity() / theta - w_mat / 2.0 + (1.0 / theta - (1.0 / std::tan(theta / 2.0)) / 2.0) * (w_mat * w_mat);
        Eigen::Matrix<Type, 3, 1> v_hat = G_inv * p;
        Eigen::Matrix<Type, 3, 1> v = v_hat * theta;
        std::cout << v << std::endl;
        return Transformation<Type>(w(0, 0), w(1, 0), w(2, 0), v(0, 0), v(1, 0), v(2, 0));
    }
}

}

template <class Type>
Transformation<Type>::Transformation(Type wx, Type wy, Type wz, Type vx, Type vy, Type vz) {
    static_assert(std::is_arithmetic<Type>::value, "Type must be numeric");
    Eigen::Matrix<Type, 3, 1> w;
    Eigen::Matrix<Type, 3, 1> v;
    w << wx, wy, wz;
    v << vx, vy, vz;
    Type mag_w = w.norm();
    Type mag_v = v.norm();
    if (mag_w > TOL) {
        theta = mag_w;
    } else if (mag_v > TOL) {
        theta = mag_v;
    } else {
        theta = 0;
    }
    if (theta > TOL) {
        w_hat = w / theta;
    } else {
        w_hat << 0, 0, 0;
    }
    rotation = Rotation<Type>(wx, wy, wz);
}



template <class Type>
Eigen::Matrix<Type, 4, 4> Transformation<Type>::se3_norm() {
    Eigen::Matrix<Type, 4, 4> V_mat;
    V_mat << rotation.so3_norm(), v_hat;
}


template <class Type>
Eigen::Matrix<Type, 4, 4> Transformation<Type>::se3() {
    return se3_norm() * theta;
}
/*

template <class Type>
Eigen::Matrix<Type, 3, 3> Transformation<Type>::SE3() {
    Eigen::Matrix<Type, 3, 3> V_hat_skew = se3_norm();
    Eigen::Matrix<Type, 3, 3> R = Eigen::Matrix<Type, 3, 3>::Identity() + std::sin(theta) * w_hat_skew + (1 - std::cos(theta)) * (w_hat_skew * w_hat_skew);
    return R;
}

// TODO 
template <class Type>
Transformation<Type> Transformation<Type>::from_array(Type* array) {
    return Rotation<Type>(0.0f, 0.0f, 0.0f);
}


template <class Type>
Transformation<Type> Transformation<Type>::identity() {
    return Rotation<Type>(0.0f, 0.0f, 0.0f);
}

template <class Type>
Eigen::Matrix <Type, 6, 6> Transformation<Type>::ad(Type* array) {
    return Transformation<Type>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f).se3();
}


template <class Type>
Eigen::Matrix <Type, 6, 6> Transformation<Type>::Ad(Type* array) {
    return Transformation<Type>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f).SE3();
}
*/
#endif