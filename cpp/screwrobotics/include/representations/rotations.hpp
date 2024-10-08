#ifndef ROTATIONS_HPP
#define ROTATIONS_HPP

#include <iostream>
#include <Eigen/Dense>
#include <representations/aliases.hpp>

template <class Type>
class Rotation{
public:
    Rotation(Type wx, Type wy, Type wz);
    Type theta;
    Vec3<Type> w_hat;
    Rotation<Type> identity();
    Rotation<Type>();
    Rotation<Type> from_array(Type* array);
    Mat3x3<Type> ad(Type* array);
    Mat3x3<Type> Ad(Type* array);
    void operator=(Rotation<Type> const& other);
    Mat3x3<Type> so3_norm();
    Mat3x3<Type> so3();
    Mat3x3<Type> SO3();
    
};

namespace repr {

template <class Type>
Rotation<Type> from_SO3(const Mat3x3<Type> &R);

template <class Type>
Rotation<Type> from_SO3(const Mat3x3<Type> &R) {
    if (std::abs((R - Mat3x3<Type>::Identity()).norm()) < TOL) {
        return Rotation<Type>(0, 0, 0);
    }
    if (std::abs(R.trace()) + 1.0 < TOL) {
        Type wx = (M_PI / std::sqrt(2 * (1 + R(2, 2)))) * R(0, 2);
        Type wy = (M_PI / std::sqrt(2 * (1 + R(2, 2)))) * R(1, 2);
        Type wz = (M_PI / std::sqrt(2 * (1 + R(2, 2)))) * (1.0 + R(2, 2));
        return Rotation<Type>(wx, wy, wz);
    }
    Type theta = std::acos((R.trace() - 1.0) / 2.0);
    Mat3x3<Type> w_hat_so3 = (R - R.transpose()) / (2 * std::sin(theta));
    Type wx = theta * w_hat_so3(2, 1);
    Type wy = theta * w_hat_so3(0, 2);
    Type wz = theta * w_hat_so3(1, 0);
    return Rotation<Type>(wx, wy, wz);
}


}
template <class Type>
Rotation<Type>::Rotation(Type wx, Type wy, Type wz) {
    static_assert(std::is_arithmetic<Type>::value, "Type must be numeric");
    w_hat << wx, wy, wz;
    theta = w_hat.norm();
    if (theta > TOL) {
        w_hat = w_hat / theta;
    }
    else {
        theta = 0;
        w_hat << 0, 0, 0;
    }
}


template <class Type>
Rotation<Type>::Rotation() {

}

template <class Type>
void Rotation<Type>::operator=(Rotation<Type> const& other) {
    theta = other.theta;
    w_hat = other.w_hat;
}

template <class Type>
Mat3x3<Type> Rotation<Type>::so3_norm() {
    Mat3x3<Type> w_hat_skew;
    w_hat_skew <<           0, -w_hat(2, 0),  w_hat(1, 0),
                  w_hat(2, 0),            0, -w_hat(0, 0),
                  w_hat(1, 0),  w_hat(0, 0),            0;
    return w_hat_skew;
}

template <class Type>
Mat3x3<Type> Rotation<Type>::so3() {
    return so3_norm() * theta;
}

template <class Type>
Mat3x3<Type> Rotation<Type>::SO3() {
    Mat3x3<Type> w_hat_skew = so3_norm();
    Mat3x3<Type> R = Mat3x3<Type>::Identity() + std::sin(theta) * w_hat_skew + (1 - std::cos(theta)) * (w_hat_skew * w_hat_skew);
    return R;
}

/* TODO */
template <class Type>
Rotation<Type> Rotation<Type>::from_array(Type* array) {
    return Rotation<Type>(0.0f, 0.0f, 0.0f);
}


template <class Type>
Rotation<Type> Rotation<Type>::identity() {
    return Rotation<Type>(0.0f, 0.0f, 0.0f);
}

template <class Type>
Eigen::Matrix <Type, 3, 3> Rotation<Type>::ad(Type* array) {
    return Rotation<Type>(0.0f, 0.0f, 0.0f).so3();
}


template <class Type>
Eigen::Matrix <Type, 3, 3> Rotation<Type>::Ad(Type* array) {
    return Rotation<Type>(0.0f, 0.0f, 0.0f).SO3();
}

#endif