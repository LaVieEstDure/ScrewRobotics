#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include <iostream>
#include <Eigen/Dense>
#include <representations/transformations.hpp>
#include <representations/aliases.hpp>

template <class Type>
class Transformation{
public:
    Transformation(Type wx, Type wy, Type wz, Type vx, Type vy, Type vz);
    Type theta;
    Vec3<Type> w_hat;
    Vec3<Type> v_hat;
    Rotation<Type> rotation;
    Transformation<Type> identity();
    Transformation<Type> from_array(Type* array);
    Mat6x6<Type> ad();
    Mat6x6<Type> Ad();
    void operator=(Transformation<Type> const& other);
    Mat4x4<Type> se3_norm();
    Mat4x4<Type> se3();
    Mat4x4<Type> SE3();
};


template <class Type>
void Transformation<Type>::operator=(Transformation<Type> const& other) {
    theta = other.theta;
    w_hat = other.w_hat;
    v_hat = other.v_hat;
    rotation = other.rotation;
}


namespace repr {

template <class Type>
Transformation<Type> from_SE3(const Mat4x4<Type> &T);

template <class Type>
Transformation<Type> from_SE3(const Mat4x4<Type> &T) {
    Mat3x3<Type> R = T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Vec3<Type> p = T(Eigen::seq(0, 2), 3);
    Rotation<Type> new_rotation = repr::from_SO3(R);
    Vec3<Type> w_hat = new_rotation.w_hat;
    if (w_hat.norm() < TOL) {
        return Transformation<Type>(0, 0, 0, p(0, 0), p(1, 0), p(2, 0));
    } else {
        Type theta = new_rotation.theta;
        Vec3<Type> w = w_hat * theta;
        Mat3x3<Type> w_mat = new_rotation.so3_norm();
        Mat3x3<Type> G_inv = Mat3x3<Type>::Identity() / theta - w_mat / 2.0 + (1.0 / theta - (1.0 / std::tan(theta / 2.0)) / 2.0) * (w_mat * w_mat);
        Vec3<Type> v_hat = G_inv * p;
        Vec3<Type> v = v_hat * theta;
        return Transformation<Type>(w(0, 0), w(1, 0), w(2, 0), v(0, 0), v(1, 0), v(2, 0));
    }
}

}

template <class Type>
Transformation<Type>::Transformation(Type wx, Type wy, Type wz, Type vx, Type vy, Type vz) {
    static_assert(std::is_arithmetic<Type>::value, "Type must be numeric");
    Vec3<Type> w;
    Vec3<Type> v;
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
        v_hat = v / theta;
    } else {
        w_hat << 0, 0, 0;
        v_hat << 0, 0, 0;
    }
    rotation = Rotation<Type>(wx, wy, wz);
}



template <class Type>
Mat4x4<Type> Transformation<Type>::se3_norm() {
    Mat4x4<Type> V_mat;
    V_mat << rotation.so3_norm(), v_hat;
}


template <class Type>
Mat4x4<Type> Transformation<Type>::se3() {
    return se3_norm() * theta;
}



template <class Type>
Transformation<Type> Transformation<Type>::identity() {
    return Transformation<Type>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
}


template <class Type>
Mat4x4<Type> Transformation<Type>::SE3() {
        Mat3x3<Type> R;
        Mat3x3<Type> G;
        if (rotation.w_hat.norm() < TOL) {
            R = Mat3x3<Type>::Identity();
            G = Mat3x3<Type>::Identity() * theta;
        } else {
            Mat3x3<Type> w_hat_skew = rotation.so3_norm();
            R = Mat3x3<Type>::Identity() + std::sin(theta) * w_hat_skew
                    + (1 - std::cos(theta)) * (w_hat_skew * w_hat_skew);
            G = Mat3x3<Type>::Identity() * theta + (1 - std::cos(theta)) * w_hat_skew
                    + (theta - sin(theta)) * (w_hat_skew * w_hat_skew);
        }
        Vec3<Type> p = G * v_hat;
        Mat4x4<Type> T;
        T << R,       p,
             0, 0, 0, 1;
        return T;
}


// TODO 
template <class Type>
Transformation<Type> Transformation<Type>::from_array(Type* array) {
    return Transformation<Type>(0.0f, 0.0f, 0.0f);
}

template <class Type>
Mat6x6<Type> Transformation<Type>::ad() {
    Vec3<Type> w = w_hat * theta;
    Vec3<Type> v = v_hat * theta;
    Mat3x3<Type> w_bracket = Rotation(w(0, 0), w(1, 0), w(2, 0)).so3();
    Mat3x3<Type> v_bracket = Rotation(v(0, 0), v(1, 0), v(2, 0)).so3();
    Mat6x6<Type> ad_V;
    ad_V << w_bracket, Mat3x3<Type>::Zero(),
            v_bracket,            w_bracket;
    return ad_V;
}


/*


        """
        Returns the lie bracket
        """
        w = self.w_hat * self.theta
        v = self.v_hat * self.theta
        w_bracket = Rotation(w[0, 0], w[1, 0], w[2, 0]).so3
        v_bracket = Rotation(v[0, 0], v[1, 0], v[2, 0]).so3
        ad_V = np.vstack((np.hstack((w_bracket, np.zeros((3, 3)))), np.hstack((v_bracket, w_bracket))))
        return ad_V

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
