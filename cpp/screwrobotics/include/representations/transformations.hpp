#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP

#include <iostream>
#include <Eigen/Dense>
#include <representations/transformations.hpp>

template <class Type>
class Transformation{
public:
    Transformation(Type wx, Type wy, Type wz, Type vx, Type vy, Type vz);
    Type theta;
    Eigen::Matrix<Type, 3, 1> w_hat;
    Eigen::Matrix<Type, 3, 1> v_hat;
    Rotation<Type> rotation;
    Transformation<Type> identity();
    Transformation<Type> from_array(Type* array);
    Eigen::Matrix<Type, 6, 6> ad();
    Eigen::Matrix<Type, 6, 6> Ad();
    void operator=(Transformation<Type> const& other);
    Eigen::Matrix<Type, 4, 4> se3_norm();
    Eigen::Matrix<Type, 4, 4> se3();
    Eigen::Matrix<Type, 4, 4> SE3();
};


template <class Type>
void Transformation<Type>::operator=(Transformation<Type> const& other) {
    theta = other.theta;
    w_hat = other.w_hat;
    v_hat = other.v_hat;
    rotation = other.rotation;
}


namespace trans {

template <class Type>
Transformation<Type> from_SE3(const Eigen::Matrix<Type, 4, 4> &T);

template <class Type>
Transformation<Type> from_SE3(const Eigen::Matrix<Type, 4, 4> &T) {
    Eigen::Matrix<Type, 3, 3> R = T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Eigen::Matrix<Type, 3, 1> p = T(Eigen::seq(0, 2), 3);
    Rotation<Type> new_rotation = rot::from_SO3(R);
    Eigen::Matrix<Type, 3, 1> w_hat = new_rotation.w_hat;
    if (w_hat.norm() < TOL) {
        return Transformation<Type>(0, 0, 0, p(0, 0), p(1, 0), p(2, 0));
    } else {
        Type theta = new_rotation.theta;
        Eigen::Matrix<Type, 3, 1> w = w_hat * theta;
        Eigen::Matrix<Type, 3, 3> w_mat = new_rotation.so3_norm();
        Eigen::Matrix<Type, 3, 3> G_inv = Eigen::Matrix<Type, 3, 3>::Identity() / theta - w_mat / 2.0 + (1.0 / theta - (1.0 / std::tan(theta / 2.0)) / 2.0) * (w_mat * w_mat);
        Eigen::Matrix<Type, 3, 1> v_hat = G_inv * p;
        Eigen::Matrix<Type, 3, 1> v = v_hat * theta;
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
        v_hat = v / theta;
    } else {
        w_hat << 0, 0, 0;
        v_hat << 0, 0, 0;
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



template <class Type>
Transformation<Type> Transformation<Type>::identity() {
    return Transformation<Type>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
}


template <class Type>
Eigen::Matrix<Type, 4, 4> Transformation<Type>::SE3() {
        Eigen::Matrix<Type, 3, 3> R;
        Eigen::Matrix<Type, 3, 3> G;
        if (rotation.w_hat.norm() < TOL) {
            R = Eigen::Matrix<Type, 3, 3>::Identity();
            G = Eigen::Matrix<Type, 3, 3>::Identity() * theta;
        } else {
            Eigen::Matrix<Type, 3, 3> w_hat_skew = rotation.so3_norm();
            R = Eigen::Matrix<Type, 3, 3>::Identity() + std::sin(theta) * w_hat_skew
                    + (1 - std::cos(theta)) * (w_hat_skew * w_hat_skew);
            G = Eigen::Matrix<Type, 3, 3>::Identity() * theta + (1 - std::cos(theta)) * w_hat_skew
                    + (theta - sin(theta)) * (w_hat_skew * w_hat_skew);
        }
        Eigen::Matrix<Type, 3, 1> p = G * v_hat;
        Eigen::Matrix<Type, 4, 4> T;
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
Eigen::Matrix <Type, 6, 6> Transformation<Type>::ad() {
    Eigen::Matrix<Type, 3, 1> w = w_hat * theta;
    Eigen::Matrix<Type, 3, 1> v = v_hat * theta;
    Eigen::Matrix<Type, 3, 3> w_bracket = Rotation(w(0, 0), w(1, 0), w(2, 0)).so3();
    Eigen::Matrix<Type, 3, 3> v_bracket = Rotation(v(0, 0), v(1, 0), v(2, 0)).so3();
    Eigen::Matrix<Type, 6, 6> ad_V;
    ad_V << w_bracket, Eigen::Matrix<Type, 3, 3>::Zero(),
            v_bracket,                          w_bracket;
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