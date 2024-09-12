#ifndef ALIASES_HPP
#define ALIASES_HPP


template <class Type>
using Vec2 = Eigen::Matrix<Type, 2, 1>;

template <class Type>
using Vec3 = Eigen::Matrix<Type, 3, 1>;

template <class Type>
using Vec4 = Eigen::Matrix<Type, 4, 1>;

template <class Type>
using Vec6 = Eigen::Matrix<Type, 6, 1>;

template <class Type>
using Mat2x2 = Eigen::Matrix<Type, 2, 2>;

template <class Type>
using Mat3x3 = Eigen::Matrix<Type, 3, 3>;

template <class Type>
using Mat4x4 = Eigen::Matrix<Type, 4, 4>;

template <class Type>
using Mat6x6 = Eigen::Matrix<Type, 4, 4>;


#endif