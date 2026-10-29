#include "rigid_motion/frame.h"
#include "rigid_motion/screw.h"

#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>

namespace ais4104::rigid_motion {

//TASK: 1d
//REFERENCE: Section 3.2.1.2, page 72, MR 2019
Eigen::Matrix3d rotate_x(double radians)
{
    Eigen::Matrix3d R;
    R<< 1, 0, 0,
    0, std::cos(radians) , -std::sin(radians) ,
    0, std::sin(radians)  , std::cos(radians);
    return R;
}

//TASK: 1e
//REFERENCE:  Section 3.2.1.2, page 72, MR 2019
Eigen::Matrix3d rotate_y(double radians)
{
    Eigen::Matrix3d R;
    R<< std::cos(radians), 0, std::sin(radians),
    0, 1 , 0 ,
    -std::sin(radians), 0  , std::cos(radians);
    return R;
}

//TASK: 1f
//REFERENCE:Section 3.2.1.2, page 72, MR 2019
Eigen::Matrix3d rotate_z(double radians)
{
    Eigen::Matrix3d R;
    R<< std::cos(radians) , -std::sin(radians), 0,
    std::sin(radians)  , std::cos(radians) , 0 ,
    0, 0 , 1;
    return R;
}

//TASK: 1g
//REFERENCE: Section 3.2.1.2, Figure 3.7, page 70, MR2019
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
{
    Eigen::Matrix3d ReferenceR; ReferenceR<<
        x(0) , y(0), z(0),
    x(1)  ,y(1) , z(1) ,
    x(2), y(2) , z(2);
    return ReferenceR;
}

//TASK: 1h
//REFERENCE: Section 3.2.3.2, Equation (3.51), page 82, MR2019
Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
{
    Eigen::Vector3d w = axis.normalized();
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d W = skew_symmetric(w);
    Eigen::Matrix3d R = I + std::sin(radians) * W + (1 - std::cos(radians)) * (W * W);
    return R;
}

//TASK: 1i Implement the special cases for Euler ZYX and XYZ using rotate_x, rotate_y, and rotate_z
//REFERENCE: Appendix B.1, Euler Angles, pp. 577–580, MR 2019
Eigen::Matrix3d rotation_matrix_from_euler(const Eigen::Vector3d &e, praxis::axis_order order)
{
    if(order == praxis::axis_order::zyx)
    {
        return rotate_z(e[0])*rotate_y(e[1])*rotate_x(e[2]);
    }
    else if(order == praxis::axis_order::xyz)
    {
        return rotate_x(e[0])*rotate_y(e[1])*rotate_z(e[2]);
    }

    const std::array<std::uint8_t, 3> axes = praxis::axis_indices(order);
    return rotation_matrix_from_axis_angle(Eigen::Vector3d::Unit(axes[0]), e[0]) * 
           rotation_matrix_from_axis_angle(Eigen::Vector3d::Unit(axes[1]), e[1]) *
           rotation_matrix_from_axis_angle(Eigen::Vector3d::Unit(axes[2]), e[2]);
}

//TASK: 1j -- Implement the special case for calculating the Euler ZYX.
//REFERENCE: Appendix B.1, Euler Angles, pp. 577–580, MR 2019
Eigen::Vector3d euler_from_rotation_matrix(const Eigen::Matrix3d &r, praxis::axis_order order)
{
    if(order == praxis::axis_order::zyx)
    {
        const double alpha = std::atan2(r(1, 0), r(0, 0));

        const double beta = std::atan2(
            -r(2, 0),
            std::sqrt(r(0, 0) * r(0, 0) + r(1, 0) * r(1, 0))
        );

        const double gamma = std::atan2(r(2, 1), r(2, 2));

        return Eigen::Vector3d(alpha, beta, gamma);
    }

    const std::array<std::uint8_t, 3> axes = praxis::axis_indices(order);

    return r.eulerAngles(axes[0], axes[1], axes[2]);
}

//TASK: 2b
//REFERENCE: Definition 3.13, Equation (3.62), page 87, MR 2019
Eigen::Matrix3d rotation_matrix_from_transform(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R = tf.block<3, 3>(0, 0);
    return R;
}

//TASK: 2c
//REFERENCE: Definition 3.13, Equation (3.62), page 87, MR 2019
Eigen::Matrix4d transformation_matrix_from_rotation_position(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    pose.block<3, 3>(0, 0) = r;
    pose.block<3, 1>(0, 3) = p;
    return pose;
}

//TASK: 2d
//REFERENCE: Definition 3.13, Equation (3.62), page 87, MR 2019
Eigen::Matrix4d transformation_matrix_from_position(const Eigen::Vector3d &p)
{
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    pose.block<3, 1>(0, 3) = p;
    return pose;
}

//TASK: 2e
//REFERENCE: Definition 3.13, Equation (3.62), page 87, MR 2019
Eigen::Matrix4d transformation_matrix_from_rotation(const Eigen::Matrix3d &r)
{
    Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();

    pose.block<3, 3>(0, 0) = r;
    return pose;
}

//TASK: 2f
//REFERENCE: Proposition 3.15, Equation (3.64), page 88, MR 2019
Eigen::Matrix4d inverse(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);

    Eigen::Matrix3d rt = r.transpose();

    Eigen::Matrix4d inversepose = Eigen::Matrix4d::Identity();

    inversepose.block<3, 3>(0, 0) = rt;
    inversepose.block<3, 1>(0, 3) = -rt * p;

    return inversepose;
}

}
