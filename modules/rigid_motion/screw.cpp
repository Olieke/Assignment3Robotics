#include "rigid_motion/frame.h"
#include "rigid_motion/screw.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>

namespace ais4104::rigid_motion {

//TASK: 1b
//REFERENCE:Definition 3.7, Equation (3.30), page 75, MR 2019
Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v)
{
    Eigen::Matrix3d ω;
    ω<< 0 , -v(2), v(1),
    v(2)  , 0 , -v(0) ,
    -v(1), v(0) , 0;
    return ω;
}

//TASK: 1c
//REFERENCE:Definition 3.7, Equation (3.30), page 75, MR 2019
Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d &m)
{
    Eigen::Vector3d ω;
    ω << m(2, 1), m(0, 2), m(1,0);
    return ω;
}

//TASK: 3b
//REFERENCE:Definition 3.24, page 102, MR 2019
praxis::expected<Eigen::Vector6d, praxis::refusal> screw_axis_from_point_direction_pitch(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::Vector3d sq = s.cross(q);
    Eigen::Vector3d v = -sq + h * s;
    Eigen::Vector6d S;
    S.head<3>()=s;
    S.tail<3>()=v;
    return S;
}

//TASK: 3c
//REFERENCE:Definition 3.24, page 102, MR 2019
Eigen::Vector6d screw_axis_from_angular_linear(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::Vector6d S;

    if (w.norm() > 1e-10)
    {
        double length = w.norm();
        S.head<3>() = w / length;
        S.tail<3>() = v / length;
    }
    else
    {
        S.head<3>() = Eigen::Vector3d::Zero();
        S.tail<3>() = v.normalized();
    }

    return S;
}

//TASK: 3d
//REFERENCE: Definition 3.24, page 102, MR 2019
Eigen::Vector6d twist_from_angular_linear(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::Vector6d S;
    S.head<3>()=w;
    S.tail<3>()=v;
    return S;
}

//TASK: 3e
//REFERENCE: Section 3.3.2.2, page 102, case (a), MR 2019
praxis::expected<Eigen::Vector6d, praxis::refusal> twist_from_screw(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h, double angular_velocity)
{
    Eigen::Vector3d sq = s.cross(q);
    Eigen::Vector3d v = -sq + h * s;
    Eigen::Vector6d S;
    S.head<3>()=s;
    S.tail<3>()=v;
    Eigen::Vector6d V = S * angular_velocity;
    return V;
}

//TASK: 3f
//REFERENCE:Definition 3.20, page 98, MR 2019
praxis::expected<Eigen::Matrix6d, praxis::refusal> adjoint_matrix_from_rotation_position(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix6d Adt = Eigen::Matrix6d::Zero();

    Adt.block<3, 3>(0, 0) = r;
    Adt.block<3, 3>(3, 3) = r;

    Eigen::Matrix3d p_skew = skew_symmetric(p);

    Adt.block<3, 3>(3, 0) = p_skew * r;

    return Adt;
}

//TASK: 3g
//REFERENCE:Definition 3.20, page 98, MR 2019
praxis::expected<Eigen::Matrix6d, praxis::refusal> adjoint_matrix_from_transform(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);

    Eigen::Matrix6d Adt = Eigen::Matrix6d::Zero();

    Adt.block<3, 3>(0, 0) = r;
    Adt.block<3, 3>(3, 3) = r;

    Eigen::Matrix3d p_skew = skew_symmetric(p);

    Adt.block<3, 3>(3, 0) = p_skew * r;

    return Adt;
}

//TASK: 3h
//REFERENCE: Definition 3.20, page 98, MR 2019
praxis::expected<Eigen::Vector6d, praxis::refusal> adjoint_map(const Eigen::Vector6d &t, const Eigen::Matrix4d &tf) {
    auto Adt = adjoint_matrix_from_transform(tf);

    if (!Adt)
        return praxis::unexpected(Adt.error());

    return Adt.value() * t;
}

//TASK: 3i
//REFERENCE:Equation (3.85), page 102, MR 2019
Eigen::Matrix4d twist_matrix_from_angular_linear(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::Matrix4d V = Eigen::Matrix4d::Zero();
    V.block<3, 3>(0, 0) = skew_symmetric(w);
    V.block<3, 1>(0, 3) = v;
    return V;
}

//TASK: 3j
//REFERENCE: Equation (3.85), page 102, MR 2019
Eigen::Matrix4d twist_matrix_from_twist(const Eigen::Vector6d &t)
{
    Eigen::Vector3d w = t.head<3>();
    Eigen::Vector3d v = t.tail<3>();

    Eigen::Matrix4d V = Eigen::Matrix4d::Zero();

    V.block<3, 3>(0, 0) = skew_symmetric(w);
    V.block<3, 1>(0, 3) = v;

    return V;
}

//TASK: 3k
//REFERENCE: Proposition 3.11, Equation (3.51), page 82, MR 2019
Eigen::Matrix3d matrix_exponential_so3(const Eigen::Vector3d &w, double theta_radians)
{
    Eigen::Matrix3d R = rotation_matrix_from_axis_angle(w, theta_radians);
    return R;
}

//TASK: 3l
//REFERENCE: Equations (3.87) and (3.88), page 103, MR 2019
Eigen::Matrix4d matrix_exponential_se3(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta_radians)
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d W = skew_symmetric(w);

    T.block<3, 3>(0, 0) = matrix_exponential_so3(w, theta_radians);

    T.block<3, 1>(0, 3) =
        (Eigen::Matrix3d::Identity() * theta_radians
        + (1 - cos(theta_radians)) * W
        + (theta_radians - sin(theta_radians)) * W * W) * v;

    return T;
}

//TASK: 3m
//REFERENCE: Proposition 3.25, Equation (3.88), page 103, MR 2019
Eigen::Matrix4d matrix_exponential_screw(const Eigen::Vector6d &s, double theta_radians)
{
    Eigen::Vector3d w = s.head<3>();
    Eigen::Vector3d v = s.tail<3>();
    return matrix_exponential_se3(w, v, theta_radians);
}

//TASK: 3n
//REFERENCE:Equations (3.58)–(3.61), pages 85–86, MR 2019
praxis::expected<std::pair<Eigen::Vector3d, double>, praxis::refusal> matrix_logarithm_so3(const Eigen::Matrix3d &r)
{
    if (r.isApprox(Eigen::Matrix3d::Identity()))
    {
        return std::make_pair(Eigen::Vector3d::UnitZ(), 0.0);
    }

    if (std::abs(r.trace() + 1.0) < 1e-8)
    {
        Eigen::Vector3d w;

        if (1.0 + r(2, 2) > 1e-8)
        {
            w << r(0, 2), r(1, 2), 1.0 + r(2, 2);
        }
        else if (1.0 + r(1, 1) > 1e-8)
        {
            w << r(0, 1), 1.0 + r(1, 1), r(2, 1);
        }
        else
        {
            w << 1.0 + r(0, 0), r(1, 0), r(2, 0);
        }

        return std::make_pair(w.normalized(), std::acos(-1.0));
    }

    double theta = std::acos((r.trace() - 1.0) / 2.0);

    Eigen::Matrix3d W =
        (r - r.transpose()) / (2.0 * std::sin(theta));

    Eigen::Vector3d w;
    w << W(2, 1), W(0, 2), W(1, 0);

    return std::make_pair(w, theta);
}

//TASK: 3o
//REFERENCE: Section 3.3.3.2, Eq. (3.92), MR 2019
praxis::expected<std::pair<Eigen::Vector6d, double>, praxis::refusal> matrix_logarithm_se3_rp(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{

    auto result = matrix_logarithm_so3(r);

    if (!result)
    {
        return praxis::unexpected(result.error());
    }

    Eigen::Vector3d w = result.value().first;
    double theta = result.value().second;

    Eigen::Matrix3d W = skew_symmetric(w);

    Eigen::Matrix3d G_inv =
        Eigen::Matrix3d::Identity() / theta
        - 0.5 * W
        + (1.0 / theta
           - 0.5 / std::tan(theta / 2.0)) * W * W;

    Eigen::Vector3d v = G_inv * p;

    Eigen::Vector6d S;
    S.head<3>() = w;
    S.tail<3>() = v;

    return std::make_pair(S, theta);
}

//TASK: 3p
//REFERENCE: Equation (3.90), Section 3.3.3.2, MR 2019
praxis::expected<std::pair<Eigen::Vector6d, double>, praxis::refusal> matrix_logarithm_se3(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);
    return matrix_logarithm_se3_rp(r, p);
}

}
