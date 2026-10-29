//Assignment 2
#include "velocity_kinematics/jacobian.h"

#include <span>
#include <vector>
#include <cstddef>

namespace ais4104::manipulator {

//TASK: 4a
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 5.1.1, eq. (5.11).
// Modern Robotics, Section 3.3.2, Definition 3.20, p. 98 (adjoint representation).
// Modern Robotics, Section 3.3.3, Proposition 3.25, Eq. (3.88), p. 104 (matrix exponential).
praxis::expected<praxis::manipulator::jacobian, praxis::refusal> space_jacobian(const praxis::rigid_motion::screw_ops &screw, std::span<const Eigen::Vector6d> space_screws,
                                                                                const Eigen::VectorXd &theta)
{
    const std::size_t n = space_screws.size();

    praxis::manipulator::jacobian J(
        6, static_cast<Eigen::Index>(n));

    if (n == 0)
        return J;

    // The first column is the first screw axis.
    J.col(0) = space_screws[0];

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    for (std::size_t i = 1; i < n; ++i)
    {
        const Eigen::Matrix4d exp_transform =
            screw.matrix_exponential_screw(
                space_screws[i - 1],
                theta(static_cast<Eigen::Index>(i - 1)));

        T = T * exp_transform;

        const auto column =
            screw.adjoint_map(space_screws[i], T);

        if (!column)
            return praxis::unexpected(column.error());

        J.col(static_cast<Eigen::Index>(i)) = column.value();
    }

    return J;
}

//TASK: 4b
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 5.1.2, eq. (5.18).
// Modern Robotics, Ch. 4.1.1, Eq. (4.10), p. 138:
// Space-frame product-of-exponentials formula for forward kinematics.
praxis::expected<praxis::manipulator::jacobian, praxis::refusal> body_jacobian(const praxis::rigid_motion::screw_ops &screw, const praxis::rigid_motion::frame_ops &frames,
                                                                               const praxis::manipulator::forward_kinematics_ops &forward, const Eigen::Matrix4d &m,
                                                                               std::span<const Eigen::Vector6d> space_screws, const Eigen::VectorXd &theta)
{
    (void)frames;
    (void)forward;

    // Number of joints
    const std::size_t n = space_screws.size();

    // Calculate the space Jacobian
    const auto Js = space_jacobian(screw, space_screws, theta);

    if (!Js)
        return praxis::unexpected(Js.error());

    // Calculate T(theta)
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    for (std::size_t i = 0; i < n; ++i)
    {
        T = T * screw.matrix_exponential_screw(
            space_screws[i],
            theta(static_cast<Eigen::Index>(i)));
    }

    T = T * m;

    // Inverse transformation
    const Eigen::Matrix4d T_inv = T.inverse();

    // Construct the body Jacobian
    praxis::manipulator::jacobian Jb(
        6, static_cast<Eigen::Index>(n));

    for (std::size_t i = 0; i < n; ++i)
    {
        const Eigen::Vector6d space_column =
            Js.value().col(static_cast<Eigen::Index>(i));

        const auto body_column =
            screw.adjoint_map(space_column, T_inv);

        if (!body_column)
            return praxis::unexpected(body_column.error());

        Jb.col(static_cast<Eigen::Index>(i)) = body_column.value();
    }

    return Jb;
}
} // namespace ais4104::manipulator


