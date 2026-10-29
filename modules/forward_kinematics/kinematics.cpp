//Assignment 2
#include "forward_kinematics/kinematics.h"

#include <span>
#include <vector>

namespace ais4104::manipulator {

//TASK: 2a
//REFERENCE: Modern Robotics, Chapter 4.1.1, Eq. (4.14).
praxis::expected<Eigen::Matrix4d, praxis::refusal> forward_kinematics(
    const praxis::rigid_motion::screw_ops &screw,
    const Eigen::Matrix4d &m,
    std::span<const Eigen::Vector6d> space_screws,
    const Eigen::VectorXd &theta)
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    const std::size_t n = space_screws.size();

    for (std::size_t i = 0; i < n; ++i)
    {
        // Convert the index to Eigen's index type.
        const Eigen::Matrix4d exp_i =
            screw.matrix_exponential_screw(
                space_screws[i],
                theta(static_cast<Eigen::Index>(i)));

        T = T * exp_i;
    }

    T = T * m;

    return T;
}


//TASK: 2b
//REFERENCE: Modern Robotics, Chapter 4.1.3, Eq. (4.16), p. 147.
// Convert space-frame screw axes to body-frame screw axes.
praxis::expected<std::vector<Eigen::Vector6d>, praxis::refusal> body_screws_from_space(
    const praxis::rigid_motion::screw_ops &screw,
    const praxis::rigid_motion::frame_ops &frames,
    const Eigen::Matrix4d &m,
    std::span<const Eigen::Vector6d> space_screws)
{
    // B_i = [Ad_{M^-1}] S_i
    const Eigen::Matrix4d m_inv = m.inverse();

    std::vector<Eigen::Vector6d> body_screws;
    body_screws.reserve(space_screws.size());

    for (const Eigen::Vector6d &space_screw : space_screws)
    {
        const auto body_screw = screw.adjoint_map(space_screw, m_inv);
        if (!body_screw)
            return praxis::unexpected(body_screw.error());

        body_screws.push_back(body_screw.value());
    }

    return body_screws;
}

//TASK: 2c
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 4.1.3, eq. (4.16), p. 147. body-frame PoE formula.
praxis::expected<Eigen::Matrix4d, praxis::refusal> body_forward_kinematics(const praxis::rigid_motion::screw_ops &screw, const praxis::rigid_motion::frame_ops &frames,
                                                                           const Eigen::Matrix4d &m, std::span<const Eigen::Vector6d> space_screws, const Eigen::VectorXd &theta)
{
    auto result = body_screws_from_space(
        screw, frames, m, space_screws);

    const auto body_screws = result.value();
    Eigen::Matrix4d transform = m;

    for (std::size_t joint = 0; joint < body_screws.size(); ++joint)
    {
        const auto angle = theta(
            static_cast<Eigen::Index>(joint));

        const Eigen::Matrix4d joint_transform =
            screw.matrix_exponential_screw(
                body_screws[joint], angle);

        transform *= joint_transform;
    }

    return transform;
}

}
