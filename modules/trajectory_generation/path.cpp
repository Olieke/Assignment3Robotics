//Assignment 3
#include "trajectory_generation/path.h"

#include <utility>

namespace ais4104::trajectory {

//TASK: T.3a
//REFERENCE: // Reference: Modern Robotics (2019), Sec. 9.2.1, Eq. (9.3)
praxis::expected<Eigen::VectorXd, praxis::refusal> joint_straight_line(const Eigen::VectorXd &start, const Eigen::VectorXd &end, double s)
{
    if (start.size() != end.size()) {
        return praxis::unexpected(praxis::refusal::unsupported_input);
    }

    return start + s * (end - start);
}

//TASK: T.3b
//REFERENCE: // Reference: Modern Robotics (2019), Sec. 9.2.1, Eq. (9.6)
praxis::expected<Eigen::Matrix4d, praxis::refusal> screw(const praxis::rigid_motion::screw_ops &motion, const praxis::rigid_motion::frame_ops &frames, const Eigen::Matrix4d &start,
                                                         const Eigen::Matrix4d &end, double s)
{
    auto relative_pose = frames.inverse(start) * end;

    auto logarithm = motion.matrix_logarithm_se3(relative_pose);
    if (!logarithm.has_value()) {
        return praxis::unexpected(logarithm.error());
    }

    auto [axis, angle] = logarithm.value();

    auto motion_step = motion.matrix_exponential_screw(axis, s * angle);

    return start * motion_step;
}

//TASK: T.3c
//REFERENCE Modern Robotics (2019), Sec. 9.2.1, Eqs. (9.3), (9.8).
praxis::expected<Eigen::Matrix4d, praxis::refusal> decoupled(const praxis::rigid_motion::screw_ops &motion, const praxis::rigid_motion::frame_ops &frames, const Eigen::Matrix4d &start,
                                                             const Eigen::Matrix4d &end, double s)
{
    auto start_rotation = frames.rotation_matrix_from_transform(start);
    auto end_rotation = frames.rotation_matrix_from_transform(end);

    auto relative_rotation = start_rotation.transpose() * end_rotation;

    auto logarithm = motion.matrix_logarithm_so3(relative_rotation);
    if (!logarithm.has_value()) {
        return praxis::unexpected(logarithm.error());
    }

    auto [axis, angle] = logarithm.value();

    auto rotation_step = motion.matrix_exponential_so3(axis, s * angle);
    auto rotation = start_rotation * rotation_step;

    Eigen::Vector3d start_position = start.block<3, 1>(0, 3);
    Eigen::Vector3d end_position = end.block<3, 1>(0, 3);

    Eigen::Vector3d position = start_position + s * (end_position - start_position);

    return frames.transformation_matrix_from_rotation_position(rotation, position);
}

}
