//Assignment 3
#include "inverse_kinematics/end_effector_frames.h"

namespace ais4104::manipulator {

//TASK: T.2a
praxis::expected<Eigen::Matrix4d, praxis::refusal> tool_pose_from_flange_pose(const Eigen::Matrix4d &flange_pose, const Eigen::Matrix4d &tool_offset)

    {
        return flange_pose * tool_offset;
    }

//TASK: T.2b
praxis::expected<Eigen::Matrix4d, praxis::refusal> flange_pose_from_tool_pose(const praxis::rigid_motion::frame_ops &frames, const Eigen::Matrix4d &tool_pose,
                                                                              const Eigen::Matrix4d &tool_offset)
{
    return tool_pose * frames.inverse(tool_offset);
}

//TASK: T.2c
praxis::expected<Eigen::Vector3d, praxis::refusal> position_from_pose(const Eigen::Matrix4d &pose)
{
    return pose.block<3, 1>(0, 3);
}

//TASK: T.2c
// Reference: Modern Robotics (2019), Sec. 3.3, Eq. (3.62)–(3.63)
praxis::expected<Eigen::Matrix3d, praxis::refusal> orientation_from_pose(const Eigen::Matrix4d &pose)
{
    return pose.block<3, 3>(0, 0);
}

//TASK: T.2d
// Reference: Modern Robotics (2019), Sec. 6.2.1–6.2.2.
praxis::expected<Eigen::VectorXd, praxis::refusal> ik_solve_flange_pose(const praxis::manipulator::kinematics &solver, const Eigen::Matrix4d &flange_pose, const Eigen::VectorXd &j0)
{
    // Check that the initial guess has one value per joint
    if (j0.size() != solver.joint_count()) {
        return praxis::unexpected(praxis::refusal::unsupported_input);
    }

    // Solve inverse kinematics using the default parameters
    return solver.ik_solve(
        flange_pose,
        j0,
        praxis::manipulator::solver_parameters()
    );
}

//TASK: T.2e
// Reference: Modern Robotics (2019), Sec. 3.3 and 6.2.1–6.2.2.
praxis::expected<Eigen::VectorXd, praxis::refusal> ik_solve_pose(const praxis::manipulator::kinematics &solver, const praxis::rigid_motion::frame_ops &frames,
                                                                 const Eigen::Matrix4d &tool_pose, const Eigen::VectorXd &j0, const Eigen::Matrix4d &tool_offset)
{
    // Convert the desired tool pose to a flange pose
    auto flange_pose = flange_pose_from_tool_pose(frames, tool_pose, tool_offset);
    if (!flange_pose.has_value()) {
        return praxis::unexpected(flange_pose.error());
    }

    // Solve inverse kinematics for the flange pose
    return ik_solve_flange_pose(solver, flange_pose.value(), j0);
}

}
