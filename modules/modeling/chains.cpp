//Assignment 2
#include "modeling/chains.h"

#include <praxis/rigid_motion/angles.h>

namespace ais4104::modeling {


//TASK: 1a
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 4.1.1.
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 3.3.2 (Screw Axes).
//REFERENCE: Lynch & Park, Modern Robotics, Sec. 3.3.1, Definition 3.13, Eq. (3.62) (Homogeneous Transformation Matrices).
std::pair<Eigen::Matrix4d, std::vector<Eigen::Vector6d>>
hardcoded_kr6r900sixx_chain(
    const praxis::rigid_motion::screw_ops &screw)
{
    Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
    M.block<3, 3>(0, 0) << 0, 0, 1,
                       0, 1, 0,
                      -1, 0, 0;

    M(0, 3) = 0.980;
    M(1, 3) = 0.000;
    M(2, 3) = 0.435;

    std::vector<Eigen::Vector3d> q = {
        {0.000, 0.000, 0.000},
        {0.025, 0.000, 0.400},
        {0.480, 0.000, 0.400},
        {0.515, 0.000, 0.435},
        {0.900, 0.000, 0.435},
        {0.980, 0.000, 0.435}
    };

    std::vector<Eigen::Vector3d> w = {
        { 0,  0, -1},
        { 0,  1,  0},
        { 0,  1,  0},
        {-1,  0,  0},
        { 0,  1,  0},
        {-1,  0,  0}
    };

    std::vector<Eigen::Vector6d> S;
    S.reserve(6);

    for (int i = 0; i < 6; ++i)
    {
        auto axis =
            screw.screw_axis_from_point_direction_pitch(
                q[i], w[i], 0.0);

        if (!axis)
        {
            return {Eigen::Matrix4d::Zero(), {}};
        }

        S.push_back(axis.value());
    }

    return {M, S};
}


//TASK: 1b
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 4.1.1.
//REFERENCE: Lynch & Park, Modern Robotics, Ch. 3.3.2 (Screw Axes).
//REFERENCE: Lynch & Park, Modern Robotics, Sec. 3.3.1, Definition 3.13, Eq. (3.62) (Homogeneous Transformation Matrices).
std::pair<Eigen::Matrix4d, std::vector<Eigen::Vector6d>>
hardcoded_lbr_med14_r820_chain(
    const praxis::rigid_motion::screw_ops &screw)
{
    Eigen::Matrix4d M = Eigen::Matrix4d::Identity();

    M(0, 3) = 0.0;
    M(1, 3) = 0.0;
    M(2, 3) = 1.3055;

    std::vector<Eigen::Vector3d> q = {
        {0.000, 0.000, 0.000},
        {0.000, 0.000, 0.3595},
        {0.000, 0.000, 0.000},
        {0.000, 0.000, 0.7795},
        {0.000, 0.000, 0.000},
        {0.000, 0.000, 1.1795},
        {0.000, 0.000, 1.306}
    };

    std::vector<Eigen::Vector3d> w = {
        {0, 0, 1},
        {0, 1, 0},
        {0, 0, 1},
        {0, -1, 0},
        {0, 0, 1},
        {0, 1, 0},
        {0, 0, 1}
    };

    std::vector<Eigen::Vector6d> S;
    S.reserve(7);

    for (int i = 0; i < 7; ++i)
    {
        auto axis =
            screw.screw_axis_from_point_direction_pitch(
                q[i], w[i], 0.0);

        if (!axis)
            return {Eigen::Matrix4d::Zero(), {}};

        S.push_back(axis.value());
    }

    return {M, S};
}



//TASK: 3a
//REFERENCE: Modern Robotics, Chapter 3.3.1, Definition 3.13, p. 87.
Eigen::Matrix4d hardcoded_kr6r900sixx_curved_rod_tool_offset(const praxis::rigid_motion::frame_ops &frames)
{
    Eigen::Matrix4d Tft = Eigen::Matrix4d::Identity();

    Tft.block<3, 3>(0, 0) << 0, 0, 1,
                          0, 1, 0,
                         -1, 0, 0;


    Tft(0, 3) =  0.082;
    Tft(1, 3) =  0.000;
    Tft(2, 3) =  0.132;

    return Tft;
}

//TASK: 3b
//REFERENCE: Section 3.2 page 72, MR pre-print 2019          (Rx, Ry, Rz)
//REFERENCE: Equation (B.1) page 578, MR pre-print 2019      (Rz*Ry*Rx)
Eigen::Matrix4d hardcoded_kr6r900sixx_humerus_tool_offset(const praxis::rigid_motion::frame_ops &frames)
{
    const double deg = std::numbers::pi / 180.0;

    const double alpha = -180.0 * deg;
    const double beta = 75.0 * deg;
    const double gamma = 90.0 * deg;

    Eigen::Matrix3d Rz;
    Rz << std::cos(alpha), -std::sin(alpha), 0,
          std::sin(alpha),  std::cos(alpha), 0,
          0,                0,               1;

    Eigen::Matrix3d Ry;
    Ry << std::cos(beta),  0, std::sin(beta),
          0,               1, 0,
          -std::sin(beta), 0, std::cos(beta);

    Eigen::Matrix3d Rx;
    Rx << 1, 0,                0,
          0, std::cos(gamma), -std::sin(gamma),
          0, std::sin(gamma),  std::cos(gamma);

    Eigen::Matrix3d Rwf = Rz * Ry * Rx;

    Eigen::Matrix4d Tft = Eigen::Matrix4d::Identity();
    Tft.block<3, 3>(0, 0) = Rwf.transpose();

    Tft(0, 3) = -0.1519;
    Tft(1, 3) = -0.0193;
    Tft(2, 3) =  0.0421;

    return Tft;
}

}
