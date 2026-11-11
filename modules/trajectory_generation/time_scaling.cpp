//Assignment 3
#include "trajectory_generation/time_scaling.h"

#include <cmath>

namespace ais4104::trajectory {

//TASK: T.4a
//REFERENCE: Modern Robotics (2019), Sec. 9.2.2.1, Eqs. (9.9)–(9.10).
praxis::expected<praxis::trajectory::scaling_sample, praxis::refusal> cubic(double t, double duration)
{
    if (!std::isfinite(duration) || duration <= 0.0) {
        return praxis::unexpected(praxis::refusal::unsupported_input);
    }

    const double T2 = duration * duration;
    const double T3 = T2 * duration;

    const double s = 3.0 * t * t / T2 - 2.0 * t * t * t / T3;
    const double s_dot = 6.0 * t / T2 - 6.0 * t * t / T3;
    const double s_ddot = 6.0 / T2 - 12.0 * t / T3;

    return praxis::trajectory::scaling_sample{s, s_dot, s_ddot};
}

//TASK: T.4b
//REFERENCE: Modern Robotics (2019), Sec. 9.2.2.1, Fifth-Order Polynomials
praxis::expected<praxis::trajectory::scaling_sample, praxis::refusal> quintic(double t, double duration)
{
    if (!std::isfinite(duration) || duration <= 0.0) {
        return praxis::unexpected(praxis::refusal::unsupported_input);
    }

    const double u = t / duration;
    const double T2 = duration * duration;

    const double s = 10.0 * std::pow(u, 3)
                   - 15.0 * std::pow(u, 4)
                   + 6.0 * std::pow(u, 5);

    const double s_dot = (30.0 * std::pow(u, 2)
                        - 60.0 * std::pow(u, 3)
                        + 30.0 * std::pow(u, 4)) / duration;

    const double s_ddot = (60.0 * u
                         - 180.0 * std::pow(u, 2)
                         + 120.0 * std::pow(u, 3)) / T2;

    return praxis::trajectory::scaling_sample{s, s_dot, s_ddot};
}

}
