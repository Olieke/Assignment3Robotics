//Assignment 3
#include "trajectory_generation/via_points.h"

namespace ais4104::trajectory {

//TASK: T.5a
//REFERENCE: Modern Robotics (2019), Sec. 9.2.2.1, Third-Order Polynomials, Eq. (9.9)–(9.10)
praxis::expected<cubic_segment, praxis::refusal> segment_coefficients(const praxis::trajectory::via_point &from, const praxis::trajectory::via_point &to)

    {
        const double duration = to.time - from.time;
        if (!std::isfinite(duration) || duration <= 0.0) {
            return praxis::unexpected(praxis::refusal::unsupported_input);
        }

        const Eigen::VectorXd delta = to.position - from.position;

        const Eigen::VectorXd a0 = from.position;
        const Eigen::VectorXd a1 = from.velocity;
        const Eigen::VectorXd a2 = 3.0 * delta / (duration * duration)
                                   - (2.0 * from.velocity + to.velocity) / duration;
        const Eigen::VectorXd a3 = -2.0 * delta / (duration * duration * duration)
                                   + (from.velocity + to.velocity) / (duration * duration);

        return cubic_segment{a0, a1, a2, a3};
    }


//TASK: T.5b
//REFERENCE: Modern Robotics (2019), Sec. 9.3, Eqs. (9.26)–(9.29)
praxis::expected<praxis::trajectory::trajectory_sample, praxis::refusal> segment_sample(const cubic_segment &segment, double dt)
{
    const Eigen::VectorXd position = segment.a0 + segment.a1 * dt
                                      + segment.a2 * dt * dt
                                      + segment.a3 * dt * dt * dt;

    const Eigen::VectorXd velocity = segment.a1 + 2.0 * segment.a2 * dt
                                     + 3.0 * segment.a3 * dt * dt;

    const Eigen::VectorXd acceleration = 2.0 * segment.a2
                                         + 6.0 * segment.a3 * dt;

    return praxis::trajectory::trajectory_sample{position, velocity, acceleration};
}

}
