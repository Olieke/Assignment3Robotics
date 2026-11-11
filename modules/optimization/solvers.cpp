//Assignment 3
#include "optimization/solvers.h"

#include <Eigen/LU>

#include <cmath>

namespace ais4104::optimization {

//TASK: T.1a
//REFERENCES: Modern Robotics (2019), Sec. 6.2.1–6.2.2; Skrede, Sec. 2.6.2.

praxis::expected<Eigen::VectorXd, praxis::refusal> newton_raphson(const praxis::optimization::problem &p, const Eigen::VectorXd &x0, std::uint32_t max_iterations,
                                                                  std::vector<Eigen::VectorXd> &iterates)
{
    //append x0 to iterates
    iterates.push_back(x0);

    //for each iteration:
    Eigen::VectorXd x = x0;
    for (std::uint32_t i = 0; i < max_iterations; ++i) {

        //    evaluate the error e(x) with p.error
        auto error_result = p.error(x);
        if (!error_result.has_value()) {
            return praxis::unexpected(error_result.error());
        }
        const Eigen::VectorXd &e = error_result.value();

        //    return x when p.converged holds for e(x)
        if (p.converged(e)) {
            return x;
        }

        //    make the step
        auto jacobian_result = p.jacobian(x);
        if (!jacobian_result.has_value()) {
            return praxis::unexpected(jacobian_result.error());
        }
        const Eigen::MatrixXd &J = jacobian_result.value();

        // dx = J^+ * e, med pseudoinversen avhengig av formen på J
        Eigen::VectorXd dx;
        if (J.rows() < J.cols()) {
            dx = J.transpose() * (J * J.transpose()).partialPivLu().solve(e);
        } else if (J.rows() > J.cols()) {
            dx = (J.transpose() * J).partialPivLu().solve(J.transpose() * e);
        } else {
            dx = J.partialPivLu().solve(e);
        }
        x += dx;

        //    append the new x to iterates
        iterates.push_back(x);
    }

    //past max_iterations: return praxis::unexpected(praxis::refusal::no_solution)
    return praxis::unexpected(praxis::refusal::no_solution);
}
//Newton-Raphson test: The function (g(x)=x^3-2x-5) converged after 7 iterations, finding a root at approximately 2.09455.
//The function (g(x)=\arctan(x)) failed with a degenerate result. Although it has a root at (x=0),
//Newton-Raphson failed because the initial guess (x_0=-3) caused large steps away from the root.
//This shows that Newton-Raphson depends on the initial guess and may not always converge.



//TASK: T.1b
//REFERENCE: Skrede, Sec. 2.6.3.
praxis::expected<Eigen::VectorXd, praxis::refusal> gradient_descent(const praxis::optimization::problem &p, const Eigen::VectorXd &x0, double step_length,
                                                                    std::uint32_t max_iterations, std::vector<Eigen::VectorXd> &iterates)
{
    //append x0 to iterates
    iterates.push_back(x0);

    //for each iteration:
    Eigen::VectorXd x = x0;
    for (std::uint32_t i = 0; i < max_iterations; ++i) {

        //    evaluate the gradient of f at x with p.gradient
        auto gradient_result = p.gradient(x);
        if (!gradient_result.has_value()) {
            return praxis::unexpected(gradient_result.error());
        }
        const Eigen::VectorXd &g = gradient_result.value();

        //    return x when p.stationary holds for the gradient
        if (p.stationary(g)) {
            return x;
        }

        //    make the step
        x -= step_length * g;

        //    append the new x to iterates
        iterates.push_back(x);
    }

    //past max_iterations: return praxis::unexpected(praxis::refusal::no_solution)
    return praxis::unexpected(praxis::refusal::no_solution);
}


// Gradient Descent test:
// The function f(x) = x^4 - 3x^2 + x failed with a degenerate result.
// The large gradient at the initial guess x0 = 3.5 caused large steps away from the minimum.
//
// The function f(x) = (x - 1)^2 converged after 132 iterations,
// finding the minimum at x = 1.
//
// This shows that Gradient Descent depends on the step length and the shape of the function.

}\

