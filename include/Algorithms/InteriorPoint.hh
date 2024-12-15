#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Functions/InteriorPointProblem.hh>
#include <Algorithms/NewtonMethods.hh>
#include <Utils/OptimizationStatistic.hh>
#include <iostream>
#include <vector>
#include <memory>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class InteriorPoint {
    public:
        // LA typedefs
        using Vec = FunctionBaseSparse::Vec;

        static Vec solve(FunctionBaseSparse *_obj, const Vec& _initial_x, const std::vector<FunctionBaseSparse*>& _constraints,
                const double _eps = 1e-4, const double _mu = 10.0, const int _max_iters = 1000) {
            std::cerr << "******** Interior Point ********" << std::endl;

            // Construct log-barrier problem
            InteriorPointProblem problem(_obj, _constraints);
            auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(&problem);

            // Count number of iterations
            int iter(0);

            // Current barrier parameter
            double t = 1.0;

            // Number of constraints
            double m = _constraints.size();

            // Points
            Vec x = _initial_x;

            // Setup
            bool converged = false;
            
            while (iter < _max_iters) {
                problem.t() = t; // Update barrier parameter

                // Centering step using Newton's method with projected Hessian
                x = NewtonMethods::solve_with_projected_hessian(opt_st.get(), converged, x, 10., _eps, 100000);

                // Stopping criterion: check duality gap
                if (m / t < _eps) {
                    break;
                }

                // Increase barrier parameter
                t *= _mu;
                ++iter;
            }

            return x;
        }
    };
    //=============================================================================

}
