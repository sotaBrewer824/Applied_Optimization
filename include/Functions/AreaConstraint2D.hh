#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Eigen/Dense>
#include <cmath>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class AreaConstraint2D : public FunctionBaseSparse {
    public:
        // Area constraint: 1/2*det(v_01 | V02) >= eps
        // f = -1/2*((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0)) + eps <= 0
        // constructor
        AreaConstraint2D(const int _n, const int _idx0, const int _idx1, const int _idx2, const double _eps = 1e-10)
                : FunctionBaseSparse(), n_(_n), idx0_(_idx0), idx1_(_idx1), idx2_(_idx2), eps_(_eps) {}

        // number of unknowns
        inline virtual int n_unknowns() override { return n_; }

        // function evaluation
        // _x stores the coordinates of all nodes
        inline virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            // Extract coordinates
            double x0 = _x[2*idx0_];     // x-coordinate of node 0
            double y0 = _x[2*idx0_+1]; // y-coordinate of node 0
            double x1 = _x[2*idx1_];
            double y1 = _x[2*idx1_+1];
            double x2 = _x[2*idx2_];
            double y2 = _x[2*idx2_+1];

            // Compute the area constraint function value
            double f = -0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) + eps_;
            return f;
            //------------------------------------------------------//
        }

        // gradient evaluation
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.setZero();
            //------------------------------------------------------//
            // Extract coordinates
            double x0 = _x[2*idx0_];
            double y0 = _x[2*idx0_+1];
            double x1 = _x[2*idx1_];
            double y1 = _x[2*idx1_+1];
            double x2 = _x[2*idx2_];
            double y2 = _x[2*idx2_+1];

            // Partial derivatives
            _g[2*idx0_]     = 0.5 * ((y2 - y1)); // x0
            _g[2*idx0_+1] = 0.5 * ((x1 - x2)); // y0
            _g[2*idx1_]     = 0.5 * ((y0 - y2)); // x1
            _g[2*idx1_+1] = 0.5 * ((x2 - x0)); // y1
            _g[2*idx2_]     = 0.5 * ((y1 - y0)); // x2 
            _g[2*idx2_+1] = 0.5 * ((x0 - x1)); // y2
            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        inline virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            _h.setZero();
            //------------------------------------------------------//
            // The Hessian matrix is zero because the second derivatives are constant
            _h.coeffRef(2*idx0_, 2*idx1_+1) = -0.5;
            _h.coeffRef(2*idx0_, 2*idx2_+1) = 0.5;
            _h.coeffRef(2*idx0_+1, 2*idx1_) = 0.5;
            _h.coeffRef(2*idx0_+1, 2*idx2_) = -0.5;
            _h.coeffRef(2*idx1_, 2*idx0_+1) = 0.5;
            _h.coeffRef(2*idx1_, 2*idx2_+1) = -0.5;
            _h.coeffRef(2*idx1_+1, 2*idx0_) = -0.5;
            _h.coeffRef(2*idx1_+1, 2*idx2_) = 0.5;
            _h.coeffRef(2*idx2_, 2*idx0_+1) = -0.5;
            _h.coeffRef(2*idx2_, 2*idx1_+1) = 0.5;
            _h.coeffRef(2*idx2_+1, 2*idx0_) = 0.5;
            _h.coeffRef(2*idx2_+1, 2*idx1_) = -0.5;
            //------------------------------------------------------//
        }

    private:
        int n_;
        // index of the nodes
        int idx0_;
        int idx1_;
        int idx2_;

        double eps_;
    };

//=============================================================================
}
