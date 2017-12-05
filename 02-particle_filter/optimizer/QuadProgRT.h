/** \file
 * \copyright
 * TaCo - Task-Space Control Library.
 * Copyright (2016) Brian Soe <bsoe@stanford.edu>,
 * Copyright (2010) Gael Guennebaud,
 * Copyright (2008) Angelo Furfaro,
 * Copyright (2006) Luca Di Gaspero
 *
 * This is a modified of eigenquadprog, to work with TaCo.
 * eigenquadprog is itself a port made by Gael Guennebaud of uQuadProg++, to work with Eigen data structures.
 * uQuadProg++ is itself a port made by Angelo Furfaro of QuadProg++, to work with ublas data structures.
 * QuadProg++ is originally developed by Luca Di Gaspero
 *
 * -----
 *
 * This file is a porting of QuadProg++ routine, originally developed
 * by Luca Di Gaspero, exploiting uBlas data structures for vectors and
 * matrices instead of native C++ array.
 *
 * uquadprog is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * uquadprog is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with uquadprog; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * -----
 *
 * Author: Luca Di Gaspero
 * DIEGM - University of Udine, Italy
 * l.digaspero@uniud.it
 * http://www.diegm.uniud.it/digaspero/
 *
 * LICENSE
 *
 * This file is part of QuadProg++: a C++ library implementing
 * the algorithm of Goldfarb and Idnani for the solution of a (convex)
 * Quadratic Programming problem by means of an active-set dual method.
 * Copyright (C) 2007-2009 Luca Di Gaspero.
 * Copyright (C) 2009 Eric Moyer.
 *
 * QuadProg++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QuadProg++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with QuadProg++. If not, see <http://www.gnu.org/licenses/>.
<br>
Author: Brian Soe <bsoe@stanford.edu>
 */




#ifndef TACO_OPTIMIZER_QUADPROGRT_H_
#define TACO_OPTIMIZER_QUADPROGRT_H_

#include "QP.h"

#include <Eigen/Dense>

namespace taco {

/** \brief Quadratic programming which implements the Goldfarb-Idnani active-set dual method.
 *
 * Limited to solution of stricty convex quadratic programs.
 *
 * This is a modified of eigenquadprog, to work with TaCo.
 * eigenquadprog is itself a port made by Gael Guennebaud of uQuadProg++, to work with Eigen data structures.
 * uQuadProg++ is itself a port made by Angelo Furfaro of QuadProg++, to work with ublas data structures.
 * QuadProg++ is originally developed by Luca Di Gaspero.
 *
 * GNU General Public License
 *
 * \note The author will be grateful if the researchers using this software will
 * acknowledge the contribution of this function in their research papers.
 * References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
 *             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.
 *
 * \ingroup optimizer
 */
class QuadProgRT : public QP
{
public:

    QuadProgRT();
    virtual ~QuadProgRT();


    virtual void setObjective(const Eigen::MatrixXd& A,
                              const Eigen::VectorXd& B,
                              bool condition = false);

    virtual void setEqualities(const Eigen::MatrixXd& E,
                               const Eigen::VectorXd& E0,
                               bool condition = false);

    virtual void setInequalities(const Eigen::MatrixXd& I,
                                 const Eigen::VectorXd& I0,
                                 bool condition = false);

    virtual void setConditionNumber(double epsilon);

    virtual double solve(Eigen::VectorXd& x);

    /** \brief Not implemented. Function calls solve(x) instead. */
    virtual double solve(Eigen::VectorXd& x, const Eigen::VectorXd x_initial);

    /** \brief Check if matrices are compatible with QuadProg.
     * Throws an error if dimensions are inconsistent. */
    void checkProblem();

private:

    void resizeEI();

    void printWarning(const std::string& message);

    // Original problem
    Eigen::MatrixXd obj_A_;
    Eigen::VectorXd obj_B_;

    double epsilon_ = 0.0001;

    // Problem variables
    Eigen::MatrixXd G;
    Eigen::VectorXd g0;
    Eigen::MatrixXd CE;
    Eigen::VectorXd ce0;
    Eigen::MatrixXd CI;
    Eigen::VectorXd ci0;
    //Eigen::VectorXd x_;
    Eigen::VectorXd x0_; // if none provided

    // Infinity
    double inf_;

    // Preallocate
    Eigen::MatrixXd R_, J_;
    Eigen::VectorXd s_, z_, r_, d_, np_, u_, x_old_, u_old_;
    Eigen::VectorXi A_, A_old_, iai_;
    std::vector<bool> iaexcl_;

    // Utility functions for updating some data needed by the solution method
    void compute_d(Eigen::VectorXd& d, const Eigen::MatrixXd& J, const Eigen::VectorXd& np);
    void update_z(Eigen::VectorXd& z, const Eigen::MatrixXd& J, const Eigen::VectorXd& d, int iq);
    void update_r(const Eigen::MatrixXd& R, Eigen::VectorXd& r, const Eigen::VectorXd& d, int iq);
    bool add_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXd& d, int& iq, double& rnorm);
    void delete_constraint(Eigen::MatrixXd& R, Eigen::MatrixXd& J, Eigen::VectorXi& A, Eigen::VectorXd& u, int p, int& iq, int l);

    // Utility functions for computing the euclidean distance between two numbers
    double distance(double a, double b);

};

}

#include "QuadProgRT.cpp"

#endif // #define TACO_OPTIMIZER_QUADPROGRT_H_
