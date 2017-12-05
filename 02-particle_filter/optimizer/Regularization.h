/** \file
* \copyright
* TaCo - Task-Space Control Library.<br>
* Copyright (c) 2015-2016
*<br>
* TaCo is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*<br>
* TaCo is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*<br>
* You should have received a copy of the Lesser GNU Lesser General Public License
* along with TaCo.  If not, see <http://www.gnu.org/licenses/>.
<br>
Author: Brian Soe <bsoe@stanford.edu>
*/


#ifndef TACO_OPTIMIZER_REGULARIZATION_H_
#define TACO_OPTIMIZER_REGULARIZATION_H_

#include <eigen3/Eigen/Dense>

namespace taco {

/** \brief Solver for regularization problems.
 *
 * Problem is in the form:
 *
 *      Minimize   norm(A x - B)^2 + x W x
 *      Subject to E x - E0 == 0
 *                 I x - I0 <= 0
 *
 * \ingroup optimizer
 */

class Regularization
{
public:

    Regularization() {}
    virtual ~Regularization() {}

    /** \brief Set objective quadratic function, given by norm(A x - B)^2 + x W x
     * \param A Quadratic cost
     * \param B Linear cost
     * \param W Regularization cost
     * \param condition If true, transform the problem such that xA'x+B'x=0 has only solution x=0
     */
    virtual void setObjective(const Eigen::MatrixXd& A,
                              const Eigen::VectorXd& B,
                              const Eigen::MatrixXd& W,
                              bool condition = false) = 0;

    /** \brief Set objective quadratic function, given by norm(A x - B)^2 + x W x
     * \param A Quadratic cost
     * \param B Linear cost
     */
    virtual void setObjective(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B) = 0;

    /** \brief Set objective quadratic function, given by norm(A x - B)^2 + x W x
     * \param W Regularization cost
     */
    virtual void setWeighting(const Eigen::MatrixXd& W) = 0;

    /** \brief Set affine equalities, given by Ex-E0=0
     * \param E Equality matrix
     * \param E0 Equality constant
     * \param condition If true, transform the problem such that rows of E' = rank of E'
     */
    virtual void setEqualities(const Eigen::MatrixXd& E,
                               const Eigen::VectorXd& E0,
                               bool condition = false) = 0;

    /** \brief Set affine inequalities, given by Ix-I0=0
     * \param I Inequality matrix
     * \param I0 Inequality constant
     * \param condition If true, transform the problem such that rows of I' = rank of I'
     */
    virtual void setInequalities(const Eigen::MatrixXd& I,
                                 const Eigen::VectorXd& I0,
                                 bool condition = false) = 0;

    /** \brief Set the conditioning threshold
     * \param epsilon Singular values below epsilon are considered ill-conditioned
     */
    virtual void setConditionNumber(double epsilon) = 0;

    /** \brief Solve the Regularization.
     * \param x_optimal The solution will be returned here. Undefined if no solution is found.
     * \return The cost of the solution or std::numeric_limits::infinity() if the problem is infeasible.
     */
    virtual double solve(Eigen::VectorXd& x_optimal) = 0;

    /** \brief Solve the Regularization with an initial guess.
     * \param x_optimal The solution will be returned here. Undefined if no solution is found.
     * \param x_initial The initial guess of the solution.
     * \return The cost of the solution or std::numeric_limits::infinity() if the problem is infeasible.
     */
    virtual double solve(Eigen::VectorXd& x_optimal, const Eigen::VectorXd x_initial) = 0;

};


} /* namespace taco */
#endif 
