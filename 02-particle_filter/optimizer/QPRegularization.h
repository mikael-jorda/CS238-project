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


#ifndef TACO_OPTIMIZER_QPREGULARIZATION_H_
#define TACO_OPTIMIZER_QPREGULARIZATION_H_

#include <taco/optimizer/Regularization.h>
#include <taco/optimizer/QP.h>
#include <eigen3/Eigen/Dense>
#include <memory>
#include <stdexcept>

namespace taco {

/** \brief Solver for a regularization problems using the template QP solver.
 *
 * \ingroup optimizer
 */

template <class QPT>
class QPRegularization : public Regularization
{
public:

    QPRegularization() {qp_.reset(new QPT());}
    virtual ~QPRegularization() {}

    virtual void setObjective(const Eigen::MatrixXd& A,
                              const Eigen::VectorXd& B,
                              const Eigen::MatrixXd& W,
                              bool condition = false)
    {
        Eigen::MatrixXd A_qp = A.transpose()*A + W;
        Eigen::VectorXd B_qp = -2*B.transpose()*A;
        qp_->setObjective(A_qp, B_qp,condition);
    }

    virtual void setObjective(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B)
    {
        throw std::runtime_error("QPRegularization. setObjective(A,B). Use setObjective(A,B,W) instead.");
    }

    virtual void setWeighting(const Eigen::MatrixXd& W)
    {
        throw std::runtime_error("QPRegularization. setWeighting(). Use setObjective(A,B,W) instead.");
    }

    virtual void setEqualities(const Eigen::MatrixXd& E,
                               const Eigen::VectorXd& E0,
                               bool condition = false)
    {
        qp_->setEqualities(E,E0,condition);
    }

    virtual void setInequalities(const Eigen::MatrixXd& I,
                                 const Eigen::VectorXd& I0,
                                 bool condition = false)
    {
        qp_->setInequalities(I,I0,condition);
    }

    virtual void setConditionNumber(double epsilon)
    {
        qp_->setConditionNumber(epsilon);
    }

    virtual double solve(Eigen::VectorXd& x_optimal)
    {
        return qp_->solve(x_optimal);
    }

    virtual double solve(Eigen::VectorXd& x_optimal, const Eigen::VectorXd x_initial)
    {
        return qp_->solve(x_optimal, x_initial);
    }

    /** \brief The QP used to solve the problem
     * \return The QP
     */
    QPT* qp(){ return dynamic_cast<QPT*>( qp_.get() ); }

protected:

    /// \brief Equivalent QP problem.
    std::unique_ptr<QP> qp_;

};


} /* namespace taco */
#endif 
