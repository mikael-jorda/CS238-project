/**
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


#ifndef TACO_OPTIMIZER_QUADPROGRT_CPP_
#define TACO_OPTIMIZER_QUADPROGRT_CPP_

#include "QuadProgRT.h"
#include <taco/controller/math.h>
#include <iostream>

using namespace Eigen;


namespace taco {

QuadProgRT::QuadProgRT()
{
    if (std::numeric_limits<double>::has_infinity) {
        inf_ = std::numeric_limits<double>::infinity();
    } else {
        inf_ = 1.0E300;
        printWarning("No infinity on this machine. Using ["+std::to_string(inf_)+"].");
    }
}

QuadProgRT::~QuadProgRT(){}

void QuadProgRT::setObjective(const Eigen::MatrixXd& A,
                              const Eigen::VectorXd& B,
                              bool condition)
{
    if (condition){
        printWarning("setObjective. No conditioning implemented. Using original problem.");
        obj_A_ = A;
        obj_B_ = B;
    } else {
        obj_A_ = A;
        obj_B_ = B;
    }

    G = 2*obj_A_;
    g0 = obj_B_;

    int n = G.cols();
    R_.resize(n, n);
    J_.resize(n, n);
    z_.resize(n);
    d_.resize(n);
    np_.resize(n);
    x_old_.resize(n);
}

void QuadProgRT::setEqualities(const Eigen::MatrixXd& E,
                               const Eigen::VectorXd& E0,
                               bool condition)
{
    if (condition){
        taco::reducedLeastSquares(E.transpose(), -E0, CE, ce0, epsilon_);
    } else {
        CE = E.transpose();
        ce0 = -E0;
    }
    resizeEI();
}

void QuadProgRT::setInequalities(const Eigen::MatrixXd& I,
                                 const Eigen::VectorXd& I0,
                                 bool condition)
{
    if (condition){
        taco::reducedLeastSquares(-I.transpose(), I0, CI, ci0, epsilon_);
    } else {
        CI = -I.transpose();
        ci0 = I0;
    }
    resizeEI();
}

void QuadProgRT::resizeEI()
{
    int p = CE.cols();
    int m = CI.cols();
    s_.resize(m + p);
    r_.resize(m + p);
    u_.resize(m + p);
    u_old_.resize(m + p);
    A_.resize(m + p);
    A_old_.resize(m + p);
    iai_.resize(m + p);
    iaexcl_.resize(m + p);
}

void QuadProgRT::setConditionNumber(double epsilon)
{
    epsilon_ = epsilon;
}

void QuadProgRT::checkProblem()
{
    std::ostringstream msg;
    {
      //Ensure that the dimensions of the matrices and ublas::vectors can be
      //safely converted from unsigned int into to int without overflow.
      unsigned mx = std::numeric_limits<int>::max();
      if(G.cols() >= mx || G.rows() >= mx ||
         CE.rows() >= mx || CE.cols() >= mx ||
         CI.rows() >= mx || CI.cols() >= mx ||
         ci0.size() >= mx || ce0.size() >= mx || g0.size() >= mx){
        msg << "The dimensions of one of the input matrices or ublas::vectors were "
        << "too large." << std::endl
        << "The maximum allowable size for inputs to solve_quadprog is:"
        << mx << std::endl;
        throw std::logic_error(msg.str());
      }
    }
    int n = G.cols(), p = CE.cols(), m = CI.cols();
    if ((int)G.rows() != n)
    {
      msg << "The ublas::matrix G is not a square ublas::matrix (" << G.rows() << " x "
      << G.cols() << ")";
      throw std::logic_error(msg.str());
    }
    if ((int)CE.rows() != n)
    {
      msg << "The ublas::matrix CE is incompatible (incorrect number of rows "
      << CE.rows() << " , expecting " << n << ")";
      throw std::logic_error(msg.str());
    }
    if ((int)ce0.size() != p)
    {
      msg << "The ublas::vector ce0 is incompatible (incorrect dimension "
      << ce0.size() << ", expecting " << p << ")";
      throw std::logic_error(msg.str());
    }
    if ((int)CI.rows() != n)
    {
      msg << "The ublas::matrix CI is incompatible (incorrect number of rows "
      << CI.rows() << " , expecting " << n << ")";
      throw std::logic_error(msg.str());
    }
    if ((int)ci0.size() != m)
    {
      msg << "The ublas::vector ci0 is incompatible (incorrect dimension "
      << ci0.size() << ", expecting " << m << ")";
      throw std::logic_error(msg.str());
    }
}


void QuadProgRT::printWarning(const std::string &message)
{
    std::cout << "WARNING. QuadProgRT. " << message << '\n';
}


// The Solving function, implementing the Goldfarb-Idnani method
double QuadProgRT::solve(Eigen::VectorXd& x, const Eigen::VectorXd x_initial)
{
    printWarning("solve(x, x_initial). Not implemented, using solve(x) instead.");
    return solve(x);
}


// The Solving function, implementing the Goldfarb-Idnani method
double QuadProgRT::solve(Eigen::VectorXd& x)
{
    // reset because solve changes G and g0
    G = 2*obj_A_;
    g0 = obj_B_;

    int i, j, k, l; /* indices */
    int ip, me, mi;
    int n=g0.size();  int p=ce0.size();  int m=ci0.size();
    //MatrixXd R(G.rows(),G.cols()), J(G.rows(),G.cols());

    LLT<MatrixXd,Lower> chol(G.cols());

    //VectorXd s(m+p), z(n), r(m + p), d(n),  np(n), u(m + p);
    //VectorXd x_old(n), u_old(m + p);
    double f_value, psi, c1, c2, sum, ss, R_norm;
    //const double inf = std::numeric_limits<double>::infinity();
    double t, t1, t2; /* t is the step length, which is the minimum of the partial step length t1
                       * and the full step length t2 */
    //VectorXi A(m + p), A_old(m + p), iai(m + p);
    int q;
    int iq, iter = 0;
    //bool iaexcl[m + p];

    me = p; /* number of equality constraints */
    mi = m; /* number of inequality constraints */
    q = 0;  /* size of the active set A (containing the indices of the active constraints) */

    /*
     * Preprocessing phase
     */

    /* compute the trace of the original matrix G */
    c1 = G.trace();

      /* decompose the matrix G in the form LL^T */
    chol.compute(G);

    /* initialize the matrix R */
    d_.setZero();
    R_.setZero();
      R_norm = 1.0; /* this variable will hold the norm of the matrix R */

      /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
    // J = L^-T
    J_.setIdentity();
    J_ = chol.matrixU().solve(J_);
      c2 = J_.trace();
  #ifdef TRACE_SOLVER
   print_matrix("J", J_, n);
  #endif

      /* c1 * c2 is an estimate for cond(G) */

      /*
     * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g0 x
     * this is a feasible point in the dual space
       * x = G^-1 * g0
     */
    x = chol.solve(g0);
    x = -x;
      /* and compute the current solution value */
      f_value = 0.5 * g0.dot(x);
  #ifdef TRACE_SOLVER
    std::cerr << "Unconstrained solution: " << f_value << std::endl;
    print_vector("x", x, n);
  #endif

      /* Add equality constraints to the working set A */
    iq = 0;
      for (i = 0; i < me; i++)
      {
      np_ = CE.col(i);
      compute_d(d_, J_, np_);
          update_z(z_, J_, d_,  iq);
          update_r(R_, r_, d_,  iq);
  #ifdef TRACE_SOLVER
          print_matrix("R", R_, iq);
          print_vector("z", z_, n);
          print_vector("r", r_, iq);
          print_vector("d", d_, n);
  #endif

      /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
        becomes feasible */
      t2 = 0.0;
      if (::std::abs(z_.dot(z_)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
        t2 = (-np_.dot(x) - ce0(i)) / z_.dot(np_);

      x += t2 * z_;

      /* set u = u+ */
      u_(iq) = t2;
      u_.head(iq) -= t2 * r_.head(iq);

      /* compute the new solution value */
      f_value += 0.5 * (t2 * t2) * z_.dot(np_);
      A_(i) = -i - 1;

      if (!add_constraint(R_, J_, d_, iq, R_norm))
      {
        // FIXME: it should raise an error
        // Equality constraints are linearly dependent
        return f_value;
      }
    }

      /* set iai = K \ A */
      for (i = 0; i < mi; i++)
          iai_(i) = i;

  l1:	iter++;
  #ifdef TRACE_SOLVER
    print_vector("x", x, n);
  #endif
    /* step 1: choose a violated constraint */
      for (i = me; i < iq; i++)
      {
        ip = A_(i);
          iai_(ip) = -1;
      }

      /* compute s(x) = ci^T * x + ci0 for all elements of K \ A */
      ss = 0.0;
      psi = 0.0; /* this value will contain the sum of all infeasibilities */
      ip = 0; /* ip will be the index of the chosen violated constraint */
      for (i = 0; i < mi; i++)
      {
          iaexcl_[i] = true;
          sum = CI.col(i).dot(x) + ci0(i);
          s_(i) = sum;
          psi += std::min(0.0, sum);
      }
  #ifdef TRACE_SOLVER
    print_vector("s", s, mi);
  #endif


      if (::std::abs(psi) <= mi * std::numeric_limits<double>::epsilon() * c1 * c2* 100.0)
      {
      /* numerically there are not infeasibilities anymore */
      q = iq;
          return f_value;
    }

    /* save old values for u, x and A */
     u_old_.head(iq) = u_.head(iq);
     A_old_.head(iq) = A_.head(iq);
     x_old_ = x;

  l2: /* Step 2: check for feasibility and determine a new S-pair */
      for (i = 0; i < mi; i++)
      {
          if (s_(i) < ss && iai_(i) != -1 && iaexcl_[i])
          {
              ss = s_(i);
              ip = i;
          }
      }
    if (ss >= 0.0)
    {
      q = iq;
      return f_value;
    }

    /* set np = n(ip) */
    np_ = CI.col(ip);
    /* set u = (u 0)^T */
    u_(iq) = 0.0;
    /* add ip to the active set A */
    A_(iq) = ip;

  #ifdef TRACE_SOLVER
      std::cerr << "Trying with constraint " << ip << std::endl;
      print_vector("np", np, n);
  #endif

  l2a:/* Step 2a: determine step direction */
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    compute_d(d_, J_, np_);
    update_z(z_, J_, d_, iq);
    /* compute N* np (if q > 0): the negative of the step direction in the dual space */
    update_r(R_, r_, d_, iq);
  #ifdef TRACE_SOLVER
    std::cerr << "Step direction z" << std::endl;
          print_vector("z", z_, n);
          print_vector("r", r_, iq + 1);
      print_vector("u", u_, iq + 1);
      print_vector("d", d_, n);
      print_ivector("A", A_, iq + 1);
  #endif

    /* Step 2b: compute step length */
    l = 0;
    /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
    t1 = inf_; /* +inf */
    /* find the index l s.t. it reaches the minimum of u+(x) / r */
    for (k = me; k < iq; k++)
    {
      double tmp;
      if (r_(k) > 0.0 && ((tmp = u_(k) / r_(k)) < t1) )
      {
        t1 = tmp;
        l = A_(k);
      }
    }
    /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
    if (::std::abs(z_.dot(z_))  > std::numeric_limits<double>::epsilon()) // i.e. z != 0
      t2 = -s_(ip) / z_.dot(np_);
    else
      t2 = inf_; /* +inf */

    /* the step is chosen as the minimum of t1 and t2 */
    t = std::min(t1, t2);
  #ifdef TRACE_SOLVER
    std::cerr << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
  #endif

    /* Step 2c: determine new S-pair and take step: */

    /* case (i): no step in primal or dual space */
    if (t >= inf_)
    {
      /* QPP is infeasible */
      // FIXME: unbounded to raise
      q = iq;
      return inf_;
    }
    /* case (ii): step in dual space */
    if (t2 >= inf_)
    {
      /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
      u_.head(iq) -= t * r_.head(iq);
      u_(iq) += t;
      iai_(l) = l;
      delete_constraint(R_, J_, A_, u_, p, iq, l);
  #ifdef TRACE_SOLVER
      std::cerr << " in dual space: "
        << f_value << std::endl;
      print_vector("x", x, n);
      print_vector("z", z_, n);
          print_ivector("A", A_, iq + 1);
  #endif
      goto l2a;
    }

    /* case (iii): step in primal and dual space */

    x += t * z_;
    /* update the solution value */
    f_value += t * z_.dot(np_) * (0.5 * t + u_(iq));

    u_.head(iq) -= t * r_.head(iq);
    u_(iq) += t;
  #ifdef TRACE_SOLVER
    std::cerr << " in both spaces: "
      << f_value << std::endl;
      print_vector("x", x, n);
      print_vector("u", u_, iq + 1);
      print_vector("r", r_, iq + 1);
      print_ivector("A", A_, iq + 1);
  #endif

    if (t == t2)
    {
  #ifdef TRACE_SOLVER
      std::cerr << "Full step has taken " << t << std::endl;
      print_vector("x", x, n);
  #endif
      /* full step has taken */
      /* add constraint ip to the active set*/
          if (!add_constraint(R_, J_, d_, iq, R_norm))
          {
              iaexcl_[ip] = false;
              delete_constraint(R_, J_, A_, u_, p, iq, ip);
  #ifdef TRACE_SOLVER
        print_matrix("R", R_, n);
        print_ivector("A", A_, iq);
  #endif
              for (i = 0; i < m; i++)
                  iai_(i) = i;
              for (i = 0; i < iq; i++)
              {
                  A_(i) = A_old_(i);
                  iai_(A_(i)) = -1;
                  u_(i) = u_old_(i);
              }
              x = x_old_;
        goto l2; /* go to step 2 */
          }
      else
        iai_(ip) = -1;
  #ifdef TRACE_SOLVER
      print_matrix("R", R_, n);
      print_ivector("A", A_, iq);
  #endif
      goto l1;
    }

    /* a patial step has taken */
  #ifdef TRACE_SOLVER
    std::cerr << "Partial step has taken " << t << std::endl;
    print_vector("x", x, n);
  #endif
    /* drop constraint l */
      iai_(l) = l;
      delete_constraint(R_, J_, A_, u_, p, iq, l);
  #ifdef TRACE_SOLVER
    print_matrix("R", R_, n);
    print_ivector("A", A_, iq);
  #endif

    s_(ip) = CI.col(ip).dot(x) + ci0(ip);

  #ifdef TRACE_SOLVER
    print_vector("s", s_, mi);
  #endif
    goto l2a;
}


// ============================================================================================
void QuadProgRT::compute_d(VectorXd &d, const MatrixXd& J, const VectorXd& np)
{
  d = J.adjoint() * np;
}

void QuadProgRT::update_z(VectorXd& z, const MatrixXd& J, const VectorXd& d,  int iq)
{
  z = J.rightCols(z.size()-iq) * d.tail(d.size()-iq);
}

void QuadProgRT::update_r(const MatrixXd& R, VectorXd& r, const VectorXd& d, int iq)
{
  r.head(iq)= R.topLeftCorner(iq,iq).triangularView<Upper>().solve(d.head(iq));
}


bool QuadProgRT::add_constraint(MatrixXd& R, MatrixXd& J, VectorXd& d, int& iq, double& R_norm)
{
    int n=J.rows();
   #ifdef TRACE_SOLVER
     std::cerr << "Add constraint " << iq << '/';
   #endif
    int i, j, k;
    double cc, ss, h, t1, t2, xny;

    /* we have to find the Givens rotation which will reduce the element
    d(j) to zero.
    if it is already zero we don't have to do anything, except of
    decreasing j */
    for (j = n - 1; j >= iq + 1; j--)
    {
       /* The Givens rotation is done with the matrix (cc cs, cs -cc).
       If cc is one, then element (j) of d is zero compared with element
       (j - 1). Hence we don't have to do anything.
       If cc is zero, then we just have to switch column (j) and column (j - 1)
       of J. Since we only switch columns in J, we have to be careful how we
       update d depending on the sign of gs.
       Otherwise we have to apply the Givens rotation to these columns.
       The i - 1 element of d has to be updated to h. */
       cc = d(j - 1);
       ss = d(j);
       h = distance(cc, ss);
       if (h == 0.0)
           continue;
       d(j) = 0.0;
       ss = ss / h;
       cc = cc / h;
       if (cc < 0.0)
       {
           cc = -cc;
           ss = -ss;
           d(j - 1) = -h;
       }
       else
           d(j - 1) = h;
       xny = ss / (1.0 + cc);
       for (k = 0; k < n; k++)
       {
           t1 = J(k,j - 1);
           t2 = J(k,j);
           J(k,j - 1) = t1 * cc + t2 * ss;
           J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
       }
   }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R */
    R.col(iq-1).head(iq) = d.head(iq);
    #ifdef TRACE_SOLVER
    std::cerr << iq << std::endl;
    #endif

    if (::std::abs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
       return false; // problem degenerate
    R_norm = std::max<double>(R_norm, ::std::abs(d(iq - 1)));
    return true;
}

void QuadProgRT::delete_constraint(MatrixXd& R, MatrixXd& J, VectorXi& A, VectorXd& u, int p, int& iq, int l)
{
    int n = R.rows();
  #ifdef TRACE_SOLVER
    std::cerr << "Delete constraint " << l << ' ' << iq;
  #endif
    int i, j, k, qq;
    double cc, ss, h, xny, t1, t2;

    /* Find the index qq for active constraint l to be removed */
    for (i = p; i < iq; i++)
    if (A(i) == l)
    {
      qq = i;
      break;
    }

    /* remove the constraint from the active set and the duals */
    for (i = qq; i < iq - 1; i++)
    {
      A(i) = A(i + 1);
      u(i) = u(i + 1);
      R.col(i) = R.col(i+1);
    }

    A(iq - 1) = A(iq);
    u(iq - 1) = u(iq);
    A(iq) = 0;
    u(iq) = 0.0;
    for (j = 0; j < iq; j++)
      R(j,iq - 1) = 0.0;
    /* constraint has been fully removed */
    iq--;
  #ifdef TRACE_SOLVER
    std::cerr << '/' << iq << std::endl;
  #endif

    if (iq == 0)
        return;

    for (j = qq; j < iq; j++)
    {
        cc = R(j,j);
        ss = R(j + 1,j);
        h = distance(cc, ss);
        if (h == 0.0)
            continue;
        cc = cc / h;
        ss = ss / h;
        R(j + 1,j) = 0.0;
        if (cc < 0.0)
        {
            R(j,j) = -h;
            cc = -cc;
            ss = -ss;
        }
        else
        R(j,j) = h;

        xny = ss / (1.0 + cc);
        for (k = j + 1; k < iq; k++)
        {
            t1 = R(j,k);
            t2 = R(j + 1,k);
            R(j,k) = t1 * cc + t2 * ss;
            R(j + 1,k) = xny * (t1 + R(j,k)) - t2;
        }
        for (k = 0; k < n; k++)
        {
            t1 = J(k,j);
            t2 = J(k,j + 1);
            J(k,j) = t1 * cc + t2 * ss;
            J(k,j + 1) = xny * (J(k,j) + t1) - t2;
        }
    }
}

double QuadProgRT::distance(double a, double b)
{
  double a1, b1, t;
  a1 = ::std::abs(a);
  b1 = ::std::fabs(b);
  if (a1 > b1) 
  {
    t = (b1 / a1);
    return a1 * ::std::sqrt(1.0 + t * t);
  }
  else
    if (b1 > a1)
    {
      t = (a1 / b1);
      return b1 * ::std::sqrt(1.0 + t * t);
    }
  return a1 * ::std::sqrt(2.0);
}


}

#endif // #define TACO_OPTIMIZER_QUADPROGRT_CPP_
