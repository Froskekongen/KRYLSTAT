// This file is part of the KRYLSTAT function library
//
// Copyright (C) 2011 Erlend Aune <erlenda@math.ntnu.no>
//
// The KRYLSTAT library is free software; 
// you can redistribute it and/or modify it under the terms 
// of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the 
// License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// The KRYLSTAT library is distributed in the 
// hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU Lesser General Public License
// or the GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// the KRYLSTAT library. If not, see 
// <http://www.gnu.org/licenses/>.

#ifndef MONITOR_H
#define MONITOR_H


#include <limits>
#include <iostream>
#include <iomanip>

// #include <Eigen/Sparse>
// using namespace Eigen;
typedef unsigned int int_type;

template <typename Real>
class default_monitor
{
  public:
    template <typename Vector>
    default_monitor(const Vector& b, int_type iteration_limit = 500, Real relative_tolerance = 1e-5, Real absolute_tolerance = 0)
      : b_norm(b.norm()),
	r_norm(std::numeric_limits<Real>::max()),
        iteration_limit_(iteration_limit),
        iteration_count_(0),
        relative_tolerance_(relative_tolerance),
        absolute_tolerance_(absolute_tolerance)
    {}
    
    void operator++(void) {  ++iteration_count_; }
    
    void reset()
    {
      iteration_count_=0;
      r_norm=std::numeric_limits<Real>::max();
    }
    
      template <typename Vector>
      bool finished(const Vector& r)
      {
	r_norm = r.norm();
	return converged() || iteration_count() >= iteration_limit();
      }
    
      bool converged() const
      {
	return residual_norm() <= tolerance();
      }
    
      Real residual_norm() const { return r_norm; };
    
    
      /*! number of iterations
      */
      int iteration_count() const { return iteration_count_; }

      /*! maximum number of iterations
      */
      int iteration_limit() const { return iteration_limit_; }

      /*! relative tolerance
      */
      Real relative_tolerance() const { return relative_tolerance_; }
    
      /*! absolute tolerance
      */
      Real absolute_tolerance() const { return absolute_tolerance_; }
   
      /*! tolerance
      *
      *  Equal to absolute_tolerance() + relative_tolerance() * ||b||
      *
      */ 
      Real tolerance() const { return absolute_tolerance() + relative_tolerance() * b_norm; }

    protected:
      Real r_norm;
      Real b_norm;
      Real relative_tolerance_;
      Real absolute_tolerance_;

      int_type iteration_limit_;
      int_type iteration_count_;
};


template <typename Real>
class verbose_monitor : public default_monitor<Real>
{
    typedef default_monitor<Real> super;

    public:
    /*! Construct a \p verbose_monitor for a given right-hand-side \p b
     *
     *  The \p verbose_monitor terminates iteration when the residual norm
     *  satisfies the condition
     *       ||b - A x|| <= absolute_tolerance + relative_tolerance * ||b||
     *  or when the iteration limit is reached.
     *
     *  \param b right-hand-side of the linear system A x = b
     *  \param iteration_limit maximum number of solver iterations to allow
     *  \param relative_tolerance determines convergence criteria
     *  \param absolute_tolerance determines convergence criteria
     *
     *  \tparam VectorType vector
     */
    template <typename Vector>
    verbose_monitor(const Vector& b, int_type iteration_limit = 500, Real relative_tolerance = 1e-5, Real absolute_tolerance = 0)
        : super(b, iteration_limit, relative_tolerance, absolute_tolerance)
    {
        std::cout << "Solver will continue until ";
        std::cout << "residual norm " << super::tolerance() << " or reaching ";
        std::cout << super::iteration_limit() << " iterations " << std::endl;
        std::cout << "  Iteration Number  | Residual Norm" << std::endl;
    }
    
    template <typename Vector>
    bool finished(const Vector& r)
    {
        super::r_norm = r.norm();

        std::cout << "       "  << std::setw(10) << super::iteration_count();
        std::cout << "       "  << std::setw(10) << std::scientific << super::residual_norm() << std::endl;

        if (super::converged())
        {
            std::cout << "Successfully converged after " << super::iteration_count() << " iterations." << std::endl;
            return true;
        }
        else if (super::iteration_count() >= super::iteration_limit())
        {
            std::cout << "Failed to converge after " << super::iteration_count() << " iterations." << std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }
};
    
    


#endif
