// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 [Released on march 05, 2007].

 Copyright [c] 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  [the
 "Software"], to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de   $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de   $
//   Date:                $Date:                      Decem 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef NURBS_SURFACE_H
#define NURBS_SURFACE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "utilities/math_utils.h"

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "control_point.h"
#include "b_spline_utilities.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

class NURBSSurface
{
public:
	///@name Type Definitions
	///@{

	// For matrix vector operations
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	typedef typename SparseSpaceType::VectorType VectorType;
	typedef std::vector<int> IntVector;
	typedef std::vector<double> DoubleVector;
	typedef boost::python::extract<double> takeDouble;
	typedef boost::python::extract<int> takeInt;
	typedef boost::python::extract<bool> takeBool;
	typedef std::vector<ControlPoint> ControlPointVector;

	///@}

	/// Pointer definition of NURBSSurface
	//    KRATOS_CLASS_POINTER_DEFINITION[NURBSSurface];

	/// Default constructor.
	NURBSSurface()
	{
	}

	NURBSSurface(DoubleVector knot_vector_u, DoubleVector knot_vector_v, int p, int q, ControlPointVector control_points) 
	: m_knot_vector_u(knot_vector_u),
	  m_knot_vector_v(knot_vector_v),
	  m_p(p),
	  m_q(q),
	  m_control_points(control_points)
	{
		m_n_u = m_knot_vector_u.size() - m_p - 1;
		m_n_v = m_knot_vector_v.size() - m_q - 1;

		if (m_control_points.size() != m_n_u * m_n_v)
		{
			std::cout << "Invalid NURBSSurface" << std::endl;
		}
	}

	/// Destructor.
	virtual ~NURBSSurface()
	{
	}

	//  #####################################################################################
	// #######################################################################################
	///
	///  \details    returns the cartesian coordinates (global) for a specific point 
	///              located on the NURBS surface S(u=fixed and v=fixed)
	///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	///              Algorithm A4.3
	///
	/// ======================================================================================
	///  \param[in]  rSurfacePoint   	evaluated point
	///  \param[in]  _uoi    			local parameter in u-direction
	///  \param[in]  _voi    			local parameter in v-direction
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void EvaluateSurfacePoint(Point<3>& rSurfacePoint, double u, double v)
	{
		rSurfacePoint[0] = 0;
		rSurfacePoint[1] = 0;
		rSurfacePoint[2] = 0;

		int span_u=find_Knot_Span(m_knot_vector_u,u,m_p,m_n_u);
		int span_v=find_Knot_Span(m_knot_vector_v,v,m_q,m_n_v);

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				Matrix R;
				EvaluateNURBSFunctions(span_u, span_v, u, v, R);

				rSurfacePoint[0] += R(b,c) * m_control_points[control_point_index].getX();
				rSurfacePoint[1] += R(b,c) * m_control_points[control_point_index].getY();
				rSurfacePoint[2] += R(b,c) * m_control_points[control_point_index].getZ();
			}
		}
	}

	//  #####################################################################################
	// #######################################################################################
	//
	///  \details    returns the basis functions of NURBS basis function w.r.t. u,v
	///              span_u, span_v are the knot span indices. if unknown, insert 0!
	///
	/// ======================================================================================
	///  \param[in]  span_u     knotspan index in u-direction
	///  \param[in]  span_v     knotspan index in v-direction
	///  \param[in]  _u         local parameter in u-direction
	///  \param[in]  _v         local parameter in v-direction
	///  \param[out] R         basis func
	///
	/// ======================================================================================
	///  \author     Daniel Baumgärtner (12/2016)
	//
	//########################################################################################
	void EvaluateNURBSFunctions(int span_u, int span_v, double _u, double _v, Matrix& R)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		Vector N;
		Vector M;

		R.resize( m_p+1,m_q+1 );
		noalias(R) = ZeroMatrix( m_p+1 ,m_q+1 );

		// Evaluate basis functions with derivatives
		eval_nonzero_basis_function(N, m_knot_vector_u, _u, span_u, m_p);
		eval_nonzero_basis_function(M, m_knot_vector_v, _v, span_v, m_q);

		double sum = 0.0;

		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Evaluate basis function
				R(b,c) = N(b)*M(c)*m_control_points[control_point_index].getWeight();
				sum +=R(b,c);
			}
		}

		// divide by sum only required in terms of rational basis functions
		//if (std::abs(sum-weight)> cepsilon) //Breitenberger 18.06.2014
		double inv_sum = 1/sum;
		// divide through by sum
		for(int c=0;c<=m_q;c++)
			for(int b=0;b<=m_p;b++)
				R(b,c) = inv_sum*R(b,c);
	}

	//  #####################################################################################
	// #######################################################################################
	///
	///  \details    returns the basis fucntions and the first derivative of NURBS basis function w.r.t. u,v
	///              _i,_j are the knot span indices. if unknown, insert 0!
	///
	/// ======================================================================================
	///  \param[in]  _i         knotspan index in u-direction
	///  \param[in]  _u         local parameter in u-direction
	///  \param[in]  _j         knotspan index in v-direction
	///  \param[in]  _v         local parameter in v-direction
	///  \param[out] _R         basis func
	///  \param[out] _dR        1st derivatives
	///
	/// ======================================================================================
	///  \author     from M.Breitenberger in Carat (12/2009)
	//
	//########################################################################################
	void EvaluateNURBSFunctionsAndDerivative( int span_u, int span_v, double _u, double _v, Matrix& R, std::vector<Matrix>& dR)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);
		
		Matrix N_matrix;
		Matrix M_matrix;
		eval_nonzero_basis_function_with_derivatives(N_matrix, m_knot_vector_u, _u, span_u, m_p, 1);
		eval_nonzero_basis_function_with_derivatives(M_matrix, m_knot_vector_v, _v, span_v, m_q, 1);
		double sum = 0.0;
		double dsum1 = 0.0;
		double dsum2 = 0.0;
		double weight;
		
		R.resize(m_p+1,m_q+1);
		dR.resize(2);
		dR[0].resize(m_p+1,m_q+1);
		dR[1].resize(m_p+1,m_q+1);
		
		for (int  c=0;c<=m_q;c++)
		{
			for (int  b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Evaluate basis function
				weight = m_control_points[control_point_index].getWeight();				
				R(b,c) = N_matrix(0,b)*M_matrix(0,c)*weight;
				sum +=R(b,c);
			
				//First derivatives
				dR[0](b,c) = N_matrix(1,b)*M_matrix(0,c)*weight;
				dsum1 += dR[0](b,c);
				dR[1](b,c) = N_matrix(0,b)*M_matrix(1,c)*weight;
				dsum2 += dR[1](b,c);
			}
		}
		
		// divide by sum only required in terms of rational basis functions
		double inv_sum = 1.0/sum;
		// divide through by sum
		for(int  c=0;c<=m_q;c++)
		{
			for(int  b=0;b<=m_p;b++)
			{
				R(b,c) = inv_sum*R(b,c);
				dR[0](b,c) = inv_sum*dR[0](b,c) - R(b,c)*dsum1*inv_sum;
				dR[1](b,c) = inv_sum*dR[1](b,c) - R(b,c)*dsum2*inv_sum;
			}
		}
	}	

	// #######################################################################################
	///
	///  \details    returns the first and second derivative of NURBS basis function w.r.t. u,v
	///              span_u,span_v are the knot span indices. if unknown, insert 0!
	///
	/// ======================================================================================
	///  \param[in]  span_u     knotspan index in u-direction
	///  \param[in]  _u     	local parameter in u-direction
	///  \param[in]  span_v 	knotspan index in v-direction
	///  \param[in]  _v         local parameter in v-direction
	///  \param[out] _dR        1st derivatives
	///  \param[out] _ddR       2nd derivatives
	///
	/// ======================================================================================
	///  \author     from M.Breitenberger in Carat (12/2009)
	//
	//########################################################################################
	
	void EvaluateNURBSFunctionsDerivatives(int span_u,int span_v, double _u, double _v, Matrix& _dR, Matrix& _ddR)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);
	
		int ne = (m_p+1)*(m_q+1); // Control Points per element
		Matrix N;              // Basisfunc at _u
		Matrix M;              // Basisfunc at _v
		eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, _u, span_u, m_p, 2);
		eval_nonzero_basis_function_with_derivatives(M, m_knot_vector_v, _v, span_v, m_q, 2);
		
		vector<double> r(ne);
		r.clear();
		_dR.resize(ne,2);
		_ddR.resize(ne,3);
		double sum = 0.0;
		Vector dsum = ZeroVector(2);
		Vector ddsum = ZeroVector(3);
		double weight;
		
		int k=0;
		for(int c=0;c<=m_q;c++)
		{
			for(int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Evaluate basis function
				weight = m_control_points[control_point_index].getWeight();

				r[k] = N(0,b)*M(0,c)*weight;
				sum +=r[k];     
				//First derivatives
				_dR(k,0) = N(1,b)*M(0,c)*weight;
				dsum[0] +=_dR(k,0);
				_dR(k,1) = N(0,b)*M(1,c)*weight;
				dsum(1) +=_dR(k,1);  
				//Second derivatives  1-du^2, 2-dv^2, 3-dudv
				_ddR(k,0) = N(2,b)*M(0,c)*weight;
				ddsum(0)  = ddsum(0) + _ddR(k,0);
				_ddR(k,1) = N(0,b)*M(2,c)*weight;
				ddsum(1)  = ddsum(1) + _ddR(k,1);
				_ddR(k,2) = N(1,b)*M(1,c)*weight;
				ddsum(2)  = ddsum(2) + _ddR(k,2);
				k++;
			}
		}
		double sum_2 = pow(sum,2);
		double sum_3 = pow(sum,3);
		// divide through by sum
		for(int k=0;k<ne;k++)
		{
		
			_ddR(k,0) = _ddR(k,0)/sum - 2.0*_dR(k,0)*dsum[0]/sum_2
					-r[k]*ddsum[0]/sum_2 + 2.0*r[k]*dsum[0]*dsum[0]/sum_3;
			_ddR(k,1) = _ddR(k,1)/sum - 2.0*_dR(k,1)*dsum[1]/sum_2
					-r[k]*ddsum[1]/sum_2 + 2.0*r[k]*dsum[1]*dsum[1]/sum_3;
			_ddR(k,2) = _ddR(k,2)/sum - _dR(k,0)*dsum[1]/sum_2 - _dR(k,1)*dsum[0]/sum_2
					-r[k]*ddsum[2]/sum_2 + 2.0*r[k]*dsum[0]*dsum[1]/sum_3;
			_dR(k,0) = _dR(k,0)/sum - r[k]*dsum[0]/sum_2;
			_dR(k,1) = _dR(k,1)/sum - r[k]*dsum[1]/sum_2;
		}
	}

	// #####################################################################################
	// #######################################################################################
	///
	///  \details    computes the variation of the local coordinate system with respect to the
	///              displacements in global direction
	///
	///  \return
	///  \param[in]   _u       parameter in u-direction
	///  \param[in]   _v       parameter in v-direction
	///  \param[in]   _par_g1  base vector in parameter space for the boundary curve -> t2
	///
	///  \author     from M.Breitenberger in Carat (04/2014)
	
	// ########################################################################################
	void ComputeVariationOfLocalCSY( double _u, double _v,
									 Vector&_par_g1, // entspricht t1 und t2 vom JSON File
									 Vector& _t1, //nicht relevant dummy
									 Vector& _t2, //Dein groß T2
									 Vector& _t3, //Dein T3
									 std::vector<Vector>& _t1r, //nicht relevant dummy
									 std::vector<Vector>& _t2r, //nicht relevant dummy
									 std::vector<Vector>& _t3r ) //t3,r
	{
		//number of dofs affecting _u and _v
		unsigned int number_of_affected_control_points = (m_p+1)*(m_q+1);
		unsigned int number_of_dofs = 3*number_of_affected_control_points;
		_t1r.resize(number_of_dofs);
		_t2r.resize(number_of_dofs);
		_t3r.resize(number_of_dofs);		
		
		Matrix R_matrix;
		std::vector<Matrix> dR_matrix;

		EvaluateNURBSFunctionsAndDerivative(-1,-1,_u,_v,R_matrix, dR_matrix);

		Vector dR1 = ZeroVector(number_of_affected_control_points); //derivatives of shape functions in u-direction
		Vector dR2 = ZeroVector(number_of_affected_control_points); //derivatives of shape functions in v-direction
		unsigned int counter=0;
		for(unsigned int j=0;j<R_matrix.size2();j++)
		{
			for(unsigned int i=0;i<R_matrix.size1();i++)
			{
				dR1(counter) = dR_matrix[0](i,j);
				dR2(counter) = dR_matrix[1](i,j);
				counter++;
			}	
		}	
		
		// computer base vectors
		Matrix g_matrix = GetBaseVectors(-1,-1,_u,_v);
		Vector g1_act = ZeroVector(3);
		g1_act(0) = g_matrix(0,0);
		g1_act(1) = g_matrix(1,0);
		g1_act(2) = g_matrix(2,0);
		Vector g2_act = ZeroVector(3);
		g2_act(0) = g_matrix(0,1);
		g2_act(1) = g_matrix(1,1);
		g2_act(2) = g_matrix(2,1);			
		
		//compute base vectors in actual configuration
		
		Vector tilde_t2 = g1_act*_par_g1(0) + g2_act*_par_g1(1);
		double l_t2 = norm_2(tilde_t2);
		_t2 = tilde_t2/l_t2;
		
		Vector tilde_t3 = MathUtils<double>::CrossProduct(g1_act,g2_act);
		double l_t3 = norm_2(tilde_t3);
		_t3 = tilde_t3/l_t3;
		
		Vector tilde_t1 = MathUtils<double>::CrossProduct(tilde_t2,tilde_t3);
		double l_t1 = norm_2(tilde_t1);
		_t1 = tilde_t1/l_t1;
		
		for(unsigned int  k=0;k<number_of_dofs;k++)
		{
			unsigned int  xyz_k = k%3; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
			unsigned int  i = k/3;     // index for the shape functions
		
			//variations of the base vectors
			Vector a1_r = ZeroVector(3);
			Vector a2_r = ZeroVector(3);
			a1_r(xyz_k) = dR1(i);
			a2_r(xyz_k) = dR2(i);
		
			//variation of the non normalized local vector
			Vector tilde_2_r = _par_g1(0)*a1_r + _par_g1(1)*a2_r;
			Vector tilde_3_r = MathUtils<double>::CrossProduct(a1_r,g2_act) + MathUtils<double>::CrossProduct(g1_act,a2_r);
			Vector tilde_1_r = MathUtils<double>::CrossProduct(tilde_2_r,tilde_t3) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_r);
		
			double line_t1_r = inner_prod(_t1,tilde_1_r);
			double line_t2_r = inner_prod(_t2,tilde_2_r);
			double line_t3_r = inner_prod(_t3,tilde_3_r);
		
			_t1r[k] = tilde_1_r/l_t1 - line_t1_r*_t1/l_t1;
			_t2r[k] = tilde_2_r/l_t2 - line_t2_r*_t2/l_t2;
			_t3r[k] = tilde_3_r/l_t3 - line_t3_r*_t3/l_t3;
		}
	}	

	// #######################################################################################
	///
	///  \details    computes the variation of the local coordinate system with respect to the
	///              displacements in global direction
	///
	///  \return
	///  \param[in]   _u       parameter in u-direction
	///  \param[in]   _v       parameter in v-direction
	///  \param[in]   _par_g1  base vector in parameter space for the boundary curve -> t2
	///
	///  \author     from M.Breitenberger in Carat (04/2014)
	//
	//########################################################################################
	void ComputeSecondVariationOfLocalCSY( double _u, double _v,
									 	   Vector &_par_g1, // entspricht t1 und t2 vom JSON File
									  	   Vector& _t1, 
										   Vector& _t2, 
										   Vector& _t3, 
										   Vector& _t1_der, // necessary for G2 only (not tested yet by Breitenberger!)
										   Vector& _t2_der, // necessary for G2 only (not tested yet by Breitenberger!)
										   Vector& _t3_der, // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<Vector >& _t1_r,
										   std::vector<Vector >& _t2_r,
										   std::vector<Vector >& _t3_r,
										   std::vector<Vector >& _t1_der_r, // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<Vector >& _t2_der_r, // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<Vector >& _t3_der_r, // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<std::vector<Vector > >& _t1_rs,
										   std::vector<std::vector<Vector > >& _t2_rs,
										   std::vector<std::vector<Vector > >& _t3_rs,
										   std::vector<std::vector<Vector > >& _t1_der_rs,  // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<std::vector<Vector > >& _t2_der_rs,  // necessary for G2 only (not tested yet by Breitenberger!)
										   std::vector<std::vector<Vector > >& _t3_der_rs ) // necessary for G2 only (not tested yet by Breitenberger!)
	{
		//number of dofs (just displacements)
		int number_of_affected_control_points = (m_p+1)*(m_q+1);
		int ndofs = 3*number_of_affected_control_points;	
		_t1_r.resize(ndofs);
		_t2_r.resize(ndofs);
		_t3_r.resize(ndofs);
		_t1_der_r.resize(ndofs);
		_t2_der_r.resize(ndofs);
		_t3_der_r.resize(ndofs);
		
		_t1_rs.resize(ndofs);
		_t2_rs.resize(ndofs);
		_t3_rs.resize(ndofs);
		_t1_der_rs.resize(ndofs);
		_t2_der_rs.resize(ndofs);
		_t3_der_rs.resize(ndofs);
		
		for (int i=0; i<ndofs;i++)
		{
			_t1_rs[i].resize(ndofs);
			_t2_rs[i].resize(ndofs);
			_t3_rs[i].resize(ndofs);
			_t1_der_rs[i].resize(ndofs);
			_t2_der_rs[i].resize(ndofs);
			_t3_der_rs[i].resize(ndofs);
		}
		
		Matrix R_matrix;
		std::vector<Matrix> dr_vector_format; //1st Derivatives of shape functions in vector format

		EvaluateNURBSFunctionsAndDerivative(-1,-1,_u,_v,R_matrix, dr_vector_format);

		Vector dr1 = ZeroVector(number_of_affected_control_points); // 1st derivative of shape functions in u-direction
		Vector dr2 = ZeroVector(number_of_affected_control_points); // 1st derivative of shape functions in v-direction
		unsigned int counter=0;
		for(unsigned int j=0;j<R_matrix.size2();j++)
		{
			for(unsigned int i=0;i<R_matrix.size1();i++)
			{
				dr1(counter) = dr_vector_format[0](i,j);
				dr2(counter) = dr_vector_format[1](i,j);
				counter++;
			}	
		}	

		Matrix dr; //1st Derivatives of shape functions in matrix format
		Matrix ddr; //2nd Derivatives of shape functions in vector format
		EvaluateNURBSFunctionsDerivatives(-1,-1, _u, _v, dr,ddr);
		
		// computer base vectors
		Vector g1_act = ZeroVector(3);
		Vector g2_act = ZeroVector(3);
		Vector g3_act = ZeroVector(3);
		Matrix gi_deri = ZeroMatrix(3,3);
		// computer base vectors derived ;
		Vector g1_der_1 = ZeroVector(3);
		Vector g1_der_2 = ZeroVector(3);
		Vector g2_der_1 = ZeroVector(3);
		Vector g2_der_2 = ZeroVector(3);
		// computer base vectors derived wrt tilde_theta;
		Vector g1_der_act = ZeroVector(3);
		Vector g2_der_act = ZeroVector(3);
		Vector g3_der_act = ZeroVector(3);
		
		
		ComputeBaseVectorsAndDerivatives(_u,_v,g1_act,g2_act,g3_act,gi_deri);
		for(int i = 0;i<3;i++)
		{
			g1_der_1(i)=gi_deri(i,0);
			g2_der_2(i)=gi_deri(i,1);
			g1_der_2(i)=gi_deri(i,2);
		}
		g2_der_1=g1_der_2;

		g1_der_act = g1_der_1*_par_g1(0) + g1_der_2*_par_g1(1);
		g2_der_act = g2_der_1*_par_g1(0) + g2_der_2*_par_g1(1);

		//compute base vectors in actual configuration
		
		Vector tilde_t2 = g1_act*_par_g1(0) + g2_act*_par_g1(1);
		double l_t2 = norm_2(tilde_t2);
		_t2 = tilde_t2/l_t2;
		
		Vector tilde_t3 = MathUtils<double>::CrossProduct(g1_act,g2_act);
		double l_t3 = norm_2(tilde_t3);
		_t3 = tilde_t3/l_t3;
		
		Vector tilde_t1 = MathUtils<double>::CrossProduct(tilde_t2,tilde_t3);
		double l_t1 = norm_2(tilde_t1);
		_t1 = tilde_t1/l_t1;
		
		//compute base vectors in actual configuration
		
		Vector tilde_t2_der = g1_der_act*_par_g1(0) + g2_der_act*_par_g1(1);
		_t2_der = tilde_t2_der/l_t2-tilde_t2*inner_prod(tilde_t2_der,tilde_t2)/pow(l_t2,3);
		
		Vector tilde_t3_der = MathUtils<double>::CrossProduct(g1_der_act,g2_act)+MathUtils<double>::CrossProduct(g1_act,g2_der_act);
		_t3_der = tilde_t3_der/l_t3-tilde_t3*inner_prod(tilde_t3_der,tilde_t3)/pow(l_t3,3);
		
		Vector tilde_t1_der = MathUtils<double>::CrossProduct(tilde_t2_der,tilde_t3)+MathUtils<double>::CrossProduct(tilde_t2,tilde_t3_der);
		_t1_der = tilde_t1_der/l_t1-tilde_t1*inner_prod(tilde_t1_der,tilde_t1)/pow(l_t1,3);
		
		//variations of the base vectors
		Vector a1_r = ZeroVector(3);
		Vector a2_r = ZeroVector(3);
		Vector a1_der_r = ZeroVector(3);
		Vector a2_der_r = ZeroVector(3);
		//variations of the base vectors
		Vector a1_s = ZeroVector(3);
		Vector a2_s = ZeroVector(3);
		Vector a1_der_s = ZeroVector(3);
		Vector a2_der_s = ZeroVector(3);
		
		for(int r=0;r<ndofs;r++)
		{
			int xyz_r = r%3; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
			int i = r/3;     // index for the shape functions
		
			a1_r.clear();
			a2_r.clear();
			a1_der_r.clear();
			a2_der_r.clear();
			a1_r(xyz_r) = dr1(i);
			a2_r(xyz_r) = dr2(i);
			a1_der_r(xyz_r) = ddr(i,0)*_par_g1(0)+ddr(i,2)*_par_g1(1);
			a2_der_r(xyz_r) = ddr(i,2)*_par_g1(0)+ddr(i,1)*_par_g1(1);
		
			//variation of the non normalized local vector
			Vector tilde_2_r = _par_g1(0)*a1_r + _par_g1(1)*a2_r;
			Vector tilde_3_r = MathUtils<double>::CrossProduct(a1_r,g2_act) + MathUtils<double>::CrossProduct(g1_act,a2_r);
			Vector tilde_1_r = MathUtils<double>::CrossProduct(tilde_2_r,tilde_t3) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_r);
			Vector tilde_2_der_r = _par_g1(0)*a1_der_r + _par_g1(1)*a2_der_r;
			Vector tilde_3_der_r = MathUtils<double>::CrossProduct(a1_der_r,g2_act) + MathUtils<double>::CrossProduct(g1_der_act,a2_r)+MathUtils<double>::CrossProduct(a1_r,g2_der_act) + MathUtils<double>::CrossProduct(g1_act,a2_der_r);
			Vector tilde_1_der_r = MathUtils<double>::CrossProduct(tilde_2_der_r,tilde_t3) + MathUtils<double>::CrossProduct(tilde_t2_der,tilde_3_r)+MathUtils<double>::CrossProduct(tilde_2_r,tilde_t3_der) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_der_r);
		
			double line_t1_r = inner_prod(_t1,tilde_1_r);
			double line_t2_r = inner_prod(_t2,tilde_2_r);
			double line_t3_r = inner_prod(_t3,tilde_3_r);
			double line_tilde_t1_r = inner_prod(tilde_t1,tilde_1_r);
			double line_tilde_t2_r = inner_prod(tilde_t2,tilde_2_r);
			double line_tilde_t3_r = inner_prod(tilde_t3,tilde_3_r);
			double line_tilde_t1_der_r = inner_prod(tilde_t1_der,tilde_1_r) + inner_prod(tilde_t1,tilde_1_der_r);
			double line_tilde_t2_der_r = inner_prod(tilde_t2_der,tilde_2_r) + inner_prod(tilde_t2,tilde_2_der_r);
			double line_tilde_t3_der_r = inner_prod(tilde_t3_der,tilde_3_r) + inner_prod(tilde_t3,tilde_3_der_r);
		
		
			std::vector<Vector > tilde_2_rs;
			std::vector<Vector > tilde_3_rs;
			std::vector<Vector > tilde_1_rs;
			tilde_2_rs.resize(ndofs);
			tilde_3_rs.resize(ndofs);
			tilde_1_rs.resize(ndofs);
			std::vector<Vector > tilde_2_der_rs;
			std::vector<Vector > tilde_3_der_rs;
			std::vector<Vector > tilde_1_der_rs;
			tilde_2_der_rs.resize(ndofs);
			tilde_3_der_rs.resize(ndofs);
			tilde_1_der_rs.resize(ndofs);
		
			_t1_r[r] = tilde_1_r/l_t1 - line_t1_r*_t1/l_t1;
			_t2_r[r] = tilde_2_r/l_t2 - line_t2_r*_t2/l_t2;
			_t3_r[r] = tilde_3_r/l_t3 - line_t3_r*_t3/l_t3;
		
			_t1_der_r[r] = tilde_1_der_r/l_t1 - (line_tilde_t1_r*tilde_t1_der)/pow(l_t1,3)-(tilde_1_r*inner_prod(tilde_t1,tilde_t1_der)+tilde_t1*line_tilde_t1_der_r)/pow(l_t1,3)+3*(tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r)/pow(l_t1,5);
			_t2_der_r[r] = tilde_2_der_r/l_t2 - (line_tilde_t2_r*tilde_t2_der)/pow(l_t2,3)-(tilde_2_r*inner_prod(tilde_t2,tilde_t2_der)+tilde_t2*line_tilde_t2_der_r)/pow(l_t2,3)+3*(tilde_t2*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r)/pow(l_t2,5);
			_t3_der_r[r] = tilde_3_der_r/l_t3 - (line_tilde_t3_r*tilde_t3_der)/pow(l_t3,3)-(tilde_3_r*inner_prod(tilde_t3,tilde_t3_der)+tilde_t3*line_tilde_t3_der_r)/pow(l_t3,3)+3*(tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r)/pow(l_t3,5);
		
			for (int s=0; s<ndofs; s++)
			{
				int xyz_s = s%3; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
				int j = s/3;     // index for the shape functions
			
				a1_s.clear();
				a2_s.clear();
				a1_der_s.clear();
				a2_der_s.clear();
				a1_s(xyz_s) = dr1(j);
				a2_s(xyz_s) = dr2(j);
				a1_der_s(xyz_s) = ddr(j,0)*_par_g1(0)+ddr(j,2)*_par_g1(1);
				a2_der_s(xyz_s) = ddr(j,2)*_par_g1(0)+ddr(j,1)*_par_g1(1);
			
				//variation of the non normalized local vector
				Vector tilde_2_s = _par_g1(0)*a1_s + _par_g1(1)*a2_s;
				Vector tilde_3_s = MathUtils<double>::CrossProduct(a1_s,g2_act) + MathUtils<double>::CrossProduct(g1_act,a2_s);
				Vector tilde_1_s = MathUtils<double>::CrossProduct(tilde_2_s,tilde_t3) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_s);
				Vector tilde_2_der_s = _par_g1(0)*a1_der_s + _par_g1(1)*a2_der_s;
				Vector tilde_3_der_s = MathUtils<double>::CrossProduct(a1_der_s,g2_act) + MathUtils<double>::CrossProduct(g1_der_act,a2_s)+MathUtils<double>::CrossProduct(a1_s,g2_der_act) + MathUtils<double>::CrossProduct(g1_act,a2_der_s);
				Vector tilde_1_der_s = MathUtils<double>::CrossProduct(tilde_2_der_s,tilde_t3) + MathUtils<double>::CrossProduct(tilde_t2_der,tilde_3_s)+MathUtils<double>::CrossProduct(tilde_2_s,tilde_t3_der) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_der_s);
			
				//tilde_2_rs[s]=0;
				tilde_3_rs[s] = MathUtils<double>::CrossProduct(a1_r,a2_s) + MathUtils<double>::CrossProduct(a1_s,a2_r);
				tilde_1_rs[s] = MathUtils<double>::CrossProduct(tilde_2_s,tilde_3_r) + MathUtils<double>::CrossProduct(tilde_2_r,tilde_3_s) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_rs[s]);
				//tilde_2_der_rs[s]=0;
				tilde_3_der_rs[s] = MathUtils<double>::CrossProduct(a1_der_s,a2_r) + MathUtils<double>::CrossProduct(a1_der_r,a2_s)+MathUtils<double>::CrossProduct(a1_s,a2_der_r) + MathUtils<double>::CrossProduct(a1_r,a2_der_s);
				tilde_1_der_rs[s] = MathUtils<double>::CrossProduct(tilde_2_der_s,tilde_3_r) + MathUtils<double>::CrossProduct(tilde_2_der_r,tilde_3_s) + MathUtils<double>::CrossProduct(tilde_t2_der,tilde_3_rs[s]) + MathUtils<double>::CrossProduct(tilde_2_s,tilde_3_der_r) + MathUtils<double>::CrossProduct(tilde_2_r,tilde_3_der_s) + MathUtils<double>::CrossProduct(tilde_t2,tilde_3_der_rs[s]);
			
				double line_tilde_t1_s = inner_prod(tilde_t1,tilde_1_s);
				double line_tilde_t2_s = inner_prod(tilde_t2,tilde_2_s);
				double line_tilde_t3_s = inner_prod(tilde_t3,tilde_3_s);
				double line_tilde_t1_der_s = inner_prod(tilde_t1_der,tilde_1_s) + inner_prod(tilde_t1,tilde_1_der_s);
				// double line_tilde_t2_der_s = inner_prod(tilde_t2_der,tilde_2_s) + inner_prod(tilde_t2,tilde_2_der_s);
				double line_tilde_t3_der_s = inner_prod(tilde_t3_der,tilde_3_s) + inner_prod(tilde_t3,tilde_3_der_s);
				double line_tilde_t1_rs = inner_prod(tilde_1_r,tilde_1_s);
				double line_tilde_t2_rs = inner_prod(tilde_2_r,tilde_2_s);
				double line_tilde_t3_rs = inner_prod(tilde_3_r,tilde_3_s);
				double line_tilde_t1_der_rs = inner_prod(tilde_1_der_r,tilde_1_s) + inner_prod(tilde_1_r,tilde_1_der_s);
				double line_tilde_t2_der_rs = inner_prod(tilde_2_der_r,tilde_2_s) + inner_prod(tilde_2_r,tilde_2_der_s);
				double line_tilde_t3_der_rs = inner_prod(tilde_3_der_r,tilde_3_s) + inner_prod(tilde_3_r,tilde_3_der_s);
				_t1_rs[r][s] = tilde_1_rs[s]/l_t1 - line_tilde_t1_s*tilde_1_r/pow(l_t1,3) - (line_tilde_t1_rs*tilde_t1+line_tilde_t1_r*tilde_1_s)/pow(l_t1,3) + 3*line_tilde_t1_r*line_tilde_t1_s*tilde_t1/pow(l_t1,5);
				_t2_rs[r][s] = tilde_2_rs[s]/l_t2 - line_tilde_t2_s*tilde_2_r/pow(l_t2,3) - (line_tilde_t2_rs*tilde_t2+line_tilde_t2_r*tilde_2_s)/pow(l_t2,3) + 3*line_tilde_t2_r*line_tilde_t2_s*tilde_t2/pow(l_t2,5);
				_t3_rs[r][s] = tilde_3_rs[s]/l_t3 - line_tilde_t3_s*tilde_3_r/pow(l_t3,3) - (line_tilde_t3_rs*tilde_t3+line_tilde_t3_r*tilde_3_s)/pow(l_t3,3) + 3*line_tilde_t3_r*line_tilde_t3_s*tilde_t3/pow(l_t3,5);
			
				//_t1_der_rs[r][s] = tilde_1_der_rs[s]/l_t1 - tilde_1_der_r*line_tilde_t1_s/pow(l_t1,3) - (line_tilde_t1_rs*tilde_t1_der + line_tilde_t1_r*tilde_1_der_s)/pow(l_t1,3) + 3*(line_tilde_t1_r*tilde_t1_der*line_tilde_t1_s)/pow(l_t1,5) -(tilde_1_rs[s]*inner_prod(tilde_t1,tilde_t1_der)+tilde_1_s*line_tilde_t1_der_r+tilde_1_r*(inner_prod(tilde_1_s,tilde_t1_der)+inner_prod(tilde_t1,tilde_1_der_s))+tilde_t1*line_tilde_t1_der_rs)/pow(l_t1,3) + 3* (tilde_1_r*inner_prod(tilde_t1,tilde_t1_der)+tilde_t1*line_tilde_t1_der_r)*line_tilde_t2_s/pow(l_t1,5) +3*(tilde_1_s*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r+tilde_t1*(inner_prod(tilde_1_s,tilde_t1_der)+inner_prod(tilde_t1,tilde_1_der_s))*line_tilde_t1_r+tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_rs)/pow(l_t1,5)-15*(tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r*line_tilde_t1_s)/pow(l_t1,7);
				_t2_der_rs[r][s] = tilde_2_der_rs[s]/l_t2 - tilde_2_der_r*line_tilde_t2_s/pow(l_t2,3) - (line_tilde_t2_rs*tilde_t2_der + line_tilde_t2_r*tilde_2_der_s)/pow(l_t2,3) + 3*(line_tilde_t2_r*tilde_t2_der*line_tilde_t2_s)/pow(l_t2,5) -(tilde_2_rs[s]*inner_prod(tilde_t2,tilde_t2_der)+tilde_2_s*line_tilde_t2_der_r+tilde_2_r*(inner_prod(tilde_2_s,tilde_t2_der)+inner_prod(tilde_t2,tilde_2_der_s))+tilde_t2*line_tilde_t2_der_rs)/pow(l_t2,3) + 3* (tilde_2_r*inner_prod(tilde_t2,tilde_t2_der)+tilde_t2*line_tilde_t2_der_r)*line_tilde_t2_s/pow(l_t2,5) +3*(tilde_2_s*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r+tilde_t2*(inner_prod(tilde_2_s,tilde_t2_der)+inner_prod(tilde_t2,tilde_2_der_s))*line_tilde_t2_r+tilde_t2*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_rs)/pow(l_t2,5)-25*(tilde_t2*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r*line_tilde_t2_s)/pow(l_t2,7);
				//_t3_der_rs[r][s] = tilde_3_der_rs[s]/l_t3 - tilde_3_der_r*line_tilde_t3_s/pow(l_t3,3) - (line_tilde_t3_rs*tilde_t3_der + line_tilde_t3_r*tilde_3_der_s)/pow(l_t3,3) + 3*(line_tilde_t3_r*tilde_t3_der*line_tilde_t3_s)/pow(l_t3,5) /*-(tilde_3_rs[s]*inner_prod(tilde_t3,tilde_t3_der)+tilde_3_s*line_tilde_t3_der_r+tilde_3_r*(inner_prod(tilde_3_s,tilde_t3_der)+inner_prod(tilde_t3,tilde_3_der_s))+tilde_t3*line_tilde_t3_der_rs)/pow(l_t3,3) + 3* (tilde_3_r*inner_prod(tilde_t3,tilde_t3_der)+tilde_t3*line_tilde_t3_der_r)*line_tilde_t2_s/pow(l_t3,5) +3*(tilde_3_s*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r+tilde_t3*(inner_prod(tilde_3_s,tilde_t3_der)+inner_prod(tilde_t3,tilde_3_der_s))*line_tilde_t3_r+tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_rs)/pow(l_t3,5)-35*(tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r*line_tilde_t3_s)/pow(l_t3,7)*/;
				_t1_der_rs[r][s] = tilde_1_der_rs[s]/l_t1 - tilde_1_der_r*line_tilde_t1_s/pow(l_t1,3) - (tilde_1_der_s*line_tilde_t1_r + tilde_t1_der * (inner_prod(tilde_1_s,tilde_1_r)+inner_prod(tilde_t1,tilde_1_rs[s])))/pow(l_t1,3)+3*(tilde_t1_der*line_tilde_t1_r*line_tilde_t1_s)/pow(l_t1,5)
									- (tilde_1_rs[s] * inner_prod(tilde_t1,tilde_t1_der)+tilde_1_r*line_tilde_t1_der_s+tilde_1_s*line_tilde_t1_der_r)/pow(l_t1,3)
									- (tilde_t1*(inner_prod(tilde_1_rs[s],tilde_t1_der)+inner_prod(tilde_1_der_rs[s],tilde_t1)+line_tilde_t1_der_rs))*pow(l_t1,3)
									+ 3*((tilde_1_r * inner_prod(tilde_t1,tilde_t1_der)+tilde_t1*line_tilde_t1_der_r)*line_tilde_t1_s)/pow(l_t1,5)
									+3*(tilde_1_s * inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r + tilde_t1*line_tilde_t1_der_s*line_tilde_t1_r+tilde_t1 * inner_prod(tilde_t1,tilde_t1_der)*(inner_prod(tilde_1_s, tilde_1_r)+inner_prod(tilde_t1,tilde_1_rs[s])))/pow(l_t1,5)
									- 15* (tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r*line_tilde_t1_s)/pow(l_t1,7);
				_t3_der_rs[r][s] = tilde_3_der_rs[s]/l_t3 - tilde_3_der_r*line_tilde_t3_s/pow(l_t3,3) - (tilde_3_der_s*line_tilde_t3_r + tilde_t3_der * (inner_prod(tilde_3_s,tilde_3_r)+inner_prod(tilde_t3,tilde_3_rs[s])))/pow(l_t3,3)+3*(tilde_t3_der*line_tilde_t3_r*line_tilde_t3_s)/pow(l_t3,5)
									- (tilde_3_rs[s] * inner_prod(tilde_t3,tilde_t3_der)+tilde_3_r*line_tilde_t3_der_s+tilde_3_s*line_tilde_t3_der_r)/pow(l_t3,3)
									- (tilde_t3*(inner_prod(tilde_3_rs[s],tilde_t3_der)+inner_prod(tilde_3_der_rs[s],tilde_t3)+line_tilde_t3_der_rs))*pow(l_t3,3)
									+ 3*((tilde_3_r * inner_prod(tilde_t3,tilde_t3_der)+tilde_t3*line_tilde_t3_der_r)*line_tilde_t3_s)/pow(l_t3,5)
									+3*(tilde_3_s * inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r + tilde_t3*line_tilde_t3_der_s*line_tilde_t3_r+tilde_t3 * inner_prod(tilde_t3,tilde_t3_der)*(inner_prod(tilde_3_s, tilde_3_r)+inner_prod(tilde_t3,tilde_3_rs[s])))/pow(l_t3,5)
									- 15* (tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r*line_tilde_t3_s)/pow(l_t3,7);
			}
		}
	}	

	// #######################################################################################
	///
	///  \details    evaluate Hessian and Gradiend modifying the input objects
	///
	/// ======================================================================================
	///  \param[in]  QminP    	 	Distance Vector
	///  \param[in]  H		     	Hessian reference	
	///  \param[in]  Gradient    	Gradient reference
	///  \param[in]  v    			parameter
	///  \param[in]  u 				parameter 
	///
	/// ======================================================================================
	///  \author     Giovanni Filomeno (1/2017) && Massimo Sferza (1/2017)
	//
	//########################################################################################	

	void EvaluateGradientsForClosestPointSearch(Vector QminP, Matrix& Hessian, Vector& Gradient , double u, double v)
	{
		// The derivatives of the basis functions are evaluated
		Matrix dR;
		Matrix ddR;
		EvaluateNURBSFunctionsDerivatives(-1,-1, u, v, dR,ddR);

		// The derivatives of Q(u,v) are evaluated
		Vector dQdu = ZeroVector(3);
		Vector dQdv = ZeroVector(3);
		Vector dQdudu = ZeroVector(3);
		Vector dQdvdv = ZeroVector(3);
		Vector dQdudv = ZeroVector(3);

		int span_u = find_Knot_Span(m_knot_vector_u,u,m_p,m_n_u);
		int span_v = find_Knot_Span(m_knot_vector_v,v,m_q,m_n_v);		

		int k=0;
		for( int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;	

				double cp_x = m_control_points[control_point_index].getX();			
				double cp_y = m_control_points[control_point_index].getY();
				double cp_z = m_control_points[control_point_index].getZ();

				dQdu(0) += dR(k,0) * cp_x;
				dQdu(1) += dR(k,0) * cp_y;
				dQdu(2) += dR(k,0) * cp_z;

				dQdv(0) += dR(k,1) * cp_x;
				dQdv(1) += dR(k,1) * cp_y;
				dQdv(2) += dR(k,1) * cp_z;

				dQdudu(0) += ddR(k,0) * cp_x;
				dQdudu(1) += ddR(k,0) * cp_y;
				dQdudu(2) += ddR(k,0) * cp_z;

				dQdvdv(0) += ddR(k,1) * cp_x;
				dQdvdv(1) += ddR(k,1) * cp_y;
				dQdvdv(2) += ddR(k,1) * cp_z;

				dQdudv(0) += ddR(k,2) * cp_x;
				dQdudv(1) += ddR(k,2) * cp_y;
				dQdudv(2) += ddR(k,2) * cp_z;

				k++;
			}
		}
		// Hessian and gradient are evaluated
		Hessian(0,0) = 2*(inner_prod(dQdudu,QminP) +  inner_prod(dQdu, dQdu));
		Hessian(0,1) = 2*(inner_prod(dQdudv,QminP) +  inner_prod(dQdu, dQdv)); 
		Hessian(1,0) = 2*(inner_prod(dQdudv,QminP) +  inner_prod(dQdu, dQdv));
		Hessian(1,1) = 2*(inner_prod(dQdvdv,QminP) +  inner_prod(dQdv, dQdv));

		Gradient(0) = 2*inner_prod(dQdu, QminP);
		Gradient(1) = 2*inner_prod(dQdv, QminP);
	}	

	//  #####################################################################################
	// #######################################################################################
	//#
	///   \details     returns the base vectors at point _u, _v               
	///
	/// ======================================================================================
	///  \param[in]  _i         knotspan index in u-direction
	///  \param[in]  _u         parameter in u-direction
	///  \param[in]  _j         knotspan index in v-direction
	///  \param[in]  _v         parameter in u-direction
	///
	/// ======================================================================================
	///  \author     M.Breitenberger (01/2017)
	//
	//########################################################################################
	Matrix GetBaseVectors(int span_u, int span_v, double _u, double _v)
	{
		if(span_u==-1) span_u = find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v = find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);
		Matrix N;
		Matrix M;
		eval_nonzero_basis_function_with_derivatives(N,m_knot_vector_u,_u,span_u,m_p,1);
		eval_nonzero_basis_function_with_derivatives(M,m_knot_vector_v,_v,span_v,m_q,1);
		double sum = 0.0;
		vector<double> dsum(2);
		dsum.clear();
		Matrix R(m_p+1,m_q+1);
		Matrix dR1(m_p+1,m_q+1);
		Matrix dR2(m_p+1,m_q+1);
		double weight;
		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				weight=m_control_points[control_point_index].getWeight();
				R(b,c) = N(0,b)*M(0,c)*weight;
				sum +=R(b,c);
						
				// derivatives
				dR1(b,c) = N(1,b)*M(0,c)*weight;
				dsum[0]+=dR1(b,c);
				dR2(b,c)=N(0,b)*M(1,c)*weight;
				dsum[1]+=dR2(b,c);
			}
		}
		
		// divide by sum only required in terms of rational basis functions
		double inv_sum = 1.0/sum;  //(Breitenberger 17.06.2014 if condition removed)
		// divide through by sum
		for(int c=0;c<=m_q;c++)
		{
			for(int b=0;b<=m_p;b++)
			{
				R(b,c) = inv_sum*R(b,c);
				dR1(b,c) = inv_sum*dR1(b,c) - R(b,c)*dsum[0]*inv_sum; //  Breitenberger 18.09.2013
				dR2(b,c) = inv_sum*dR2(b,c) - R(b,c)*dsum[1]*inv_sum; //  Breitenberger 18.09.2013
			}
		}
		Matrix g_matrix(3,3);
		g_matrix.clear();
		
		for(int c=0;c<=m_q;c++)
		{
			for(int b=0;b<=m_p;b++)
			{

				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				weight=m_control_points[control_point_index].getWeight();			
		
				g_matrix(0,0) += dR1(b,c)*m_control_points[control_point_index].getX();
				g_matrix(1,0) += dR1(b,c)*m_control_points[control_point_index].getY();
				g_matrix(2,0) += dR1(b,c)*m_control_points[control_point_index].getZ();
				g_matrix(0,1) += dR2(b,c)*m_control_points[control_point_index].getX();
				g_matrix(1,1) += dR2(b,c)*m_control_points[control_point_index].getY();
				g_matrix(2,1) += dR2(b,c)*m_control_points[control_point_index].getZ();
			}
		}

		// Compute g3
		Vector g1_act = ZeroVector(3);
		g1_act(0) = g_matrix(0,0);
		g1_act(1) = g_matrix(1,0);
		g1_act(2) = g_matrix(2,0);
		Vector g2_act = ZeroVector(3);
		g2_act(0) = g_matrix(0,1);
		g2_act(1) = g_matrix(1,1);;
		g2_act(2) = g_matrix(2,1);;	
		Vector g3_act = MathUtils<double>::CrossProduct(g1_act,g2_act);
		g3_act /= norm_2(g3_act);
		g_matrix(0,2) = g3_act(0);
		g_matrix(1,2) = g3_act(1);
		g_matrix(2,2) = g3_act(2);	

		return g_matrix;
	}

	//  #####################################################################################
	// #######################################################################################
	//#
	///  \details    returns the basis vectors in the reference configuration
	///
	///  \return
	///  \param[out]  _g     basis vector 1 and basis vector 2 and the hessian matrix
	///
	///  \author     M.Breitenberger (04/2014)
	//
	//########################################################################################
	void ComputeBaseVectorsAndDerivatives( double _u, 
										   double _v, 
										   Vector& _g1, 
										   Vector& _g2, 
										   Vector& _g3, 
										   Matrix& _h )
	{
		_g1.resize(3);
		_g2.resize(3);
		_g3.resize(3);
		_h.resize(3,3);
		_g1.clear();
		_g2.clear();
		_g3.clear();
		_h.clear();

		int span_u = find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		int span_v = find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		Matrix dR; //1st Derivatives of shape functions in matrix format
		Matrix ddR; //2nd Derivatives of shape functions in vector format
		EvaluateNURBSFunctionsDerivatives(span_u,span_v, _u, _v, dR,ddR);

		int k = 0;
		for(int c=0;c<=m_q;c++)
		{
			for(int b=0;b<=m_p;b++)
			{

				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				_g1[0] += dR(k,0)*m_control_points[control_point_index].getX();
				_g2[0] += dR(k,1)*m_control_points[control_point_index].getX();
				
				_g1[1] += dR(k,0)*m_control_points[control_point_index].getY();
				_g2[1] += dR(k,1)*m_control_points[control_point_index].getY();
				
				_g1[2] += dR(k,0)*m_control_points[control_point_index].getZ();
				_g2[2] += dR(k,1)*m_control_points[control_point_index].getZ();

				_h(0,0) += ddR(k,0)*m_control_points[control_point_index].getX();
				_h(0,1) += ddR(k,1)*m_control_points[control_point_index].getX();
				_h(0,2) += ddR(k,2)*m_control_points[control_point_index].getX();

				_h(1,0) += ddR(k,0)*m_control_points[control_point_index].getY();
				_h(1,1) += ddR(k,1)*m_control_points[control_point_index].getY();
				_h(1,2) += ddR(k,2)*m_control_points[control_point_index].getY();

				_h(2,0) += ddR(k,0)*m_control_points[control_point_index].getZ();
				_h(2,1) += ddR(k,1)*m_control_points[control_point_index].getZ();
				_h(2,2) += ddR(k,2)*m_control_points[control_point_index].getZ();

				k++;

			}
		}

		//basis vector _g3
		_g3 = MathUtils<double>::CrossProduct(_g1,_g2);

		//differential area _dA
		double dA = norm_2(_g3);
		
		//normal vector _n
		_g3 = _g3/dA;
	}

	// --------------------------------------------------------------------------
	void FlagControlPointsForMapping(int span_u, int span_v, double _u, double _v)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		// Loop in the same order as for the evaluation of the NURBs functiosn
		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Flag control point as relevant for mapping
				m_control_points[control_point_index].SetRelevantForMapping();
			}
		}
	}

	// --------------------------------------------------------------------------
	void FlagControlPointsAsRelevantAndActive(int span_u, int span_v, double _u, double _v)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		// Loop in the same order as for the evaluation of the NURBs functiosn
		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Flag control point as relevant for mapping
				m_control_points[control_point_index].SetRelevantForMapping();
				m_control_points[control_point_index].SetActive();
			}
		}
	}

	// --------------------------------------------------------------------------
	matrix<unsigned int> GetMappingMatrixIds(int span_u, int span_v, double _u, double _v)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		// Initialize matrices
		matrix<unsigned int> mapping_matrix_ids = zero_matrix<unsigned int>( m_p+1 ,m_q+1 );

		// Loop in the same order as for the evaluation of the NURBs functiosn
		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Store control point Id in corresponding matrix
				mapping_matrix_ids(b,c) = m_control_points[control_point_index].GetMappingMatrixId();
			}
		}

		return mapping_matrix_ids;
	}

	// --------------------------------------------------------------------------
	matrix<int> GetRelevantControlPointsIndexes(int span_u, int span_v, double _u, double _v)
	{
		if(span_u==-1) span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		if(span_v==-1) span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		// Initialize matrices
		matrix<int> control_points_matrix = zero_matrix<int>( m_p+1 ,m_q+1 );

		// Loop in the same order as for the evaluation of the NURBs functiosn
		for (int c=0;c<=m_q;c++)
		{
			for (int b=0;b<=m_p;b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u-m_p+b;
				int vi = span_v-m_q+c;
				int control_point_index =vi*m_n_u + ui;

				// Store control point in corresponding matrix
				control_points_matrix(b,c) = control_point_index;
			}
		}

		return control_points_matrix;
	}

	// --------------------------------------------------------------------------
	void SetGrevilleAbscissae(int control_point_index, double& u_cp, double& v_cp)
	{
		// compute indeces of control point
			// control_point_index = v_index * m_n_u + u_index, hence
			int u_index = control_point_index % m_n_u;
			int v_index = control_point_index / m_n_u;

		// compute Greville abscissae
			// u_cp = (1/m_p) * (u_i + ... + u_i+p)
			// v_cp = (1/m_q) * (v_i + ... + v_i+p)	
			u_cp = 0;
			v_cp = 0;
		
			for(int i=u_index; i<=u_index+m_p; i++)
			{
				u_cp += m_knot_vector_u[i];
			}
			u_cp /= m_p;

			for(int i=v_index; i<=v_index+m_q; i++)
			{
				v_cp += m_knot_vector_v[i];
			}
			v_cp /= m_q;
	}

	// --------------------------------------------------------------------------
	IntVector GetKnotSpan(double _u, double _v)
	{
		int span_u=find_Knot_Span(m_knot_vector_u,_u,m_p,m_n_u);
		int span_v=find_Knot_Span(m_knot_vector_v,_v,m_q,m_n_v);

		IntVector span(2);
		span[0] = span_u;
		span[1] = span_v;

		return span;
	}
		
	// --------------------------------------------------------------------------
	DoubleVector& GetKnotVectorU()
	{
		return m_knot_vector_u;
	}

	// --------------------------------------------------------------------------
	DoubleVector& GetKnotVectorV()
	{
		return m_knot_vector_v;
	}

	// --------------------------------------------------------------------------
	ControlPointVector& GetControlPoints()
	{
		return m_control_points;
	}

	// --------------------------------------------------------------------------
	int& GetDegP()
	{
		return m_p;
	}

	// --------------------------------------------------------------------------
	int& GetDegQ()
	{
		return m_q;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "NURBSSurface";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "NURBSSurface";
	}

	// ==============================================================================
	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}


private:
	// ==============================================================================
	// Initialized by class constructor
	// ==============================================================================
	DoubleVector m_knot_vector_u;
	DoubleVector m_knot_vector_v;
	int m_p;
	int m_q;
	ControlPointVector m_control_points;
	unsigned int m_n_u; // number of control points in u-direction
	unsigned int m_n_v; // number of control points in v-direction

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      NURBSSurface& operator=[NURBSSurface const& rOther];

	/// Copy constructor.
	//      NURBSSurface[NURBSSurface const& rOther];

}; // Class NURBSSurface

} // namespace Kratos.

#endif // NURBS_SURFACE_H
