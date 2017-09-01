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
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef BOUNDARY_LOOP_H
#define BOUNDARY_LOOP_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "boundary_edge.h"

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

class BoundaryLoop
{
public:
	///@name Type Definitions
	///@{

	typedef std::vector<BoundaryEdge> BoundaryEdgeVector;

	///@}

	/// Pointer definition of BoundaryLoop
	//    KRATOS_CLASS_POINTER_DEFINITION[BoundaryLoop];

	/// Default constructor.
	BoundaryLoop(BoundaryEdgeVector boundary_edges, bool is_inner_loop=false) 
	: m_boundary_edges(boundary_edges),
	  m_is_inner_loop(is_inner_loop)
	{
		// Create polygon for inside / outside check
		CreatePolygon();
	}

	/// Destructor.
	virtual ~BoundaryLoop()
	{
	}

	// --------------------------------------------------------------------------
	void CreatePolygon()	
	{
		std::cout << "\n Better creation of points on polygon needed, like the one in Carat e.g. Currently only hard coded!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		unsigned int np_per_edge = 500;

		// Variables needed
		unsigned int counter=0;

		// Loop over all boundary edges creating the closed boundary loop
		for (BoundaryEdgeVector::iterator edge_i = m_boundary_edges.begin(); edge_i != m_boundary_edges.end(); ++edge_i)
		{
			DoubleVector& knot_vec_u = edge_i->GetKnotVectorU();

			double u_min = knot_vec_u[0];
			double u_max = knot_vec_u[knot_vec_u.size()-1];
			double delta_u = (u_max-u_min) / np_per_edge;

			// Resize to incude points on current edge
			m_boundary_polygon.resize(m_boundary_polygon.size()+np_per_edge);

			// Add points of edge to polygon
			double u_i = u_min;
			for(unsigned int i=0;i<np_per_edge;i++)
			{
				u_i += delta_u;
				Point<3> curve_point;

				edge_i->EvaluateCurvePoint(curve_point,u_i);

				m_boundary_polygon[counter][0]=curve_point.X();
				m_boundary_polygon[counter][1]=curve_point.Y();

				counter++;
			}			
		}
	}

	// --------------------------------------------------------------------------
	BoundaryEdgeVector& GetBoundaryEdges()
	{
		return m_boundary_edges;
	}

	// --------------------------------------------------------------------------
	std::vector<array_1d<double, 2>>& GetBoundaryPolygon()
	{
		return m_boundary_polygon;
	}	

	// --------------------------------------------------------------------------
	bool IsInnerLoop()
	{
		return m_is_inner_loop;
	}
	

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "BoundaryLoop";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "BoundaryLoop";
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
	BoundaryEdgeVector m_boundary_edges;
	bool m_is_inner_loop;
	std::vector<array_1d<double, 2>> m_boundary_polygon;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      BoundaryLoop& operator=[BoundaryLoop const& rOther];

	/// Copy constructor.
	//      BoundaryLoop[BoundaryLoop const& rOther];

}; // Class BoundaryLoop

} // namespace Kratos.

#endif // BOUNDARY_LOOP_H
