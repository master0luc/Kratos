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

#ifndef PATCH_H
#define PATCH_H

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
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "control_point.h"
#include "nurbs_surface.h"
#include "boundary_loop.h"

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

class Patch
{
public:
	///@name Type Definitions
	///@{

	// For matrix / vector operations
	typedef std::vector<ControlPoint> ControlPointVector;
	typedef std::vector<BoundaryLoop> BoundaryLoopVector;

	///@}

	/// Pointer definition of Patch
	//    KRATOS_CLASS_POINTER_DEFINITION[Patch];

	/// Default constructor.
	Patch()
	{
	}

	Patch(unsigned int patch_id, NURBSSurface nurbs_surface, BoundaryLoopVector boundary_loops ) 
	: m_patch_id(patch_id),
	  m_nurbs_surface(nurbs_surface),
	  m_boundary_loops(boundary_loops)
	{
	}
	
	/// Destructor.
	virtual ~Patch()
	{
	}

	// --------------------------------------------------------------------------
	bool CheckIfPointIsInside(array_1d<double, 2> &point_of_interest)
	{
		// Boost is used to check whether point of interest is inside given polygon or not

		// Type definitions to use boost functionalities
		typedef boost::geometry::model::d2::point_xy<double> point_type;
		typedef boost::geometry::model::polygon<point_type> polygon_type;

		// We assume point is inside, check all boundary loops if this is true. If this is not true for a single loop, then point is considered outside of this patch
		point_type poi(point_of_interest[0],point_of_interest[1]);
		bool is_inside = true;

		// Loop over all boundary loops of current patch
		for(BoundaryLoopVector::iterator loop_i =  m_boundary_loops.begin(); loop_i!=m_boundary_loops.end(); loop_i++)
		{

			// Initialize necessary variables
			polygon_type poly;
			std::vector<array_1d<double, 2>>& boundary_polygon = loop_i->GetBoundaryPolygon();

			// Prepare polygon for boost
			for(unsigned int i=0;i<boundary_polygon.size();i++)
				boost::geometry::append( boost::geometry::exterior_ring(poly), 
										boost::geometry::make<point_type>(boundary_polygon[i][0], boundary_polygon[i][1]));
			if(boundary_polygon.size()>0)
				boost::geometry::append( boost::geometry::exterior_ring(poly), 
										boost::geometry::make<point_type>(boundary_polygon[0][0], boundary_polygon[0][1]));

			// Check inside or outside
			is_inside = boost::geometry::within(poi, poly);
			
			// Boost does not consider the polygon direction, it always assumes inside as in the interior of the closed polygon.
			// If the CAD loop is an inner loop, however, the area which is considered inner is the unbounded side of the closed polygon.
			// So we toggle the results from the within search to be correct.
			if(loop_i->IsInnerLoop())
				is_inside = !is_inside;

			// If a point is considered outside in one loop, it is outside in general, hence we can break the loop
			if(!is_inside) 
				break;
		}
		

		return is_inside;
	}
	
	bool CheckIfPointIsInside(double& u, double& v) 
	{	
		// Boost is used to check whether point of interest is inside given polygon or not

		// Type definitions to use boost functionalities
		typedef boost::geometry::model::d2::point_xy<double> point_type;
		typedef boost::geometry::model::polygon<point_type> polygon_type;

		point_type poi(u, v);

		// same as: bool CheckIfPointIsInside(array_1d<double, 2> &point_of_interest)
			// We assume point is inside, check all boundary loops if this is true. If this is not true for a single loop, then point is considered outside of this patch
			bool is_inside = true;

			// Loop over all boundary loops of current patch
			for(BoundaryLoopVector::iterator loop_i =  m_boundary_loops.begin(); loop_i!=m_boundary_loops.end(); loop_i++)
			{

				// Initialize necessary variables
				polygon_type poly;
				std::vector<array_1d<double, 2>>& boundary_polygon = loop_i->GetBoundaryPolygon();

				// Prepare polygon for boost
				for(unsigned int i=0;i<boundary_polygon.size();i++)
					boost::geometry::append( boost::geometry::exterior_ring(poly), 
											boost::geometry::make<point_type>(boundary_polygon[i][0], boundary_polygon[i][1]));
				if(boundary_polygon.size()>0)
					boost::geometry::append( boost::geometry::exterior_ring(poly), 
											boost::geometry::make<point_type>(boundary_polygon[0][0], boundary_polygon[0][1]));

				// Check inside or outside
				is_inside = boost::geometry::within(poi, poly);
				
				// Boost does not consider the polygon direction, it always assumes inside as in the interior of the closed polygon.
				// If the CAD loop is an inner loop, however, the area which is considered inner is the unbounded side of the closed polygon.
				// So we toggle the results from the within search to be correct.
				if(loop_i->IsInnerLoop())
					is_inside = !is_inside;

				// If a point is considered outside in one loop, it is outside in general, hence we can break the loop
				if(!is_inside) 
					break;
			}

		return is_inside;			
	}

	// --------------------------------------------------------------------------
	unsigned int GetId()
	{
		return m_patch_id;
	}

	// --------------------------------------------------------------------------
	NURBSSurface& GetSurface()
	{
		return m_nurbs_surface;
	}

	// --------------------------------------------------------------------------
	BoundaryLoopVector& GetBoundaryLoops()
	{
		return m_boundary_loops;
	}

	// ==============================================================================
	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "Patch";
	}

	// ==============================================================================
	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "Patch";
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
	unsigned int m_patch_id;
	NURBSSurface m_nurbs_surface;
	BoundaryLoopVector m_boundary_loops;
	std::vector<array_1d<double, 2>> m_closed_boundary_polygon;

	// ==============================================================================
	// General working arrays
	// ==============================================================================
	/// Assignment operator.
	//      Patch& operator=[Patch const& rOther];

	/// Copy constructor.
	//      Patch[Patch const& rOther];

}; // Class Patch

} // namespace Kratos.

#endif // PATCH_H
