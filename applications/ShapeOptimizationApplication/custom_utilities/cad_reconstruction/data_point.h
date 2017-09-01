// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Giovanni Filomeno
                     giovanni.filomeno@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
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
//   Created by:          $Author:    giovanni.filomeno@tum.de   $
//   Last modified by:    $Co-Author: giovanni.filomeno@tum.de   $
//   Date:                $Date:                      Decem 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef DATA_POINT_H_
#define DATA_POINT_H_

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
// #include <boost/python.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "patch.h"
// ==============================================================================

namespace Kratos
{
class DataPoint
{
public:

    // ==========================================================================
    // Type definitions
    // ==========================================================================
    typedef std::vector<int> IntVector;
	typedef Node<3> NodeType;

    /// Pointer definition of DataPoint
    KRATOS_CLASS_POINTER_DEFINITION(DataPoint);

    /// Constructor 1
    DataPoint(NodeType::Pointer node_ptr)
    {
        m_node_ptr = node_ptr;
    }
    
    /// Constructor 2
    DataPoint(double X, double Y, double Z)
    {
        m_node_ptr = NodeType::Pointer(new NodeType(-1, X, Y, Z)); // node id is set to -1
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_X, 0);
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_Y, 0);
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_Z, 0);
    }

    /// Constructor 3
    DataPoint(double X, double Y, double Z,
              double updated_X, double updated_Y, double updated_Z)
    {
        m_node_ptr = NodeType::Pointer(new NodeType(-1, X, Y, Z)); // node id is set to -1
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_X, updated_X);
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_Y, updated_Y);
        m_node_ptr->SetValue(SHAPE_CHANGE_ABSOLUTE_Z, updated_Z);
    }

    // Getters
        NodeType::Pointer getNodePtr()
        {
            return m_node_ptr;
        }

        double getX()
        {
            return m_node_ptr->X() + m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_X);
        }

        double getY()
        {
            return m_node_ptr->Y() + m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_Y);
        }

        double getZ()
        {
            return m_node_ptr->Z() + m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_Z);
        }

        double getU()
        {
            return m_u;
        }

        double getV()
        {
            return m_v;
        }

    void setPatch(Patch& patch)
    {
        m_patch = patch;
    }

    // void setU(double u)
    // {
    //     m_u = u;
    // }

    // void setV(double v)
    // {
    //     m_v = v;
    // }

    void updateUAndV(double u, double v)
    {
        m_u = u;
        m_v = v;
        m_patch.GetSurface().EvaluateSurfacePoint(m_cad_point, m_u, m_v);
    }

    Vector getDistanceToOriginalDataPointVector()
    {
        Vector distance_vector = ZeroVector(3);
        distance_vector(0) = m_cad_point.X() - m_node_ptr->X();
        distance_vector(1) = m_cad_point.Y() - m_node_ptr->Y();
        distance_vector(2) = m_cad_point.Z() - m_node_ptr->Z();
        return distance_vector;
    }

    void optimize_parametrisation(double tol, unsigned int max_itr)
    {
        Matrix hessian = ZeroMatrix(2,2);
        Vector gradient = ZeroVector(2);
        double determinant_of_hessian = 0;
		Matrix inverse_of_hessian = ZeroMatrix(2,2);        

        double norm_deltau = 1;
        unsigned int k = 0;
        while (norm_deltau > tol)
        {
            // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
            Vector Q_minus_P = getDistanceToOriginalDataPointVector();

            // The distance is used to compute Hessian and gradient
            m_patch.GetSurface().EvaluateGradientsForClosestPointSearch(Q_minus_P,
                                                                        hessian,
                                                                        gradient,
                                                                        m_u, m_v);

            // u_k and v_k are updated
            MathUtils<double>::InvertMatrix( hessian, inverse_of_hessian, determinant_of_hessian);
            Vector deltau = prod(inverse_of_hessian, gradient);
            updateUAndV(m_u - deltau(0), m_v - deltau(1)); // Q is updated

            norm_deltau = norm_2(deltau);
            k++;

            if(k > max_itr)
            {
                std::cout << "WARNING!!! Newton-Raphson to find closest point did not converge in the following number of iterations: " << k-1 << std::endl;
                return;
            }
        }
    }

    //
        // IntVector getKnotSpan()
        // {
        //     return m_patch.GetSurface().GetKnotSpan(m_u, m_v);
        // }

        // double getX0()
        // {
        //     return m_node_ptr->X();
        // }

        // double getY0()
        // {
        //     return m_node_ptr->Y();
        // }

        // double getZ0()
        // {
        //     return m_node_ptr->Z();
        // }

        // double getdX()
        // {
        //     return m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_X);
        // }

        // double getdY()
        // {
        //     return m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_Y);
        // }

        // double getdZ()
        // {
        //     return m_node_ptr->GetValue(SHAPE_CHANGE_ABSOLUTE_Z);
        // }  

    /// Destructor.
    virtual ~DataPoint()
    {
    }
    // ==============================================================================


//    evaluat

    // ==============================================================================


        /// Turn back information as a string.
        virtual std::string Info() const
        {
            return "DataPoint";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "DataPoint";
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
        }


private:
    // ==============================================================================
    // General working arrays
    // ==============================================================================
    // FE data
    NodeType::Pointer m_node_ptr;
    // CAD point data
    Patch m_patch;
    double m_u, m_v;
    Point<3> m_cad_point; 

}; // Class DataPoint

}  // namespace Kratos.

#endif // DATA_POINT_H_
