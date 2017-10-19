// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_PROJECTION_UTILITY_H
#define CAD_PROJECTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "reconstruction_data_base.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"

// ==============================================================================

namespace Kratos
{
class CADProjectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::python::extract<int> ExtractInt;
    typedef Node<3> NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;    
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;    
    typedef std::vector<Patch> PatchVector; 
    typedef std::vector<double> DoubleVector;
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of CADProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(CADProjectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADProjectionUtility( PatchVector& patch_vector, int max_projection_iterations, double projection_tolerance )
    : mrPatchVector( patch_vector ),
      mMaxIterations( max_projection_iterations ),
      mTolerance( projection_tolerance )
    {      
    }

    /// Destructor.
    virtual ~CADProjectionUtility()
    {
    }

    // --------------------------------------------------------------------------
    void Initialize( boost::python::list ParameterResolutionForCoarseNeighborSearch )
    {
        std::cout << "\n> Initializing CAD projection..." << std::endl;           
        boost::timer timer;

        CreateCADPointCloudUsingLists( ParameterResolutionForCoarseNeighborSearch );
        CreateSearchTreeWithCADPointCloud();

		std::cout << "> Time needed initializing CAD projection: " << timer.elapsed() << " s" << std::endl;    
    }
    
    // --------------------------------------------------------------------------
    void CreateCADPointCloudUsingLists( boost::python::list ParameterResolutionForCoarseNeighborSearch )
    {
      int u_resolution = ExtractInt( ParameterResolutionForCoarseNeighborSearch[0] );
      int v_resolution = ExtractInt( ParameterResolutionForCoarseNeighborSearch[1] );

      for (auto & patch_i : mrPatchVector)
      {
            int index_in_patch_vector = &patch_i - &mrPatchVector[0];
            DoubleVector& knot_vec_u_i = patch_i.GetSurfaceKnotVectorU();
            DoubleVector& knot_vec_v_i = patch_i.GetSurfaceKnotVectorV();
            std::cout << "> Processing Patch with brep_id " << patch_i.GetId() << std::endl;
      
            double u_min = knot_vec_u_i[0];
            double u_max = knot_vec_u_i[knot_vec_u_i.size()-1];
            double v_min = knot_vec_v_i[0];
            double v_max = knot_vec_v_i[knot_vec_v_i.size()-1];
            double delta_u = (u_max-u_min) / u_resolution;
            double delta_v = (v_max-v_min) / v_resolution;

            // Loop over all u & v according to specified resolution
            array_1d<double,2> point_in_parameter_space;
            for(unsigned int i=0; i<=u_resolution; i++)
            {
                point_in_parameter_space[0] = u_min + i*delta_u;

                for(unsigned int j=0; j<=v_resolution; j++)
                {
                    point_in_parameter_space[1] = v_min + j*delta_v;

                    bool point_is_inside = patch_i.IsPointInside(point_in_parameter_space);
                    if(point_is_inside)
                    {
                        ++mNumberOfNodesInCADPointCloud;					
                        Point<3> cad_point_coordinates;
                        patch_i.EvaluateSurfacePoint( point_in_parameter_space, cad_point_coordinates );

                         NodeType::Pointer new_cad_node = Node <3>::Pointer(new Node<3>(mNumberOfNodesInCADPointCloud, cad_point_coordinates));

                        mOrderedListOfNodes.push_back(new_cad_node);
                        mOrderedListOfParameterValues.push_back(point_in_parameter_space);
                        mOrderedListOfPatchIndices.push_back(index_in_patch_vector);
                    }
                }
            }
        }      
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithCADPointCloud()
    {
      mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mOrderedListOfNodes.begin(), mOrderedListOfNodes.end(), mBucketSize));
    }     

    // --------------------------------------------------------------------------
    void DetermineNearestCADPoint( NodeType::Pointer PointOfInterest,
                                   array_1d<double,2>& parameter_values_of_nearest_point,
                                   int& patch_index_of_nearest_point )
    {
        // 1) Coarse search in the point cloud
        NodeType::Pointer nearest_point = mpSearchTree->SearchNearestPoint( *PointOfInterest );

        // 2) Detailed projection using Newton-Raphson
        
        // Recover CAD information from lists representing point cloud
        parameter_values_of_nearest_point = mOrderedListOfParameterValues[nearest_point->Id()-1];
        patch_index_of_nearest_point = mOrderedListOfPatchIndices[nearest_point->Id()-1];
        Patch& patch_of_nearest_point = mrPatchVector[patch_index_of_nearest_point];
        
				// Initialize what's needed in the Newton-Raphson iteration				
				Vector Distance = ZeroVector(3); 
				Matrix hessian = ZeroMatrix(2,2);
				Vector gradient = ZeroVector(2);
				double determinant_of_hessian = 0;
        Matrix inverse_of_hessian = ZeroMatrix(2,2);
        Point<3> current_nearest_point;
        current_nearest_point[0] = nearest_point->X();
        current_nearest_point[1] = nearest_point->Y();
        current_nearest_point[2] = nearest_point->Z();

        // Variables neeed by the Netwon Raphson algorithm
				double norm_delta_u = 100000000;
        
        // Newton-Raphson algorithm if iterations are specified
        for(unsigned int k=0; k<mMaxIterations; k++)
        {
          // KRATOS_WATCH(k);
          // KRATOS_WATCH(current_nearest_point);
          // KRATOS_WATCH(*PointOfInterest);
          // KRATOS_WATCH(parameter_values_of_nearest_point);
          // The distance between point on CAD surface point on the FE-mesh
          Distance(0) = current_nearest_point[0] - PointOfInterest->X();
          Distance(1) = current_nearest_point[1] - PointOfInterest->Y();
          Distance(2) = current_nearest_point[2] - PointOfInterest->Z();
          
          // The distance is used to compute hessian and gradient
          patch_of_nearest_point.EvaluateGradientsForClosestPointSearch( Distance, hessian, gradient , parameter_values_of_nearest_point );

          // u_k and v_k are updated
          MathUtils<double>::InvertMatrix( hessian, inverse_of_hessian, determinant_of_hessian );
          Vector delta_u = prod(inverse_of_hessian,gradient);
          parameter_values_of_nearest_point[0] -= delta_u(0);
          parameter_values_of_nearest_point[1] -= delta_u(1);

          // Point on CAD surface is udpated
          patch_of_nearest_point.EvaluateSurfacePoint( parameter_values_of_nearest_point, current_nearest_point );
          
          // Check convergence
          norm_delta_u = norm_2(delta_u);
          if(norm_delta_u<mTolerance)
            break;
          else if(k+1==mMaxIterations)
          {
            std::cout << "WARNING!!! Newton-Raphson in projection did not converge in the following number of iterations: " << k+1 << std::endl;
            KRATOS_WATCH(current_nearest_point)
            KRATOS_WATCH(*PointOfInterest)
          }
        }
    }
    
    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADProjectionUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADProjectionUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

  private:

    
    // ==============================================================================
    // Variables initialized by constructor
    // ==============================================================================
    PatchVector& mrPatchVector;
    int mMaxIterations;
    double mTolerance;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mNumberOfNodesInCADPointCloud = 0;
    NodeVector mOrderedListOfNodes;
    std::vector<array_1d<double,2>> mOrderedListOfParameterValues; 
    std::vector<int> mOrderedListOfPatchIndices;
    unsigned int mBucketSize = 100;
    KDTree::Pointer mpSearchTree;

    /// Assignment operator.
    //      CADProjectionUtility& operator=(CADProjectionUtility const& rOther);

    /// Copy constructor.
    //      CADProjectionUtility(CADProjectionUtility const& rOther);

}; // Class CADProjectionUtility
} // namespace Kratos.

#endif // CAD_PROJECTION_UTILITY_H
