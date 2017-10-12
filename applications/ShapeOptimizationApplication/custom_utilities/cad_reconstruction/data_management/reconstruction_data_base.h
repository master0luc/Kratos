// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_DATA_BASE_H
#define RECONSTRUCTION_DATA_BASE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../basic_nurbs_brep_handling/patch.h"
#include "../basic_nurbs_brep_handling/brep_element.h"
#include "cad_model_reader.h"

// ==============================================================================

namespace Kratos
{
class ReconstructionDataBase
{
public:

    typedef boost::python::extract<int> ExtractInt;    

    /// Pointer definition of ReconstructionDataBase
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionDataBase);

    /// Default constructor.
    ReconstructionDataBase(ModelPart& fe_model_part, boost::python::dict cad_geometry, boost::python::dict cad_integration_data)
	: mpFEModelPart(fe_model_part),
      mrCADGeometry(cad_geometry),
      mrCADIntegrationData(cad_integration_data)
    {     
    }

    /// Destructor.
    virtual ~ReconstructionDataBase()
    {
    }

    // ------------------------------------------------------------------------------
    void Create()
    {
        CADModelReader cad_reader(mrCADGeometry,mrCADIntegrationData);
        
        cad_reader.ReadGeometry(mPatches);
            if(len(mrCADIntegrationData.keys())>0)

        // Create map to identify patch in patch vector for a given patch_id
        for(auto & patch_i : mPatches)
        {
            unsigned int index_in_patch_vector = &patch_i - &mPatches[0];
            mPatchVectorIndexFromPatchId[patch_i.GetId()] = index_in_patch_vector;
        }   

        cad_reader.ReadIntegrationData(mBrepElements);      
    }

    // --------------------------------------------------------------------------
    void AddToControlPointDisplacements( Vector& delta_displacements )
    { 
        // Map to get control point corresponding to a given global id (needed for python update later)
        std::map<unsigned int, ControlPoint*> control_point_corresponding_to_global_id;
        
        // 1. Step: Update C++ data base
        for(auto & patch_i : mPatches) 
        {
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                if(control_point_i.IsRelevantForReconstruction())
                {
                    // Updating c++ data base
                    unsigned int cp_equation_id = control_point_i.GetEquationId();
                    control_point_i.SetdX( control_point_i.GetdX() + delta_displacements[3*cp_equation_id+0] );
                    control_point_i.SetdY( control_point_i.GetdY() + delta_displacements[3*cp_equation_id+1] );
                    control_point_i.SetdZ( control_point_i.GetdZ() + delta_displacements[3*cp_equation_id+2] );

                }
                // Filling map to be used later
                unsigned int cp_global_id = control_point_i.GetGlobalId();
                control_point_corresponding_to_global_id[cp_global_id] = &control_point_i;
            }
        }

        // 2. Step: Update pyhon data base
        for (int i = 0; i < boost::python::len(mrCADGeometry["faces"]); i++)
        {
            for (int cp_idx = 0; cp_idx < boost::python::len(mrCADGeometry["faces"][i]["surface"][0]["control_points"]); cp_idx++)
            {

                unsigned int global_id = ExtractInt(mrCADGeometry["faces"][i]["surface"][0]["control_points"][cp_idx][0]);
                ControlPoint* cp_j = control_point_corresponding_to_global_id[global_id];

                double position_x = cp_j->GetX();
                double position_y = cp_j->GetY();
                double position_z = cp_j->GetZ();

                // Update python data base which is stored as a reference (so it is also updated on python level)
                mrCADGeometry["faces"][i]["surface"][0]["control_points"][cp_idx][1][0] = position_x;
                mrCADGeometry["faces"][i]["surface"][0]["control_points"][cp_idx][1][1] = position_y;
                mrCADGeometry["faces"][i]["surface"][0]["control_points"][cp_idx][1][2] = position_z;
            }
        }                         
    }

    // ------------------------------------------------------------------------------
    ModelPart& GetFEModelPart()
    {
        return mpFEModelPart;   
    }    

    // ------------------------------------------------------------------------------
    std::vector<Patch>& GetPatchVector()
    {
        return mPatches;   
    }
    
    // ------------------------------------------------------------------------------    
    Patch& GetPatchFromPatchId( unsigned int patch_id )
    {
        return mPatches[mPatchVectorIndexFromPatchId[patch_id]];
    }

    // ------------------------------------------------------------------------------
    std::vector<BREPElement>& GetBREPElements()
    {
        return mBrepElements;   
    } 

    // ==============================================================================


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ReconstructionDataBase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ReconstructionDataBase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


private:
    
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mpFEModelPart;
    boost::python::dict mrCADGeometry;
	boost::python::dict mrCADIntegrationData;
    std::vector<Patch> mPatches;
    std::vector<BREPElement> mBrepElements;
    std::map<unsigned int, unsigned int> mPatchVectorIndexFromPatchId;

}; // Class ReconstructionDataBase

}  // namespace Kratos.

#endif // RECONSTRUCTION_DATA_BASE_H
