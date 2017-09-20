//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_DAM_APPLY_FORCE_BY_SPATIAL_POSITION_PROCESS )
#define  KRATOS_DAM_APPLY_FORCE_BY_SPATIAL_POSITION_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamApplyForceBySpatialPositionProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamApplyForceBySpatialPositionProcess);

    typedef Table<double,double> TableType;  
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamApplyForceBySpatialPositionProcess(ModelPart& rModelPart,
                                Parameters& rParameters
                                ) : Process(Flags()) , mrModelPart(rModelPart)
    {
        KRATOS_TRY
			 
        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name"   :"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id"           : 0,
                "variable_name"     : "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "value"             : 0.0,
                "table"             : 0,
                "position_x"        : 0.0,
                "position_y"        : 0.0,
                "position_z"        : 0.0
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["value"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mMeshId = rParameters["mesh_id"].GetInt();
        mVariableName = rParameters["variable_name"].GetString();
        mValue = rParameters["value"].GetDouble();

        // Getting the values of the device coordinates
        mForcePosition.resize(3,false);
        mForcePosition[0] = rParameters["position_x"].GetDouble();
        mForcePosition[1] = rParameters["position_y"].GetDouble();
        mForcePosition[2] = rParameters["position_z"].GetDouble();
        
        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();
        
        if(mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamApplyForceBySpatialPositionProcess() {}
  

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize()
    {
        
        KRATOS_TRY;

        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type var_component = KratosComponents< component_type >::Get(mVariableName);
        const int nelements = mrModelPart.GetMesh(mMeshId).Elements().size();
        bool IsInside = false;
        array_1d<double,3> LocalCoordinates;       
        Element::Pointer pSelectedElement;

        // Getting the values of table in case that it exist        
        if(mTableId != 0 )
        { 
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time/mTimeUnitConverter;
            mValue = mpTable->GetValue(time);
        }

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();
            int PointsNumber = 0;

            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                pSelectedElement = (*(it.base()));
                IsInside = pSelectedElement->GetGeometry().IsInside(mForcePosition,LocalCoordinates);
                if(IsInside) break; 
            }
            if(IsInside==false)
            {
                for(int k = 0; k<nelements; k++)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    pSelectedElement = (*(it.base()));
                    IsInside = pSelectedElement->GetGeometry().IsInside(mForcePosition,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
            {
                KRATOS_ERROR << "ERROR!!, PLEASE REPEAT THE SEARCH " << std::endl;
            }

            PointsNumber = pSelectedElement->GetGeometry().PointsNumber();
            std::vector<std::size_t> ConditionNodeIds(1);

            // The global condition index should be a ProcessInfo or sth similar
            for(int j = 0; j < PointsNumber; j++)
            {
                    pSelectedElement->GetGeometry().GetPoint(j).FastGetSolutionStepValue(var_component) = mValue;
                    ConditionNodeIds[0]= pSelectedElement->GetGeometry().GetPoint(j).Id();
                    mrModelPart.CreateNewCondition("PointLoadCondition2D1N", j+1, ConditionNodeIds, 0);
            }
            
            
        }
    
        KRATOS_CATCH("");
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamApplyForceBySpatialPositionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamApplyForceBySpatialPositionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    std::size_t mMeshId;
    std::string mVariableName;
    double mValue;
    array_1d<double,3> mForcePosition;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId; 
    

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamApplyForceBySpatialPositionProcess& operator=(DamApplyForceBySpatialPositionProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                        DamApplyForceBySpatialPositionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamApplyForceBySpatialPositionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_APPLY_FORCE_BY_SPATIAL_POSITION_PROCESS defined */

