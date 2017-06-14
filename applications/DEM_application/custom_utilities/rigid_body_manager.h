#ifndef KRATOS_RIGID_BODY_MANAGER_H_INCLUDED
#define KRATOS_RIGID_BODY_MANAGER_H_INCLUDED

#include <iostream>
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "GeometryFunctions.h"
#include "utilities/timer.h"
#include "custom_elements/rigid_body_element.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{

class RigidBodyManager {

    public:

    KRATOS_CLASS_POINTER_DEFINITION(RigidBodyManager);

    /// Default constructor

    RigidBodyManager() {}

    /// Destructor

    virtual ~RigidBodyManager();

    /// Turn back information as a string
    virtual std::string Info() const;

    /// Print information about this object
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data
    virtual void PrintData(std::ostream& rOStream) const;

    protected:

        vector<RigidBodyElement3D> mVectorOfRigidBodies;

    private:

    /// Assignment operator
    RigidBodyManager& operator=(RigidBodyManager const& rOther);

    }; // Class RigidBodyManager

} // namespace Kratos

#endif // KRATOS_RIGID_BODY_MANAGER_H_INCLUDED
