// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/point_3d.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
namespace Kratos
{
KratosStructuralMechanicsApplication::KratosStructuralMechanicsApplication():
    /* ELEMENTS */
    // Adding the beam element
    mSmallDisplacementBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    // Adding the shells elements
    mIsotropicShellElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mShellThickElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), false ),
    mShellThickCorotationalElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), true ),
    mShellThinElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), false ),
    mShellThinCorotationalElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), true ),
    // Adding the membrane element
    mMembraneElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    // Adding the SPRISM element
    mSprismElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    // Adding the nodal concentrated element
    mNodalConcentratedElement2D1N( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1 ) ) ) ),
    mNodalConcentratedElement3D1N( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1 ) ) ) ),
    /* CONDITIONS */
    // Beam's point moment condition
    mPointMomentCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    // Torque's point condition
    mPointTorqueCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ), //TOTAL LAGRANGIAN ELEMENTS
    mTotalLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mTotalLagrangian2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mTotalLagrangian2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mTotalLagrangian2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mTotalLagrangian3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mTotalLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mTotalLagrangian3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mTotalLagrangian3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mTotalLagrangian3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mTotalLagrangian3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),

    mKinematicLinear2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ), //TOTAL LAGRANGIAN ELEMENTS
    mKinematicLinear2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mKinematicLinear2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mKinematicLinear2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mKinematicLinear2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mKinematicLinear3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mKinematicLinear3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mKinematicLinear3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mKinematicLinear3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mKinematicLinear3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mKinematicLinear3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mKinematicLinear3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),
    mLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSurfaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSurfaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSurfaceLoadCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSurfaceLoadCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSurfaceLoadCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9 ) ) ) )


{}

void KratosStructuralMechanicsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "     KRATOS   ___|  |                   |                   |               " << std::endl;
    std::cout << "            \\___ \\  __|  __| |   |  __| __| |   |  __| _` | |               " << std::endl;
    std::cout << "                  | |   |    |   | (    |   |   | |   (   | |               " << std::endl;
    std::cout << "            _____/ \\__|_|   \\__,_|\\___|\\__|\\__,_|_|  \\__,_|_| MECHANICS     " << std::endl;


    // Generalized eigenvalue problem
    KRATOS_REGISTER_VARIABLE( BUILD_LEVEL )
    KRATOS_REGISTER_VARIABLE( EIGENVALUE_VECTOR )
    KRATOS_REGISTER_VARIABLE( EIGENVECTOR_MATRIX )

    // Geometrical
    KRATOS_REGISTER_VARIABLE( AREA )
    KRATOS_REGISTER_VARIABLE( IX )
    KRATOS_REGISTER_VARIABLE( IY )
    KRATOS_REGISTER_VARIABLE( IZ )
    KRATOS_REGISTER_VARIABLE( CROSS_AREA )
    KRATOS_REGISTER_VARIABLE( MEAN_RADIUS )
    KRATOS_REGISTER_VARIABLE( SECTION_SIDES )
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS )

    //  Shell generalized variables
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT_GLOBAL )

    // Cross section
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

    // Nodal stiffness
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_STIFFNESS )

    // CONDITIONS
    /* Moment condition */
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )
    /* Torque condition */
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_TORQUE )

//    // Orthotropy
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )

//    // Material orientation
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DX )
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DY )
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DZ )

    // Adding the SPRISM EAS variables
    KRATOS_REGISTER_VARIABLE(ALPHA_EAS);
    KRATOS_REGISTER_VARIABLE(EAS_IMP);
    KRATOS_REGISTER_VARIABLE(SPRISM_TL_UL);

    // Adding the SPRISM additional variables
    KRATOS_REGISTER_VARIABLE(ANG_ROT);

    // Adding the SPRISM number of transversal integration points
    KRATOS_REGISTER_VARIABLE(NINT_TRANS);

    // Adding the SPRISM variable to deactivate the quadratic interpolation
    KRATOS_REGISTER_VARIABLE(QUAD_ON);

    // Strain measures
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_VECTOR);
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_TENSOR);
    
    //material orientation
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DX )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DY )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DZ )
    
    //othotropic/anisotropic constants
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )

  KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS )
  KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN )
  KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_TENSOR )
  KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR )
  KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_VECTOR )
  KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS ) 
  
  KRATOS_REGISTER_VARIABLE( RAYLEIGH_ALPHA )
  KRATOS_REGISTER_VARIABLE( RAYLEIGH_BETA )
  
  
  //nodal load variables
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(  POINT_LOAD )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(  LINE_LOAD )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(  SURFACE_LOAD )
    
  //condition load variables
  KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR )
  KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR )
  KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR )
  KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURES_VECTOR )
  KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURES_VECTOR )


    // Register the beam element
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBeamElement3D2N", mSmallDisplacementBeamElement3D2N )

    //Register the shells elements
    KRATOS_REGISTER_ELEMENT( "IsotropicShellElement3D3N", mIsotropicShellElement3D3N )
    KRATOS_REGISTER_ELEMENT( "ShellThickElement3D4N", mShellThickElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThickElementCorotational3D4N", mShellThickCorotationalElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElement3D3N", mShellThinElement3D3N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElementCorotational3D3N", mShellThinCorotationalElement3D3N )
    
    //VOLUME ELEMENTS
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D3N", mTotalLagrangian2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D4N", mTotalLagrangian2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D6N", mTotalLagrangian2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D8N", mTotalLagrangian2D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D9N", mTotalLagrangian2D9N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D4N", mTotalLagrangian3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D6N", mTotalLagrangian3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D8N", mTotalLagrangian3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D10N", mTotalLagrangian3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D15N", mTotalLagrangian3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D20N", mTotalLagrangian3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D27N", mTotalLagrangian3D27N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D3N", mKinematicLinear2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D4N", mKinematicLinear2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D6N", mKinematicLinear2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D8N", mKinematicLinear2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D9N", mKinematicLinear2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D4N", mKinematicLinear3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D6N", mKinematicLinear3D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D8N", mKinematicLinear3D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D10N", mKinematicLinear3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D15N", mKinematicLinear3D15N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D20N", mKinematicLinear3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D27N", mKinematicLinear3D27N )
    
    // Register the membrane element
    KRATOS_REGISTER_ELEMENT( "MembraneElement3D3N", mMembraneElement3D3N )

    // Register the SPRISM element
    KRATOS_REGISTER_ELEMENT("SprismElement3D6N", mSprismElement3D6N);
    
    // Register the nodal concentrated element
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement2D1N", mNodalConcentratedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement3D1N", mNodalConcentratedElement3D1N);

    // Register the conditions
    // Beam's point moment condition
    KRATOS_REGISTER_CONDITION( "PointMomentCondition3D1N", mPointMomentCondition3D1N );
    // Torque moment condition
    KRATOS_REGISTER_CONDITION( "PointTorqueCondition3D1N", mPointTorqueCondition3D1N );
    
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D2N", mLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D3N", mLineLoadCondition2D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N )
    
}

}  // namespace Kratos.
