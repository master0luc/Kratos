//

#include <iostream>
#include<cmath>

#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
#include "utilities/math_utils.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "includes/ublas_interface.h"
#include "solid_mechanics_application.h"
namespace Kratos
{


NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule()
  :FlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:FlowRule(pYieldCriterion)
{
   
}


NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(const NonAssociativeExplicitPlasticFlowRule& rOther)
  :FlowRule(rOther)
{

}

FlowRule::Pointer NonAssociativeExplicitPlasticFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new NonAssociativeExplicitPlasticFlowRule(*this));
  return p_clone;
}


NonAssociativeExplicitPlasticFlowRule::~NonAssociativeExplicitPlasticFlowRule()
{
}


// THIS FUNCTION IS EQUIVALENT OF THE NEXT ONE  BUT IMPLICIT (RETURN MAPPING)
// TO BE PUT IN A DERIVED CLASS.
/*bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeformationGradientF0, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool PlasticityActive = false;
    double Tolerance = 1e-5;
//    std::cout <<  " HOLA " << 1.0-rDeltaDeformationGradient(0,0) << " 2 " << 1.0 - rDeltaDeformationGradient(1,1) << std::endl;
    InternalVariables PlasticVariables = mInternalVariables;
    PlasticVariables.EquivalentPlasticStrainOld = 0.0;

    rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

    rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
    rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen);

    Vector StressVector;
    this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);
  
    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, StressVector, PlasticVariables.EquivalentPlasticStrain);

    if ( ElasticTrialStateFunction > Tolerance)  {
        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
        rReturnMappingVariables.DeltaBeta = 0.0;
        rReturnMappingVariables.NormIsochoricStress = PlasticVariables.EquivalentPlasticStrain; 
        rReturnMappingVariables.EigenVectors = rNewElasticLeftCauchyGreen; 
        this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, StressVector, PlasticVariables.EquivalentPlasticStrain, ElasticTrialStateFunction, Tolerance);
    }

    rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain;
    if ( PlasticityActive) {
        rReturnMappingVariables.LameMu = PlasticVariables.EquivalentPlasticStrainOld; 
    }
    else {
        rReturnMappingVariables.LameMu = 0.0;
    }
    rStressMatrix = MathUtils<double>::StressVectorToTensor(StressVector);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

    return PlasticityActive;

}
*/

bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeformationGradientF0, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{

   bool PlasticityActive = false;
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

   //1- Initialize some variables
   Vector StressVector = ZeroVector(6);
   InternalVariables PlasticVariables = mInternalVariables;
   PlasticVariables.DeltaPlasticStrain = 0.0;
   
   Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen; // = ZeroMatrix(3);
   PlasticVariables.EquivalentPlasticStrainOld = 0.0;
   //PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
   
   Matrix NewElasticLeftCauchyGreen;

   double Tolerance = 1E-5;


   //2- Check for yield Condition (at the beggining)
   Vector NewStressVector;
   Vector PreviousStressVector;

   this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, PlasticVariables.EquivalentPlasticStrain);

    if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
            this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, PlasticVariables.EquivalentPlasticStrain, rReturnMappingVariables.TrialStateFunction, Tolerance);
            //THE INITIAL TENSIONAL STATE MAY BE OUTSIDE THE YIELD SURFACE??
    }

    //1. Compute an Elastic Trial State
    double StressErrorMeasure;
    double NewEquivalentPlasticStrain = 0.0;  // ONLY USED AS RETURN
    double NewPlasticShearStrain = 0.0;   // ONLY USED AS RETURN

    this->CalculateOneExplicitStep(rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, NewEquivalentPlasticStrain, NewPlasticShearStrain, false, StressErrorMeasure);

    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, NewEquivalentPlasticStrain);

    //2a. It is an elastic step
    if ( ElasticTrialStateFunction <= Tolerance)   {
        //std::cout << " COMPLETE ELASTIC " << std::endl;
        //Error does not matter because we are analitically in an hyperelastic bla bla bla
        PlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain; // NECESARY?? COPIED IN THE FUNCTION 
        PlasticityActive = false;
    }
    //2b. It is an elasto-plastic step
    else {
         if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
            //2b.1 Transition from elastic to plastic
	    //std::cout << " ELASTIC LOADING STEP " << std::endl;
            this->CalculateExplicitSolutionWithChange( rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, Tolerance);
          
            PlasticityActive = true;
         }
         else {         

            bool UnloadingCondition;
            UnloadingCondition =  this->EvaluateElastoPlasticUnloadingCondition( UnloadingCondition, PreviousElasticLeftCauchyGreen, rDeltaDeformationGradient, PlasticVariables, Tolerance); 
            if (UnloadingCondition ) {
                 //2c.1 ElastoPlastic Unloading
	         //std::cout << " ELASTIC UNLOADING STEP " << std::endl;
                 this->CalculateExplicitSolutionWithChange( rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, Tolerance);

                 PlasticityActive = true; 
            }
            else {
	         //std::cout << " PLASTIC STEP " << std::endl;
                //2c. 2 Completedly ELastoPlastic Step 
                this->CalculateExplicitSolution( rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);

                PlasticityActive = true;
            }
         }
 

    }


 //   PlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;


//    if (PlasticityActive) {
//
//       double DriftViolation;
//       DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
//
//       if ( fabs( DriftViolation) > Tolerance ) {
//            this->ReturnStressToYieldSurface( NewElasticHenckyStrain, NewStressVector, PlasticVariables.EquivalentPlasticStrain, DriftViolation, Tolerance);
//       }
//    } 

    // COMO EL TEMA CONMUTA forque f = eye(3), luego se puede aprovechar lo anterior!!! BIEN!!!!


    rReturnMappingVariables.NormIsochoricStress = PlasticVariables.EquivalentPlasticStrain;
    if (PlasticityActive) {
       double DriftViolation;
     
       DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
       rReturnMappingVariables.DeltaBeta = 0.0;
       rReturnMappingVariables.EigenVectors = NewElasticLeftCauchyGreen; 
       if ( fabs(DriftViolation) > Tolerance ) {
            //std::cout << " rDrift " << DriftViolation << std::endl;
            this->ReturnStressToYieldSurface( rReturnMappingVariables, NewElasticLeftCauchyGreen, NewStressVector, PlasticVariables.EquivalentPlasticStrain, DriftViolation, Tolerance);
       }


    }

    rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain;
    if ( PlasticityActive) {
         // SAVE THE DEVIATORIC INCREMENTAL PLASTIC DEFORMATION
        rReturnMappingVariables.LameMu = PlasticVariables.EquivalentPlasticStrainOld; 
    }
    else {
        rReturnMappingVariables.LameMu = 0.0;
    }

    if (PlasticityActive) {
        //rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain + 1.0;
    } 
    rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);
    rNewElasticLeftCauchyGreen = NewElasticLeftCauchyGreen; 

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);

   return PlasticityActive;

}


bool& NonAssociativeExplicitPlasticFlowRule::EvaluateElastoPlasticUnloadingCondition( bool& rUnloadingCondition, const Matrix& rElasticLeftCauchyGreen, const Matrix& rDeltaDeformationGradient, const InternalVariables& rPlasticVariables, const double& rTolerance)
{

     rUnloadingCondition = false;
    
     Vector StressVector;
     Vector YieldFunctionD;
     this->CalculateKirchhoffStressVector( rElasticLeftCauchyGreen, StressVector);
     mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, YieldFunctionD, rPlasticVariables.EquivalentPlasticStrain);

     Matrix ElasticMatrix;
     Vector ElasticStrainVector ;
     Vector IncrementalHencky;

     ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen);
     IncrementalHencky = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient) ));
     this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

     Vector IncrementalElasticStress;
     IncrementalElasticStress = prod(ElasticMatrix, IncrementalHencky);

     double Numerador = 0.0;
     double Denom1 = 0.0;
     double Denom2 = 0.0;

     for (unsigned int i = 0; i < 6; ++i) {
         Numerador += YieldFunctionD(i)*IncrementalElasticStress(i);
         Denom1 += pow( YieldFunctionD(i), 2.0);
         Denom2 += pow( IncrementalElasticStress(i), 2.0);
      }

      Denom1 = pow( Denom1*Denom2, 0.5);
    
      Numerador /= Denom1;

      if (Numerador < rTolerance) 
          rUnloadingCondition = true;

       return rUnloadingCondition;

}



// *************** YIELD SURFACE VIOLATION DRIFT CORRECTION **********
// *******************************************************************
            
void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rAlpha, double& rDrift, const double& rTolerance)
{

/*     Vector ActualElasticHenckyStrain;
     ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

     this->ReturnStressToYieldSurface( ActualElasticHenckyStrain, rStressVector, rAlpha, rDrift, rTolerance);

     rNewElasticLeftCauchyGreen = this->ConvertHenckyStrainToCauchyGreenTensor( ActualElasticHenckyStrain);
 */    // HAY ALGO MAL PORQUE DEJA DE CONVERGER LO BIEN QUE LO HACIA 
 
    AuxiliarDerivativesStructure  AuxiliarDerivatives;
    Matrix ElasticMatrix;
    double H;
    double DeltaGamma;
    Vector ElasticCorrection;
    Matrix  CorrectedLeftCauchyGreen;
    Matrix  UpdateMatrix;
    Vector ActualElasticHenckyStrain;
    Vector StressVector;
  
    rReturnMappingVariables.DeltaBeta = 0.0;
    for (unsigned int i = 0; i < 150; ++i) {

        ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

        this->UpdateDerivatives(ActualElasticHenckyStrain, AuxiliarDerivatives, rAlpha);
        this->ComputeElasticMatrix(ActualElasticHenckyStrain, ElasticMatrix);
        this->ComputePlasticHardeningParameter(ActualElasticHenckyStrain, rAlpha, H);

        DeltaGamma = rDrift;
        DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

        Vector AuxF = prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

        double denominador = 0.0;
        for (unsigned int k = 0; k < 6; ++k)
             denominador += AuxiliarDerivatives.YieldFunctionD(k) * AuxF(k);

        denominador += H;
        Vector MuyAuxiliar = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD/ 2.0;
        UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor( MuyAuxiliar);

 
        CorrectedLeftCauchyGreen =  prod((UpdateMatrix), rNewElasticLeftCauchyGreen);
        CorrectedLeftCauchyGreen =  prod( CorrectedLeftCauchyGreen, trans(UpdateMatrix));
    
        double AlphaUpdate = 0.0;
        for (unsigned int j = 0; j < 3; ++j)
            AlphaUpdate += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j);

        rAlpha += AlphaUpdate;
        rReturnMappingVariables.DeltaBeta += DeltaGamma ;
        this->CalculateKirchhoffStressVector( CorrectedLeftCauchyGreen, rStressVector);

//        std::cout << " DRIFT CORRECTION " << i << " " << rDrift << std::endl;
        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, rAlpha);

        rNewElasticLeftCauchyGreen = CorrectedLeftCauchyGreen;

//        this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);

//        std::cout << i << " Drift " <<  rDrift << " rAlpha " << rAlpha <<  "DELTA " << DeltaGamma  << std::endl;
//        std::cout << rNewElasticLeftCauchyGreen << std::endl;
//        std::cout << " stress " << StressVector << std::endl;
//        std::cout << " AlphPrev " << rAlpha-AlphaUpdate << " UPDATE " << AlphaUpdate << std::endl;
//        std::cout << " h  " << H <<  " " << MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)) << std::endl;
//        std::cout << " " << std::endl;
        if( fabs(rDrift) < rTolerance) {
//            std::cout << " DRIFT LEAVING" << rDrift << " " << rTolerance << std::endl;
            return;
        }

    }

    if ( fabs(rDrift) > 10.0*rTolerance) {
       std::cout<< " " << std::endl;
       std::cout << "Leaving drift correction WITHOUT converging " << rDrift << std::endl;
       std::cout << " StressVector " << rStressVector << std::endl;
}
/*     this->ReturnStressToYieldSurface( ActualElasticHenckyStrain, rStressVector, rAlpha, rDrift, rTolerance);

     rNewElasticLeftCauchyGreen = this->ConvertHenckyStrainToCauchyGreenTensor( ActualElasticHenckyStrain);
     // HAY ALGO MAL PORQUE DEJA DE CONVERJER LO BIEN QUE LO HACIA 
    
 */    

}


void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface(Vector& rElasticHenckyStrainVector, Vector& rStressVector, double& rAlpha, double& rDrift, const double& rTolerance)
{
    AuxiliarDerivativesStructure  AuxiliarDerivatives;
    Matrix ElasticMatrix;
    double H;
    double DeltaGamma;
    Vector ElasticCorrection;
        

    for (unsigned int i = 0; i < 15; ++i)     {


        this->UpdateDerivatives(rElasticHenckyStrainVector, AuxiliarDerivatives, rAlpha);
        this->ComputeElasticMatrix(rElasticHenckyStrainVector, ElasticMatrix);
        this->ComputePlasticHardeningParameter(rElasticHenckyStrainVector, rAlpha, H);

        DeltaGamma = rDrift;
        DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

        ElasticCorrection = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD;

        rElasticHenckyStrainVector += ElasticCorrection;

        for (unsigned int j = 0; j<3; ++j)
           rAlpha -= ElasticCorrection(j);

        this->CalculateKirchhoffStressVector( rElasticHenckyStrainVector, rStressVector);

        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, rAlpha);
    
        if (fabs(rDrift) < rTolerance)  {
            return;
        }
     }

     if (fabs(rDrift) > 10.0*rTolerance)
         std::cout << "Leaving Drift Correction whitout converging! " << rDrift << std::endl;


}

// ************** CONVERT HENCKY STRAIN TO CAUCHY-GREEN TENSOR  *********
// **********************************************************************

Vector NonAssociativeExplicitPlasticFlowRule::ConvertCauchyGreenTensorToHenckyStrain(const Matrix& rCauchyGreenMatrix)
{
    Matrix EigenVectors;
    Vector EigenValues;
    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues);
    
    Matrix Aux = ZeroMatrix(3);
    for (unsigned int i = 0; i < 3; ++i)
        Aux(i,i) = (std::log(EigenValues(i)))/2.0;

    Aux = prod(Aux, (EigenVectors));
    Aux = prod(trans(EigenVectors), Aux);
    Vector Result; 
    Result = MathUtils<double>::StrainTensorToVector(Aux, 6);
    return Result;
}

// ********************** AND THE CONTRARY ***************
// **99**88**77**66**
Matrix NonAssociativeExplicitPlasticFlowRule::ConvertHenckyStrainToCauchyGreenTensor(const Vector& rElasticHenckyStrain)
{

     Matrix HenckyStrainMatrix;
     HenckyStrainMatrix = MathUtils<double>::StrainVectorToTensor(rElasticHenckyStrain);

     Matrix EigenVectors;
     Vector EigenValues;
 
     SolidMechanicsMathUtilities<double>::EigenVectors(HenckyStrainMatrix, EigenVectors, EigenValues);

     Matrix Aux = ZeroMatrix(3);

    for (unsigned int i = 0; i < 3; ++i)
        Aux(i,i) = std::exp(2.0*EigenValues(i));

    Aux = prod(Aux, (EigenVectors));
    Aux = prod(trans(EigenVectors), Aux);

    // COMPROVACION
/*    Matrix A = ZeroMatrix(3);
    for (unsigned int i = 0; i < 3; ++ i ) 
        A(i,i) = EigenValues(i);

    A = prod(trans(EigenVectors), A);
    A = prod(A, EigenVectors);

    A = A - HenckyStrainMatrix; 
*/
    return Aux;

}


bool NonAssociativeExplicitPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables & rReturnMappingVariables)
{
   mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.DeltaGamma;
   mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.LameMu ; 
   return 0;
}


void NonAssociativeExplicitPlasticFlowRule::UpdateDerivatives(const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const double& rAlpha)
{
   Vector StressVector = ZeroVector(6);
   this->CalculateKirchhoffStressVector( rHenckyStrain, StressVector);
   this->CalculatePlasticPotentialDerivatives(StressVector, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);
   mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, rAuxiliarDerivatives.YieldFunctionD, rAlpha);
   rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;
}

//****** COMPUTE INCREMENTAL DEFORMATION GRADIENT  ********* 
// FROM BATHE NOTATION ^{rFinalConfiguration} _{rReferenceConfiguration} X
//*********************************************************
void NonAssociativeExplicitPlasticFlowRule::ComputeIncrementalDeformationGradient(const Matrix& rInitialDeformationGradient, const Matrix& rFinalDeformationGradient, const double & rReferenceConfiguration, const double & rFinalConfiguration, Matrix& rIncrementalDeformationGradient)
{

    Matrix DeformationGradientReference;
    Matrix DeformationGradientFinal;

    DeformationGradientReference = rReferenceConfiguration*rFinalDeformationGradient + ( 1.0 - rReferenceConfiguration) * rInitialDeformationGradient;
    DeformationGradientFinal     =     rFinalConfiguration*rFinalDeformationGradient + ( 1.0 -     rFinalConfiguration) * rInitialDeformationGradient;

    double Det; 
    MathUtils<double>::InvertMatrix(DeformationGradientReference, rIncrementalDeformationGradient, Det);
    rIncrementalDeformationGradient = prod( DeformationGradientFinal, rIncrementalDeformationGradient);
   
}




//********************* UPDATES ONE ELASTOPLASTIC STRAIN INCREMENT THE STRESS ****************
//******************************************************************************

void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitPlasticStep(const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, double& rNewEquivalentPlasticStrain, double& rNewPlasticShearStrain)
{
        rNewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;
        Matrix ElasticMatrix;
        Vector ElasticStrainVector ;

        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rPreviousElasticLeftCauchyGreen);
        this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

       	AuxiliarDerivativesStructure AuxiliarDerivatives;
       	this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rNewEquivalentPlasticStrain);

        double H;
        this->ComputePlasticHardeningParameter(ElasticStrainVector, rNewEquivalentPlasticStrain, H);


        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient)) );


        Vector auxVector;
        auxVector = prod(ElasticMatrix, ElasticStrainVector) ;
        double DeltaGamma = 0.0;
        for (unsigned int i = 0; i<6; ++i)
             DeltaGamma += auxVector(i)*AuxiliarDerivatives.YieldFunctionD(i);

        double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));

	DeltaGamma /= auxDenominador;

        if (DeltaGamma < 0 )
                DeltaGamma =0;

        Vector MuyAuxiliar;
        MuyAuxiliar = -DeltaGamma * AuxiliarDerivatives.PlasticPotentialD / 2.0; 

        Matrix UpdateMatrix;
        UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor (MuyAuxiliar);
        UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);

        rNewElasticLeftCauchyGreen = prod(UpdateMatrix, rPreviousElasticLeftCauchyGreen);
        rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

        double PlasticPotentialP = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
             PlasticPotentialP += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i);

        rNewEquivalentPlasticStrain +=  PlasticPotentialP ;// 3.0; 

        rNewPlasticShearStrain = 0;
        for (unsigned int i = 0; i < 3; ++i)
            rNewPlasticShearStrain += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i) - rNewEquivalentPlasticStrain, 2.0);

        for (unsigned int i = 3; i < 6; ++i)
            rNewPlasticShearStrain += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i), 2.0);

        rNewPlasticShearStrain = pow( rNewPlasticShearStrain, 1.0/2.0) ;
       

}

//********************* UPDATES ONE STRAIN INCREMENT THE STRESS ****************
//******************************************************************************

void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitStep(const Matrix& rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, double& rNewEquivalentPlasticStrain, double& rNewPlasticShearStrain, const bool & rElastoPlasticBool, double& rStressErrorMeasure)
{
        
  //double DeltaGamma = 0.0; not used JM

    if ( rElastoPlasticBool)  {


        // ONE STEP UPDATE
        this->CalculateOneExplicitPlasticStep( rDeltaDeformationGradient, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewEquivalentPlasticStrain, rNewPlasticShearStrain);

        this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector); 

        //TWO STEP UPDATE
        Matrix FinalDeformationGradient = prod(rDeltaDeformationGradient, rDeformationGradientF0);

        Matrix HalfStepIncrementalDeformation;
        ComputeIncrementalDeformationGradient(rDeformationGradientF0, FinalDeformationGradient, 0.0, 0.5, HalfStepIncrementalDeformation); 
        
        Matrix SecondApproxLeftCauchyGreen; 

        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, rPreviousElasticLeftCauchyGreen, rPlasticVariables, SecondApproxLeftCauchyGreen, rNewEquivalentPlasticStrain, rNewPlasticShearStrain);

        ComputeIncrementalDeformationGradient(rDeformationGradientF0, FinalDeformationGradient, 0.5, 1.0, HalfStepIncrementalDeformation);

        double PreviousEquivalent = rPlasticVariables.EquivalentPlasticStrain;
        double PreviousPlasticShearStrain = rNewPlasticShearStrain;
        rPlasticVariables.EquivalentPlasticStrain = rNewEquivalentPlasticStrain;

        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, SecondApproxLeftCauchyGreen, rPlasticVariables, SecondApproxLeftCauchyGreen, rNewEquivalentPlasticStrain, rNewPlasticShearStrain);

         rPlasticVariables.EquivalentPlasticStrain = PreviousEquivalent;
         rNewPlasticShearStrain += PreviousPlasticShearStrain;

         Vector SecondApproxStressVector;
         this->CalculateKirchhoffStressVector( SecondApproxLeftCauchyGreen, SecondApproxStressVector);

         rStressErrorMeasure = 0.0;
         double Denominador = 0.0;
         for (unsigned int i = 0; i<6; ++i) {
            rStressErrorMeasure += pow( rNewStressVector(i) - SecondApproxStressVector(i), 2.0);
            Denominador  += pow(rNewStressVector(i), 2.0);
         }
          rStressErrorMeasure = 1.0* pow( rStressErrorMeasure/Denominador, 0.5);

         rNewStressVector = SecondApproxStressVector;
         rNewElasticLeftCauchyGreen = SecondApproxLeftCauchyGreen;
         
   } 
   else  {
    
      rNewElasticLeftCauchyGreen = prod( rPreviousElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
      rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen );
  
      rStressErrorMeasure = 0.0;
      rNewPlasticShearStrain = 0.0;
      rNewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;   
      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);
     
   }


}


//***************** CalculateStress From LEFT CAUCHY GREEN MATRIX **************
void NonAssociativeExplicitPlasticFlowRule::CalculateKirchhoffStressVector( const Matrix& rElasticLeftCauchyGreen, Vector& rNewStressVector)
{

    Vector ElasticHenckyStrainVector;
    ElasticHenckyStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen); 
    this->CalculateKirchhoffStressVector(ElasticHenckyStrainVector, rNewStressVector);
}




void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolution( const Matrix & rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix&  rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance)
{
   double TimeStep = 0.5;
   double MinTimeStep = 1.0e-3;
   double DoneTimeStep = 0.0;
   double MaxTimeStep = 0.5;
   TimeStep = MaxTimeStep;

   Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);
  
   Matrix ActualDeformationGradient = rDeformationGradientF0;
   Matrix ActualElasticLeftCauchyGreen = rPreviousElasticLeftCauchyGreen;
   Matrix IncrementalDeformationGradientF;

   double StressErrorMeasure ;
   double NewEquivalentPlasticStrain;
   double NewPlasticShearStrain = 0.0;

   while (DoneTimeStep < 1.0)  {
      //std::cout << " TIME " << DoneTimeStep << " " << TimeStep  << std::endl;
       if (DoneTimeStep + TimeStep >= 1.0) {
           TimeStep = 1.0 - DoneTimeStep;
       }

       this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, DoneTimeStep, DoneTimeStep + TimeStep, IncrementalDeformationGradientF);
 
       this->CalculateOneExplicitStep( IncrementalDeformationGradientF, ActualDeformationGradient, ActualElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, NewPlasticShearStrain, rElastoPlasticBool, StressErrorMeasure);


       if ( StressErrorMeasure < rTolerance ) {
           // Se acepta el paso
           //std::cout << "   SUbstepping " << DoneTimeStep << " dt " << TimeStep << "Error " << StressErrorMeasure << std::endl;
           ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
           ActualDeformationGradient = prod( IncrementalDeformationGradientF, ActualDeformationGradient);
           DoneTimeStep += TimeStep;
           rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
           rPlasticVariables.EquivalentPlasticStrainOld += NewPlasticShearStrain;
       }
       else {
          if (TimeStep <= MinTimeStep) {
              if ( StressErrorMeasure > 10.0*rTolerance) {
                 std::cout << "ExplicitStressIntegrationDidNotConverged" << StressErrorMeasure << " " << TimeStep << std::endl;
                 std::cout << "    stress " << rNewStressVector << std::endl;
                 std::cout << "    PreviousEGCG" << rPreviousElasticLeftCauchyGreen << std::endl;
                 std::cout << "    Actual  EGCG" <<    ActualElasticLeftCauchyGreen << std::endl;
                 std::cout << "    New     EGCG" << rNewElasticLeftCauchyGreen << std::endl;
                 std::cout << "    rAlpha      " << NewEquivalentPlasticStrain << std::endl;
                 double Hola;
                 Hola = mpYieldCriterion->CalculateYieldCondition(Hola, rNewStressVector, NewEquivalentPlasticStrain);
                 std::cout << "    rYieldValue " << Hola << std::endl;
              }
              //RadialReturnVariables ReturnMappingVariables;
              //this->ReturnStressToYieldSurface( ReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, Hola, rTolerance);

              ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
              ActualDeformationGradient = prod( IncrementalDeformationGradientF, ActualDeformationGradient);
              DoneTimeStep += TimeStep;
              rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
              rPlasticVariables.EquivalentPlasticStrainOld += NewPlasticShearStrain;
          }
 
 
      }

      TimeStep *= 0.9* (pow( rTolerance / (StressErrorMeasure+ 1e-8), 1.0/2.0 ));
      TimeStep = std::max(TimeStep, MinTimeStep);
      TimeStep = std::min(TimeStep, MaxTimeStep);
//      TimeStep = std::min(TimeStep, 1e-1);
  }

}





// Search for the change and then continue
void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolutionWithChange(const Matrix& rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector,  const double& rTolerance)
{

   Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);

// bisecction in elastic Regime

   double InitialPosition = 0.0;
   double EndPosition = 1.0;
   double HalfPosition;

   double InitialStateFunction = -1.0;
   double EndStateFunction = 1.0;
   double HalfStateFunction;

   Matrix IncrementalDeformationGradientF;

   double StressErrorMeasure;
   double NewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain ;
   double NewPlasticShearStrain  = 0.0;
   double Measure1;
   double Measure2;

   for (unsigned int i = 0 ; i< 150; ++i) {

        HalfPosition = 0.5*(InitialPosition + EndPosition);
        this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, 0, HalfPosition, IncrementalDeformationGradientF);

        this->CalculateOneExplicitStep( IncrementalDeformationGradientF, rDeformationGradientF0, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, NewPlasticShearStrain, false, StressErrorMeasure);
        
        HalfStateFunction = mpYieldCriterion->CalculateYieldCondition( HalfStateFunction, rNewStressVector, rPlasticVariables.EquivalentPlasticStrain);

       if ( HalfStateFunction < 0)  {
           InitialStateFunction = HalfStateFunction;
           InitialPosition = HalfPosition;
       }
       else {
           EndPosition = HalfPosition;
           EndStateFunction = HalfStateFunction;
       }
       Measure1 = fabs(InitialStateFunction-EndStateFunction);
       Measure2 = fabs(InitialPosition - EndPosition);

       if ( (Measure1 < rTolerance ) && (Measure2 < rTolerance)) {
           break;
       }
 
   }
   //std::cout << " Half Position Found " << HalfPosition << " SF " << HalfStateFunction << std::endl;

   // COMPUTE ELASTIC STEP;

   this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, 0, HalfPosition, IncrementalDeformationGradientF);

   this->CalculateOneExplicitStep( IncrementalDeformationGradientF, rDeformationGradientF0, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, NewPlasticShearStrain, false, StressErrorMeasure);

   // COMPUTE ELASTOPLASTIC STEP
   Matrix HalfStepDeformationGradient = prod(IncrementalDeformationGradientF, rDeformationGradientF0);

   this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, HalfPosition, 1, IncrementalDeformationGradientF); 

   Matrix HalfStepLCG = rNewElasticLeftCauchyGreen;
   this->CalculateExplicitSolution( IncrementalDeformationGradientF, HalfStepDeformationGradient, HalfStepLCG, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  true, rTolerance);


}

   
void NonAssociativeExplicitPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix)
{

   Vector ElasticStrainVector = ZeroVector(6);
       
   ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain(rLeftCauchyGreenMatrix);

   this->ComputeElasticMatrix(ElasticStrainVector, rElasticMatrix);


   //if (rReturnMappingVariables.Control.PlasticRegion) {
   if (rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ) {


   // MATRIX 1.
      AuxiliarDerivativesStructure AuxiliarDerivatives;
      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rAlpha);

      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rAlpha, H);
 
      Vector AuxVectorF;
      Vector AuxVectorG;

      AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), rElasticMatrix);
      AuxVectorG = prod( rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

      Matrix PlasticUpdate;
      PlasticUpdate = MyCrossProduct(rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

      rElasticMatrix -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));

//      std::cout << " MATRIX 1 " << rElasticMatrix << std::endl;
      return;
      // TRY TO COMPUTE THE CONSISTENT TANGENT MATRIX.
      if ( rReturnMappingVariables.DeltaBeta == 0.0) {
         return;
      }
      // MATRIX 2

      ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain( rReturnMappingVariables.EigenVectors);

      Matrix ElasticMatrixTrial;
      this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrixTrial);
      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rReturnMappingVariables.NormIsochoricStress);

//      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rReturnMappingVariables.NormIsochoricStress, H);
 
      Matrix AuxiliarInverse = ElasticMatrixTrial;

      int Singular ;
      Singular = SolidMechanicsMathUtilities<double>::InvertMatrix( ElasticMatrixTrial, AuxiliarInverse);
  //    std::cout << Hola << " HOLA 1 : " << rElasticMatrix << std::endl;
      AuxiliarInverse += rReturnMappingVariables.DeltaBeta * AuxiliarDerivatives.PlasticPotentialDD;

//      std::cout<< " AUX INV " << AuxiliarInverse << " UPD " << rReturnMappingVariables.DeltaBeta * AuxiliarDerivatives.PlasticPotentialDD << std::endl;

      Singular = SolidMechanicsMathUtilities<double>::InvertMatrix( AuxiliarInverse,ElasticMatrixTrial); 
  //    std::cout << Hola << " HOLA 2 : " << rElasticMatrix << std::endl;

//      Vector AuxVectorF;
//      Vector AuxVectorG;

      AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), ElasticMatrixTrial);
      AuxVectorG = prod( ElasticMatrixTrial, AuxiliarDerivatives.PlasticPotentialD);

//      Matrix PlasticUpdate;
      PlasticUpdate = MyCrossProduct(ElasticMatrixTrial, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

      ElasticMatrixTrial -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));
//      rElasticMatrix = ElasticMatrixTrial; 
//      std::cout << " MATRIX 2 " << ElasticMatrixTrial << std::endl;
//      std::cout << "    DBeta " << rReturnMappingVariables.DeltaBeta << std::endl;
//      std::cout << std::endl;

  }

}


Matrix NonAssociativeExplicitPlasticFlowRule::MyCrossProduct(const Matrix& rM, const Vector& rA, const Vector& rB)
{
 

    Vector A = rA;
    Vector B = rB;

    for (unsigned int i = 3; i<6; ++i) {
//      A(i) /= 2.0;
//      B(i) /= 2.0;
    }

    A = prod(rM, A);
    B = prod(trans(B), rM); 

    Matrix Result = ZeroMatrix(6);
    for (unsigned int i = 0; i<6; ++i) {
        for ( unsigned int j = 0; j<6; ++j) {
            Result(i,j) = A(i)*B(j);
        }
     }
     return Result;

}

void NonAssociativeExplicitPlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void NonAssociativeExplicitPlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}



} //End Namepace Kratos
