// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKModifierBasic_hpp
#define Tempus_StepperExplicitRKModifierBasic_hpp


namespace Tempus {

/// MODIFIER_TYPE indicates solution variable and stepper location to modify.
enum MODIFIER_TYPE {
  STAGEX_BEGINSTAGE, ///< Modify stage x value at beginning of stage.
  X_ENDSTEP          ///< Modify solution x value at end step.
};


/** \brief Basic Modifier illustrating usage.
 *
 *  The Modifier allows applications to modify the solution and
 *  stage solutions in StepperExplicitRK::takeStep_modify() at
 *  specific locations in the ExplicitRK algorithm.  This is
 *  specified through the MODIFIER_TYPE.  See
 *  StepperExplicitRK<Scalar>::takeStep_modify() for details on
 *  algorithm.
 */
template<class Scalar>
class StepperExplicitRKModifierBasic
{
public:

  /// Constructor
  StepperExplicitRKModifierBasic(){}

  /// Destructor
  virtual ~StepperExplicitRKModifierBasic(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(RCP<Thyra::VectorBase<Scalar> > /* x */,
                      MODIFIER_TYPE /* modType */){}
};

} // namespace Tempus


#endif // Tempus_StepperExplicitRKModifierBasic_hpp
