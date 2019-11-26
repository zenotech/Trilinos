// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKModifier_hpp
#define Tempus_StepperExplicitRKModifier_hpp


namespace Tempus {



/** \brief Pure virtual Modifier class.
 *
 *  The Modifier allows applications to modify the solution and
 *  stage solutions in StepperExplicitRK::takeStep_modify() at
 *  specific locations in the ExplicitRK algorithm.  The locations
 *  are specified through the MODIFIER_TYPE.  See
 *  StepperExplicitRK<Scalar>::takeStep_modify() for details on
 *  algorithm.
 */
template<class Scalar>
class StepperExplicitRKModifier
{
public:

  /// MODIFIER_TYPE indicates solution variable and stepper location to modify.
  enum MODIFIER_TYPE {
    X_BEGINSTEP,     ///< Modify x at the beginning of the step.
    X_BEGINSTAGE,    ///< Modify x at the beginning of the stage.
    XDOT_ENDSTAGE,   ///< Modify xDot at the end of the stage.
    X_ENDSTEP        ///< Modify x at the end of the step.
  };

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > x,
                      Scalar time, Scalar dt, MODIFIER_TYPE modType) = 0;
};

} // namespace Tempus


#endif // Tempus_StepperExplicitRKModifier_hpp
