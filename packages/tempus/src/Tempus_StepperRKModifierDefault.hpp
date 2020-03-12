// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKModifierDefault_hpp
#define Tempus_StepperRKModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKModifierBase.hpp"


namespace Tempus {

/** \brief Default modifier for StepperRK.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperRKModifierBase for details on the algorithm.
 */
template<class Scalar>
class StepperRKModifierDefault
  : virtual public Tempus::StepperRKModifierBase<Scalar>
{
public:

  /// Constructor
  StepperRKModifierDefault(){}

  /// Destructor
  virtual ~StepperRKModifierDefault(){}

  /// Modify RK Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperERK<Scalar> > /* stepper */,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP:
      case StepperRKAppAction<Scalar>::BEGIN_STAGE:
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE:
      case StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperRKAppAction<Scalar>::END_STAGE:
      case StepperRKAppAction<Scalar>::END_STEP:
      {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

};

} // namespace Tempus

#endif // Tempus_StepperRKModifierDefault_hpp
