// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBEModifierDefault_hpp
#define Tempus_StepperBEModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBEModifierBase.hpp"


namespace Tempus {

/** \brief Default modifier for StepperBE.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperBEModifierBase for details on the algorithm.
 */
template<class Scalar>
class StepperBEModifierDefault
  : virtual public Tempus::StepperBEModifierBase<Scalar>
{
public:

  /// Constructor
  StepperBEModifierDefault(){}

  /// Destructor
  virtual ~StepperBEModifierDefault(){}

  /// Modify BE Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperBE<Scalar> > /* stepper */,
    const typename StepperBEAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperBEAppAction<Scalar>::BEGIN_STEP:
      case StepperBEAppAction<Scalar>::BEFORE_SOLVE:
      case StepperBEAppAction<Scalar>::AFTER_SOLVE:
      case StepperBEAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperBEModifierDefault_hpp
