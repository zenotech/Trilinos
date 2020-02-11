// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBEObserverDefault_hpp
#define Tempus_StepperBEObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBEObserverBase.hpp"


namespace Tempus {

/** \brief Default observer for StepperBE.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperBEObserverBase for details on the algorithm.
 */
template<class Scalar>
class StepperBEObserverDefault
  : virtual public Tempus::StepperBEObserverBase<Scalar>
{
public:

  /// Constructor
  StepperBEObserverDefault(){}

  /// Destructor
  virtual ~StepperBEObserverDefault(){}

  /// Observe BE Stepper at end of takeStep.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperBE<Scalar> > /* stepper */,
    const typename StepperBEAppAction<Scalar>::ACTION_LOCATION actLoc) const
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

#endif // Tempus_StepperBEObserverDefault_hpp
