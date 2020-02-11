// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBEModifierX_hpp
#define Tempus_StepperBEModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBEModifierXBase.hpp"


namespace Tempus {

/** \brief Default ModifierX for StepperBE.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperBEModifierXBase for details on the algorithm.
 */
template<class Scalar>
class StepperBEModifierXDefault
  : virtual public Tempus::StepperBEModifierXBase<Scalar>
{
public:

  /// Constructor
  StepperBEModifierXDefault(){}

  /// Destructor
  virtual ~StepperBEModifierXDefault(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
    const Scalar /* time */, const Scalar /* dt */,
    const typename StepperBEModifierXBase<Scalar>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperBEModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperBEModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperBEModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperBEModifierXBase<Scalar>::XDOT_END_STEP:
      {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown modifier type.\n");
    }
  }

};

} // namespace Tempus

#endif // Tempus_StepperBEModifierX_hpp
