// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKModifierDefault_hpp
#define Tempus_StepperExplicitRKModifierDefault_hpp

#include "Tempus_StepperExplicitRKModifier.hpp"


namespace Tempus {


/** \brief Default Modifier which is a no-op.
 *
 *  See StepperExplicitRKModifier for details.
 */
template<class Scalar>
class StepperExplicitRKModifierDefault
 : virtual public StepperExplicitRKModifier<Scalar>
{
public:

  /// Constructor
  StepperExplicitRKModifierDefault(){}

  /// Destructor
  virtual ~StepperExplicitRKModifierDefault(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
    Scalar /* time */, Scalar /* dt */,
    typename StepperExplicitRKModifier<Scalar>::MODIFIER_TYPE /* modType */){}
};

} // namespace Tempus


#endif // Tempus_StepperExplicitRKModifierDefault_hpp
