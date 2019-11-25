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


/** \brief Default Modifier illustrating usage.
 *
 *  The Modifier allows applications to modify the solution and
 *  stage solutions in StepperExplicitRK::takeStep_modify() at
 *  specific locations in the ExplicitRK algorithm.  This is
 *  specified through the MODIFIER_TYPE.  See
 *  StepperExplicitRK<Scalar>::takeStep_modify() for details on
 *  algorithm.
 */
template<class Scalar>
class StepperExplicitRKModifierDefault
 : virtual public Tempus::StepperExplicitRKModifier<Scalar>
{
public:

  /// Constructor
  StepperExplicitRKModifierDefault(){}

  /// Destructor
  virtual ~StepperExplicitRKModifierDefault(){}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
                      MODIFIER_TYPE /* modType */){}
};

} // namespace Tempus


#endif // Tempus_StepperExplicitRKModifierDefault_hpp
