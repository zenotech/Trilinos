// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBEAppActionComposite_hpp
#define Tempus_StepperBEAppActionComposite_hpp

#include "Tempus_StepperBEAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template<class Scalar>
class StepperBEAppActionComposite
  : virtual public Tempus::StepperBEAppAction<Scalar>
{
public:

  /// Default constructor
  StepperBEAppActionComposite();

  /// Destructor
  virtual ~StepperBEAppActionComposite();

  /// Execute application action for BE Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperBE<Scalar> > stepper,
    const typename StepperBEAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for(auto& a : appActions_)
      a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addBEAppAction(Teuchos::RCP<StepperBEAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearBEAppActions();
  { appActions_.clear();}

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

private:

  std::vector<Teuchos::RCP<StepperBEAppAction<Scalar > > > appActions_;

};

} // namespace Tempus
#endif // Tempus_StepperBEAppActionComposite_hpp
