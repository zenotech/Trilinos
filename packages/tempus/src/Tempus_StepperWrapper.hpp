// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperWrapper_hpp
#define Tempus_StepperWrapper_hpp

// Tempus
#include "Tempus_Stepper.hpp"


namespace Tempus {


/** \brief A generic wrapper to Stepper.
 *
 *  <b>Design Considerations</b>
 *    - A Stepper that just wraps some code into a Stepper.
 *    - It has default implementations of everything, except takeStep(),
 *      which would contain the wrapped code.
 *    - This allows quick and easy creation of a Stepper that does not
 *      integrate the solution, x, but does other related work.
 */
template<class Scalar>
class StepperWrapper : virtual public Stepper
{
public:

  void warningMessage(const std::string func)
  {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperWrapper");
    *out << "Warning -- " << func << " is unimplemented in StepperWrapper.\n"
  }

  /// \name Basic stepper methods
  //@{

    /// Only required to implement this.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& /* appModel */)
    { warningMessage(__func__); }

    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& /* appModel */)
    { warningMessage(__func__); }

    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
    { warningMessage(__func__); return Teuchos::null; }

    /// Set solver via ParameterList solver name.
    virtual void setSolver(std::string /* solverName */)
    { warningMessage(__func__); }

    /// Set solver via solver ParameterList.
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> /* solverPL */=Teuchos::null)
    { warningMessage(__func__); }

    /// Set solver.
    virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > /* solver */)
    { warningMessage(__func__); }

    /// Get solver
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver()
    { warningMessage(__func__); return Teuchos::null; }

    /// Set Observer
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > /* obs */ = Teuchos::null)
    { warningMessage(__func__); }

    /// Initialize during construction and after changing input parameters.
    virtual void initialize()
    { warningMessage(__func__); }

    /// Set initial conditions, make them consistent, and set stepper memory.
    virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */)
    { warningMessage(__func__); }

    /// Pass initial guess to Newton solver (for implicit schemes)
    virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> >
      /* initial_guess */ = Teuchos::null)
    { warningMessage(__func__); }

    virtual std::string getStepperType() const
    { return "StepperWrapper"; }

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState()
    { warningMessage(__func__); }

    virtual Scalar getOrder()    const { return Scalar(0.0); }
    virtual Scalar getOrderMin() const { return Scalar(0.0); }
    virtual Scalar getOrderMax() const { return Scalar(0.0); }

    virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
    { warningMessage(__func__); return Scalar(0.0); }

    virtual Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const
    { warningMessage(__func__); return Teuchos::null; }

    virtual bool isExplicit() const { return false; }
    virtual bool isImplicit() const { return false; }
    virtual bool isExplicitImplicit() const
    { return isExplicit() and isImplicit(); }

    virtual bool isOneStepMethod() const   { return false; }
    virtual bool isMultiStepMethod() const { return false; }

    virtual OrderODE getOrderODE() const { return OrderODE::NONE; }

    virtual void setUseFSAL(bool /* a */) { warningMessage(__func__); }
    virtual bool getUseFSAL() const
    { warningMessage(__func__); return false; }

    virtual void setICConsistency(std::string /* s */)
    { warningMessage(__func__); }
    virtual std::string getICConsistency() const
    { warningMessage(__func__); return "None"; }

    virtual void setICConsistencyCheck(bool /* c */)
    { warningMessage(__func__); }
    virtual bool getICConsistencyCheck() const
    { warningMessage(__func__); return false; }
  //@}

};
} // namespace Tempus
#endif // Tempus_StepperWrapper_hpp
