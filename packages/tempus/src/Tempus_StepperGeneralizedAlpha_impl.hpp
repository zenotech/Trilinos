// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperGeneralizedAlpha_impl_hpp
#define Tempus_StepperGeneralizedAlpha_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;

template<class Scalar>
StepperGeneralizedAlpha<Scalar>::StepperGeneralizedAlpha(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validSecondOrderODE_DAE(appModel);
  Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
    Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(appModel,
                                                      "Generalized-Alpha"));
  this->wrapperModel_ = wrapperModel;
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperGeneralizedAlpha::initialize()\n");

  xDotDotScratch_ = Thyra::createMember(this->wrapperModel_->get_x_space());
  assign(xDotDotScratch_.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

  this->setSolver();
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperGeneralizedAlpha::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperGeneralizedAlpha<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for HHTAlpha.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    RCP<const Thyra::VectorBase<Scalar> > xOld     = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > xDotOld  = currentState->getXDot();
    RCP<const Thyra::VectorBase<Scalar> > xDotDotOld=currentState->getXDotDot();

    RCP<Thyra::VectorBase<Scalar> > x       = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot    = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > xDotDot = workingState->getXDotDot();

    const Scalar rho_inf = this->stepperPL_->template get<double>("Rho_Inf");
    const Scalar alpha_m = this->stepperPL_->template get<double>("Alpha_m");
    const Scalar alpha_f = this->stepperPL_->template get<double>("Alpha_f");
    const Scalar beta    = this->stepperPL_->template get<double>("Beta"   );
    const Scalar gamma   = this->stepperPL_->template get<double>("Gamma"  );

    // Setup initial guess of x^{n+1-alpha_f_}_0, using initial x^{n+1}_0
    Thyra::V_StVpStV(x.ptr(), alpha_f, *xOld, (1.0-alpha_f), *x);

    // Get time and dt
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperGeneralizedAlphaTimeDerivative<Scalar>(
        dt, xOld, xDotOld, xDotDotOld, xDotDotScratch_,
        rho_inf, alpha_m, alpha_f, beta, gamma));

    // Setup InArgs and OutArgs
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
    MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
    inArgs.set_x        (x);
    inArgs.set_x_dot    (xDot);
    inArgs.set_x_dot_dot(xDotDot);

    const Scalar Omega = 1.0/((1.0-alpha_f)*beta*dt*dt);  // d(xDotDot)/dx
    const Scalar Alpha = gamma/(dt*beta);                 // d(xDot)/dx
    const Scalar Beta  = 1.0;                             // d(x)/dx

    if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
    if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
    if (inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
      inArgs.set_W_x_dot_dot_coeff(Omega);
    if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (Alpha);
    if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (Beta);

    this->wrapperModel_->setForSolve(timeDer, inArgs, outArgs);

    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(x);

    // Get solution at t^{n+1}
    Thyra::V_StVpStV(x.ptr(),       1.0/(1.0-alpha_f), *x,
                               -alpha_f/(1.0-alpha_f), *xOld);
    Thyra::V_StVpStV(xDot.ptr(),    1.0/(1.0-alpha_f), *xDot,
                               -alpha_f/(1.0-alpha_f), *xDotOld);
    Thyra::V_StVpStV(xDotDot.ptr(), 1.0/(1.0-alpha_m), *xDotDot,
                               -alpha_m/(1.0-alpha_m), *xDotDotOld);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
    workingState->setOrder(this->getOrder());
  }
  return;
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperGeneralizedAlpha<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperGeneralizedAlpha<Scalar>::description() const
{
  std::string name = "Generalized-Alpha";
  return(name);
}


template<class Scalar>
void StepperGeneralizedAlpha<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperGeneralizedAlpha<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (this->stepperPL_ == Teuchos::null)
      this->stepperPL_ = this->getDefaultParameters();
  } else {
    this->stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters.
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab ostab(out,1,"Generalized-Alpha");

  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  std::string stepperType = stepperPL->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Generalized-Alpha",
    std::logic_error,
       "\nError - Stepper Type is not 'Generalized-Alpha'!\n"
       << "Stepper Type = "
       << stepperPL->get<std::string>("Stepper Type") << "\n");

  auto method = stepperPL->get<std::string>("Method","Generalized-Alpha");
  double rho_inf_ = 0.0;
  double alpha_m_ = 0.0;
  double alpha_f_ = 0.0;
  double beta_    = 0.0;
  double gamma_   = 0.0;

  if (method == "Generalized-Alpha") {
    rho_inf_ = stepperPL->get<double>("Rho_Inf", double(1.0/3.0));
    TEUCHOS_TEST_FOR_EXCEPTION( (rho_inf_ > 1.0) || (rho_inf_ < 0.0),
      std::logic_error,
         "\nError in 'Generalized-Alpha' stepper: invalid value of Rho_Inf = "
         << rho_inf_ << ".  Please select 0 <= Rho_Inf <= 1. \n");

    // Values based on rho_inf_
    alpha_m_ = (2*rho_inf_ - 1.0)/(rho_inf_ + 1.0);
    alpha_f_ = rho_inf_/(rho_inf_ + 1.0);
    // Alpha values overridden by user
    alpha_m_ = stepperPL->get<double>("Alpha_m", alpha_m_);
    alpha_f_ = stepperPL->get<double>("Alpha_f", alpha_f_);
    // Values based on alphas
    beta_    = 0.25*pow((1.0 - alpha_m_ + alpha_f_), 2);
    gamma_   = 0.5 - alpha_m_ + alpha_f_;
    // Values overridden by user
    beta_    = stepperPL->get<double>("Beta", beta_);
    gamma_   = stepperPL->get<double>("Gamma", gamma_);

  } else if (method == "WBZ-Alpha") {
    rho_inf_ = stepperPL->get<double>("Rho_Inf", double(1.0/3.0));
    TEUCHOS_TEST_FOR_EXCEPTION( (rho_inf_ > 1.0) || (rho_inf_ < 0.0),
      std::logic_error,
         "\nError in 'WBZ-Alpha' stepper: invalid value of Rho_Inf = "
         << rho_inf_ << ".  Please select 0 <= Rho_Inf <= 1. \n");

    // Values based on rho_inf_
    alpha_m_ = (rho_inf_ - 1.0)/(rho_inf_ + 1.0);
    alpha_f_ = 0.0;
    // Alpha values overridden by user
    alpha_m_ = stepperPL->get<double>("Alpha_m", alpha_m_);
    if (stepperPL->isParameter("Alpha_f"))
      *out << "Warning - WBZ-Alpha ignoring Alpha_f! Set to zero." << std::endl;
    // Values based on alpha_m_
    beta_    = 0.25*pow((1.0 - alpha_m_), 2);
    gamma_   = 0.5 - alpha_m_;
    // Values overridden by user
    beta_    = stepperPL->get<double>("Beta", beta_);
    gamma_   = stepperPL->get<double>("Gamma", gamma_);

  } else if (method == "HHT-Alpha") {
    rho_inf_ = stepperPL->get<double>("Rho_Inf", double(0.5));
    TEUCHOS_TEST_FOR_EXCEPTION( (rho_inf_ > 1.0) || (rho_inf_ < 0.5),
      std::logic_error,
         "\nError in 'HHT-Alpha' stepper: invalid value of Rho_Inf = "
         << rho_inf_ << ".  Please select 0.5 <= Rho_Inf <= 1. \n");

    // Values based on rho_inf_
    alpha_m_ = 0.0;
    alpha_f_ = (1.0 - rho_inf_)/(rho_inf_ + 1.0);
    // Alpha values overridden by user
    if (stepperPL->isParameter("Alpha_m"))
      *out << "Warning - HHT-Alpha ignoring Alpha_m! Set to zero." << std::endl;
    alpha_f_ = stepperPL->get<double>("Alpha_f", alpha_f_);
    // Values based on alpha_f_
    beta_    = 0.25*pow((1.0 + alpha_f_), 2);
    gamma_   = 0.5 + alpha_f_;
    // Values overridden by user
    beta_    = stepperPL->get<double>("Beta", beta_);
    gamma_   = stepperPL->get<double>("Gamma", gamma_);

  } else if (method == "Newmark-Beta") {
    // Default values
    if (stepperPL->isParameter("Rho_Inf"))
      *out<<"Warning - Newmark-Beta ignoring Rho_inf! Set to zero."<<std::endl;
    rho_inf_ = 0.0;
    if (stepperPL->isParameter("Alpha_m"))
      *out<<"Warning - Newmark-Beta ignoring Alpha_m! Set to zero."<<std::endl;
    alpha_m_ = 0.0;
    if (stepperPL->isParameter("Alpha_f"))
      *out<<"Warning - Newmark-Beta ignoring Alpha_f! Set to zero."<<std::endl;
    alpha_f_ = 0.0;
    beta_    = stepperPL->get<double>("Beta", 0.25);
    gamma_   = stepperPL->get<double>("Gamma", 0.5);

  } else if (method == "Newmark Beta Average Acceleration") {
    // Default values
    if (stepperPL->isParameter("Rho_Inf"))
      *out << "Warning - Newmark Beta Average Acceleration ignoring Rho_inf! "
           << "Set to zero." << std::endl;
    rho_inf_ = 0.0;
    if (stepperPL->isParameter("Alpha_m"))
      *out << "Warning - Newmark Beta Average Acceleration ignoring Alpha_m! "
           << "Set to zero." << std::endl;
    alpha_m_ = 0.0;
    if (stepperPL->isParameter("Alpha_f"))
      *out << "Warning - Newmark Beta Average Acceleration ignoring Alpha_f! "
           << "Set to zero." << std::endl;
    alpha_f_ = 0.0;
    if (stepperPL->isParameter("Beta"))
      *out << "Warning - Newmark Beta Average Acceleration ignoring Beta! "
           << "Set to 0.25" << std::endl;
    beta_    = 0.25;
    if (stepperPL->isParameter("Gamma"))
      *out << "Warning - Newmark Beta Average Acceleration ignoring Gamma! "
           << "Set to 0.5" << std::endl;
    gamma_   = 0.5;

  } else if (method == "Newmark Beta Linear Acceleration") {
    // Default values
    if (stepperPL->isParameter("Rho_Inf"))
      *out << "Warning - Newmark Beta Linear Acceleration ignoring Rho_inf! "
           << "Set to zero." << std::endl;
    rho_inf_ = 0.0;
    if (stepperPL->isParameter("Alpha_m"))
      *out << "Warning - Newmark Beta Linear Acceleration ignoring Alpha_m! "
           << "Set to zero." << std::endl;
    alpha_m_ = 0.0;
    if (stepperPL->isParameter("Alpha_f"))
      *out << "Warning - Newmark Beta Linear Acceleration ignoring Alpha_f! "
           << "Set to zero." << std::endl;
    alpha_f_ = 0.0;
    if (stepperPL->isParameter("Beta"))
      *out << "Warning - Newmark Beta Linear Acceleration ignoring Beta! "
           << "Set to 0.25" << std::endl;
    beta_    = 0.25;
    if (stepperPL->isParameter("Gamma"))
      *out << "Warning - Newmark Beta Linear Acceleration ignoring Gamma! "
           << "Set to 1/6" << std::endl;
    gamma_   = 1.0/6.0;

  } else if (method == "Newmark Beta Central Difference") {
    // Default values
    if (stepperPL->isParameter("Rho_Inf"))
      *out << "Warning - Newmark Beta Central Difference ignoring Rho_inf! "
           << "Set to zero." << std::endl;
    rho_inf_ = 0.0;
    if (stepperPL->isParameter("Alpha_m"))
      *out << "Warning - Newmark Beta Central Difference ignoring Alpha_m! "
           << "Set to zero." << std::endl;
    alpha_m_ = 0.0;
    if (stepperPL->isParameter("Alpha_f"))
      *out << "Warning - Newmark Beta Central Difference ignoring Alpha_f! "
           << "Set to zero." << std::endl;
    alpha_f_ = 0.0;
    if (stepperPL->isParameter("Beta"))
      *out << "Warning - Newmark Beta Central Difference ignoring Beta! "
           << "Set to 0.0" << std::endl;
    beta_    = 0.0;
    if (stepperPL->isParameter("Gamma"))
      *out << "Warning - Newmark Beta Central Difference ignoring Gamma! "
           << "Set to 0.5" << std::endl;
    gamma_   = 0.5;

  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
       std::logic_error,
       "\nError in Tempus::StepperGeneralizedAlpha!  Invalid Scheme Name = "
       << method <<".  \n"
       <<"Valid Scheme Names are: 'Newmark Beta', 'Newmark Beta "
       <<"Average Acceleration', \n"
       <<"'Newmark Beta Linear Acceleration', and 'Newmark Beta "
       <<"Central Difference'.\n");
  }

  stepperPL->set<double>("Rho_Inf", rho_inf_);
  stepperPL->set<double>("Alpha_m", alpha_m_);
  stepperPL->set<double>("Alpha_f", alpha_f_);
  stepperPL->set<double>("Beta",    beta_   );
  stepperPL->set<double>("Gamma",   gamma_  );

  TEUCHOS_TEST_FOR_EXCEPTION( (beta_ > 1.0) || (beta_ < 0.0),
    std::logic_error,
       "\nError in 'Generalized-Alpha' stepper: invalid value of Beta = "
       << beta_ << ".  Please select Beta >= 0 and <= 1. \n");
  TEUCHOS_TEST_FOR_EXCEPTION( (gamma_ > 1.0) || (gamma_ < 0.0),
    std::logic_error,
       "\nError in 'Generalized-Alpha' stepper: invalid value of Gamma = "
       <<gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");

  if (gamma_ != 0.5 - alpha_m_ + alpha_f_) {
    *out << "Warning - may be first order\n"
         << "  gamma_ != 0.5 - alpha_m_ + alpha_f_\n"
         << "    alpha_m_ = " << alpha_m_ << "\n"
         << "    alpha_f_ = " << alpha_f_ << "\n"
         << "    beta_    = " << beta_    << "\n"
         << "    gamma_   = " << gamma_   << std::endl;
  }
  if (alpha_m_ > alpha_f_ || alpha_f_ > 0.5 ||
      beta_ < 0.25 + 0.5*(alpha_f_ - alpha_m_) ) {
    *out << "Warning - may be unstable\n"
         << "  alpha_m_ > alpha_f_ || alpha_f_ > 0.5 ||\n"
         << "  beta_ < 0.25 + 0.5*(alpha_f_ - alpha_m_\n"
         << "    alpha_m_ = " << alpha_m_ << "\n"
         << "    alpha_f_ = " << alpha_f_ << "\n"
         << "    beta_    = " << beta_    << "\n"
         << "    gamma_   = " << gamma_   << std::endl;
  }

  if (beta_ == 0.0) {
    *out << "\n \nRunning  HHT-Alpha Stepper with Beta = 0.0, which \n"
         << "specifies an explicit scheme.  WARNING: code has not been "
         << "optimized \nyet for this case (no mass lumping)\n";
  }
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::getValidParameters() const
{
  return getDefaultParameters();
}
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Solver Name", "Default Solver");

  pl->set<std::string>("Method","Generalized-Alpha");
  pl->set<double>     ("Rho_Inf", double(1.0/3.0));
  pl->set<double>     ("Alpha_m", double(-0.25));
  pl->set<double>     ("Alpha_f", double( 0.25));
  pl->set<double>     ("Beta",    double( 0.25));
  pl->set<double>     ("Gamma",   double( 0.5 ));

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperGeneralizedAlpha<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperGeneralizedAlpha_impl_hpp
