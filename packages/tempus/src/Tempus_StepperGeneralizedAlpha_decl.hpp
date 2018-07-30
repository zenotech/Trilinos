// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperGeneralizedAlpha_decl_hpp
#define Tempus_StepperGeneralizedAlpha_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorSecondOrder.hpp"

namespace Tempus {


/** \brief Generalized-Alpha time stepper.
 *
 * Here, we implement the Generalized-\f$\alpha\f$ method from
 * J. Chung and G. M. Hulbert. "A time integration algorithm for
 * structural dynamics with improved numerical dissipation:
 * The generalized-alpha method. Journal of Applied Mechanics,
 * 60:371–375, June 1993.
 *
 * The Generalized-\f$\alpha\f$ method can be written
 *  \f[
 *    \mathbf{M}\,\ddot{\mathbf{x}}^{n+1-\alpha_{m}}+\mathbf{C}\,
 *    \dot{\mathbf{x}}^{n+1-\alpha_{f}}+\mathbf{K}\,
 *     \mathbf{x}^{n+1-\alpha_{f}}=\mathbf{F}\left(t^{n+1-\alpha_{f}}\right)
 *  \f]
 *  where the vectors,
 *  \f$\psi=[\ddot{\mathbf{x}},\dot{\mathbf{x}},\mathbf{x},\mathbf{F}]\f$
 *  are evaluated by
 *  \f{eqnarray*}{
 *    \psi^{n+1-\alpha_{f}} & = & \alpha_{f}\psi^{n}+(1-\alpha_{f})\psi^{n+1}\\
 *    \psi^{n+1-\alpha_{m}} & = & \alpha_{m}\psi^{n}+(1-\alpha_{m})\psi^{n+1}
 *  \f}
 *  and the time evaluations are similarly found by
 *  \f[
 *    t^{n+1-\alpha_{f}} = \alpha_{f}t^{n}+(1-\alpha_{f})t^{n+1}
 *  \f]
 *  The solution and derivative updates are
 *  \f[
 *    \mathbf{x}^{n+1}=\mathbf{x}^{n}+\Delta t^{n}\dot{\mathbf{x}}^{n}
 *    +\left(\Delta t^{n}\right)^{2}\left[\left(\frac{1}{2}
 *    -\beta\right)\ddot{\mathbf{x}}^{n}+\beta\ddot{\mathbf{x}}^{n+1}\right]
 *  \f]
 *  \f[
 *     \dot{\mathbf{x}}^{n+1}=\dot{\mathbf{x}}^{n}
 *     +\Delta t^{n}\left[(1-\gamma)\mathbf{\ddot{x}}^{n}
 *     +\gamma\mathbf{\ddot{x}}^{n+1}\right]
 *  \f]
 *  The solution steps for Generalized-\f$\alpha\f$ method are then to
 *  implicitly solve the governing equation with the above updates.
 *  The Generalized-\f$\alpha\f$ method order of accuracy can be expressed as
 *  \f[
 *    {\cal O}=\left\{
 *    \begin{array}{ll}
 *    \Delta t^{2} & \mbox{if }\gamma=\frac{1}{2}-\alpha_{m}+\alpha_{f}\\
 *    \Delta t & \mbox{otherwise}
 *    \end{array}
 *    \right.,
 *  \f]
 *  and unconditional stability will be obtained if
 *  \f[
 *    \alpha_{m}\leq\alpha_{f}\le\frac{1}{2}\,\,\,\,\mbox{ and }\,\,\,\,
 *    \beta\ge\frac{1}{4}+\frac{1}{2}\left(\alpha_{f}-\alpha_{m}\right).
 *  \f]
 *  To maximize high-frequency dissipation and minimize the low-frequency
 *  dissipation, the Generalized-\f$\alpha\f$ method requires
 *  \f[
 *    \beta=\frac{1}{4}\left(1-\alpha_{m}+\alpha_{f}\right)^{2}.
 *  \f]
 *  These properties can be maintained and the amount of dissipation can
 *  be specified through the spectral radius in the high-frequency limit,
 *  \f$\rho_{\infty}\f$, where
 *  \f[
 *    \alpha_{m}=\frac{2\rho_{\infty}-1}{\rho_{\infty}+1}
 *    \,\,\,\,\mbox{ and }\,\,\,\,
 *    \alpha_{f}=\frac{\rho_{\infty}}{\rho_{\infty}+1}
 *  \f]
 *  with
 *  \f[
 *    0\le\rho_{\infty}\le1
 *  \f]
 *  and \f$\rho_{\infty}=0\f$ is the asymptotic annihilation case
 *  (high-frequency response is annihilated after one time step),
 *  and \f$\rho_{\infty}=1\f$ is the no dissipation case.
 *  The Generalized-\f$\alpha\f$ method can be expressed as
 *  \f{eqnarray*}{
 *  \alpha_{m} & = & \frac{2\rho_{\infty}-1}{\rho_{\infty}+1}\\
 *  \alpha_{f} & = & \frac{\rho_{\infty}}{\rho_{\infty}+1}\\
 *  \beta & = & \frac{1}{4}\left(1-\alpha_{m}+\alpha_{f}\right)^{2}\\
 *  \gamma & = & \frac{1}{2}-\alpha_{m}+\alpha_{f}
 *  \f}
 *
 *  If we are solving for \f$\mathbf{x}^{n+1}\f$, the initial guess for
 *  \f$\mathbf{x}_{0}^{n+1-\alpha_{f}}\f$ is interpolated with the initial
 *  guess of \f$\mathbf{x}_{0}^{n+1}\f$,
 *  \f[
 *    \mathbf{x}_{0}^{n+1-\alpha_{f}}=
 *      \alpha_{f}\,\mathbf{x}^{n}+(1-\alpha_{f})\,\mathbf{x}_{0}^{n+1}
 *  \f]
 *  and in subsequent iterations, \f$\mathbf{x}_{(k)}^{n+1-\alpha_{f}}\f$
 *  will be determined by the solver. Updates for
 *  \f$\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}\f$
 *  and \f$\ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}}\f$ can be found by
 *  \f{eqnarray*}{
 *    \ddot{\mathbf{x}}_{(k)}^{n+1} & = &
 *       \frac{\mathbf{x}_{(k)}^{n+1-\alpha_{f}}
 *       -\mathbf{x}^{n}}{(1-\alpha_{f})\beta\left(\Delta t^{n}\right)^{2}}
 *       -\frac{\dot{\mathbf{x}}^{n}}{\beta\Delta t^{n}}
 *       +\left(1-\frac{1}{2\beta}\right)\ddot{\mathbf{x}}^{n}\\
 *    \dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}} & = &
 *       \dot{\mathbf{x}}^{n}+(1-\alpha_{f})\,
 *       \Delta t^{n}\left[(1-\gamma)\mathbf{\ddot{x}}^{n}
 *                         +\gamma\mathbf{\ddot{x}}_{(k)}^{n+1}\right]\\
 *    \ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}} & = &
 *       \alpha_{m}\ddot{\mathbf{x}}^{n}
 *       +(1-\alpha_{m})\ddot{\mathbf{x}}_{(k)}^{n+1}
 *  \f}
 *
 *  If we are solving for \f$\dot{\mathbf{x}}^{n+1}\f$, the initial guess
 *  for \f$\mathbf{\dot{\mathbf{x}}}_{0}^{n+1-\alpha_{f}}\f$ is interpolated
 *  with the initial guess of \f$\dot{\mathbf{x}}_{0}^{n+1}\f$,
 *  \f[
 *    \dot{\mathbf{x}}_{0}^{n+1-\alpha_{f}}=
 *      \alpha_{f}\,\dot{\mathbf{x}}^{n}
 *      +(1-\alpha_{f})\,\dot{\mathbf{x}}_{0}^{n+1}
 *  \f]
 *  and in subsequent iterations, \f$\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}\f$
 *  will be determined by the solver. Updates for
 *  \f$\mathbf{x}_{(k)}^{n+1-\alpha_{f}}\f$
 *  and \f$\ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}}\f$ can be found by
 *  \f{eqnarray*}{
 *    \mathbf{\ddot{x}}_{(k)}^{n+1} & = &
 *      \frac{\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}
 *      -\dot{\mathbf{x}}^{n}}{(1-\alpha_{f})\gamma\Delta t^{n}}
 *      +(1-\frac{1}{\gamma})\mathbf{\ddot{x}}^{n}\\
 *    \mathbf{x}_{(k)}^{n+1-\alpha_{f}} & = &
 *      \mathbf{x}^{n}+(1-\alpha_{f})\,
 *      \left\{ \Delta t^{n}\dot{\mathbf{x}}^{n}+\left(\Delta t^{n}\right)^{2}
 *      \left[\left(\frac{1}{2}-\beta\right)\ddot{\mathbf{x}}^{n}
 *      +\beta\ddot{\mathbf{x}}_{(k)}^{n+1}\right]\right\} \\
 *    \ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}} & = &
 *      \alpha_{m}\ddot{\mathbf{x}}^{n}
 *      +(1-\alpha_{m})\ddot{\mathbf{x}}_{(k)}^{n+1}
 *  \f}
 *
 *  If we are solving for \f$\ddot{\mathbf{x}}^{n+1}\f$, the initial guess
 *  for \f$\mathbf{\ddot{\mathbf{x}}}_{0}^{n+1-\alpha_{m}}\f$ is interpolated
 *  with the initial guess of \f$\ddot{\mathbf{x}}_{0}^{n+1}\f$,
 *  \f[
 *    \mathbf{\ddot{\mathbf{x}}}_{0}^{n+1-\alpha_{m}}=
 *    \alpha_{m}\ddot{\mathbf{x}}^{n}+(1-\alpha_{m})\ddot{\mathbf{x}}_{0}^{n+1}
 *  \f]
 *  and in subsequent iterations,
 *  \f$\mathbf{\ddot{\mathbf{x}}}_{0}^{n+1-\alpha_{m}}\f$
 *  will be determined by the solver. Updates for
 *  \f$\mathbf{x}_{(k)}^{n+1-\alpha_{f}}\f$
 *  and \f$\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}\f$ can be found by
 *  \f{eqnarray*}{
 *    \ddot{\mathbf{x}}_{(k)}^{n+1} & = &
 *      \frac{\ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}}
 *      -\alpha_{m}\,\ddot{\mathbf{x}}^{n}}{(1-\alpha_{m})}\\
 *    \mathbf{x}_{(k)}^{n+1-\alpha_{f}} & = &
 *      \mathbf{x}^{n}+(1-\alpha_{f})\,
 *      \left\{ \Delta t^{n}\dot{\mathbf{x}}^{n}
 *      +\left(\Delta t^{n}\right)^{2}\left[\left(\frac{1}{2}-\beta\right)
 *      \ddot{\mathbf{x}}^{n}+\beta\ddot{\mathbf{x}}_{(k)}^{n+1}
 *      \right]\right\} \\
 *    \dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}} & = &
 *      \dot{\mathbf{x}}^{n}+(1-\alpha_{f})\,
 *      \Delta t^{n}\left[(1-\gamma)\mathbf{\ddot{x}}^{n}
 *      +\gamma\mathbf{\ddot{x}}_{(k)}^{n+1}\right]
 *  \f}
 *  After the solver has iterated to convergence, the solution at
 *  \f$t^{n+1}\f$ can be found via
 *  \f{eqnarray*}{
 *    \mathbf{x}^{n+1} & = &
 *      \frac{\mathbf{x}_{(k)}^{n+1-\alpha_{f}}
 *      -\alpha_{f}\,\mathbf{x}^{n}}{(1-\alpha_{f})}\\
 *    \dot{\mathbf{x}}^{n+1} & = &
 *      \frac{\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}
 *      -\alpha_{f}\,\dot{\mathbf{x}}^{n}}{(1-\alpha_{f})}\\
 *    \ddot{\mathbf{x}}^{n+1} & = &
 *      \frac{\ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}}
 *      -\alpha_{m}\,\ddot{\mathbf{x}}^{n}}{(1-\alpha_{m})}
 *  \f}
 *
 *  <table>
 *  <caption id="multi_row">Methods for second-order systems which are
 *     second-order accurate and unconditionally stable.</caption>
 *  <tr><th> Integration Methods
 *      <th> \f$\alpha_{m}\f$
 *      <th> \f$\alpha_{f}\f$
 *      <th> \f$\beta\f$
 *      <th> \f$\gamma\f$
 *  <tr><td> Generalized-\f$\alpha\f$ (\f$0\leq\rho_{\infty}\leq1\f$)
 *      <td> \f$-1\le\frac{2\rho_{\infty}-1}{\rho_{\infty}+1}\le\frac{1}{2}\f$
 *      <td> \f$0\le\frac{\rho_{\infty}}{\rho_{\infty}+1}\le\frac{1}{2}\f$
 *      <td> \f$\frac{1}{4}\left(1-\alpha_{m}+\alpha_{f}\right)^{2}\f$
 *      <td> \f$\frac{1}{2}-\alpha_{m}+\alpha_{f}\f$
 *  <tr><td> WBZ-\f$\alpha\f$ (\f$0\leq\rho_{\infty}\leq1\f$)
 *      <td> \f$-1\le\frac{\rho_{\infty}-1}{\rho_{\infty}+1}\le0\f$
 *      <td> 0
 *      <td> \f$\frac{1}{4}\left(1-\alpha_{m}\right)^{2}\f$
 *      <td> \f$\frac{1}{2}-\alpha_{m}\f$
 *  <tr><td> HHT-\f$\alpha\f$ (\f$1\ge\rho_{\infty}\ge1/2\f$)
 *      <td> 0
 *      <td> \f$0\le\frac{1-\rho_{\infty}}{1+\rho_{\infty}}\le\frac{1}{3}\f$
 *      <td> \f$\frac{1}{4}\left(1+\alpha_{f}\right)^{2}\f$
 *      <td> \f$\frac{1}{2}+\alpha_{f}\f$
 *  <tr><td> Newmark-\f$\beta\f$ (Trapezoidal rule)
 *      <td> 0
 *      <td> 0
 *      <td> \f$\frac{1}{4}\f$
 *      <td> \f$\frac{1}{2}\f$
 *  </table>
 *
 *  <table>
 *  <caption id="multi_row">Accuracy and stability requirements for
 *     second-order systems.</caption>
 *  <tr><th> Method <th> Accuracy <th> Unconditional Stability <th> Comments
 *  <tr><td> Generalized-\f$\alpha\f$
 *      <td> \f[ \left\{ \begin{array}{ll}
 *           \Delta t^{2} & \mbox{if }\gamma=\frac{1}{2}-\alpha_{m}+\alpha_{f}\\
 *           \Delta t & \mbox{otherwise}
 *           \end{array} \right. \f]
 *      <td> \f[ \begin{array}{c}
 *           \alpha_{m}\leq\alpha_{f}\le\frac{1}{2}\,\,\mbox{ and}\\
 *           \beta\ge\frac{1}{4}+\frac{1}{2}\left(\alpha_{f}-\alpha_{m}\right)
 *           \end{array} \f]
 *      <td>
 *  <tr><td> WBZ-\f$\alpha\f$
 *      <td> \f[ \left\{ \begin{array}{ll}
 *           \Delta t^{2} & \mbox{if }\gamma=\frac{1}{2}-\alpha_{m}\\
 *           \Delta t & \mbox{otherwise}
 *           \end{array} \right. \f]
 *      <td> \f[ \begin{array}{c}
 *           \alpha_{m}\leq\frac{1}{2}\,\,\mbox{ and}\\
 *           \beta\ge\frac{1}{4}-\frac{\alpha_{m}}{2}
 *           \end{array} \f]
 *      <td>
 *  <tr><td> HHT-\f$\alpha\f$
 *      <td> \f[ \left\{ \begin{array}{ll}
 *           \Delta t^{2} & \mbox{if }\gamma=\frac{1}{2}+\alpha_{f}\\
 *           \Delta t & \mbox{otherwise}
 *           \end{array} \right. \f]
 *      <td> \f[ \begin{array}{c}
 *           \alpha_{m}\leq\alpha_{f}\le\frac{1}{2}\,\,\mbox{ and}\\
 *           \beta\ge\frac{1}{4}+\frac{\alpha_{f}}{2}
 *           \end{array} \f]
 *      <td>
 *  <tr><td> Newmark-\f$\beta\f$
 *      <td> \f[ \left\{ \begin{array}{ll}
 *           \Delta t^{2} & \mbox{if }\gamma=\frac{1}{2}\\
 *           \Delta t & \mbox{otherwise}
 *           \end{array} \right. \f]
 *      <td> \f[ \begin{array}{c}
 *             \beta\ge\frac{\gamma}{2}\ge\frac{1}{4}
 *           \end{array} \f]
 *      <td> \f[ \left\{ \begin{array}{ll}
 *              \gamma=\frac{1}{2},\,\beta=\frac{1}{4} &
 *                \mbox{Constant average acceleration}\\
 *              \gamma=\frac{1}{2},\,\beta=\frac{1}{6} &
 *                \mbox{Linear acceleration}
 *           \end{array} \right. \f]
 *  </table>
 *
 */
template<class Scalar>
class StepperGeneralizedAlpha : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperGeneralizedAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null){}

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {
      if (this->stepperPL_->template get<double>("Gamma"  ) ==
          0.5 - this->stepperPL_->template get<double>("Alpha_m")
          + this->stepperPL_->template get<double>("Alpha_f")) return 2.0;
      else return 1.0;
    }
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
  //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

private:

  /// Default Constructor -- not allowed
  StepperGeneralizedAlpha();

  Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDotScratch_;
};


/** \brief Time-derivative interface for Generalized-Alpha.
 *
 *  If we are solving for \f$\mathbf{x}^{n+1}\f$, the initial guess for
 *  \f$\mathbf{x}_{0}^{n+1-\alpha_{f}}\f$ is interpolated with the initial
 *  guess of \f$\mathbf{x}_{0}^{n+1}\f$,
 *  \f[
 *    \mathbf{x}_{0}^{n+1-\alpha_{f}}=
 *      \alpha_{f}\,\mathbf{x}^{n}+(1-\alpha_{f})\,\mathbf{x}_{0}^{n+1}
 *  \f]
 *  and in subsequent iterations, \f$\mathbf{x}_{(k)}^{n+1-\alpha_{f}}\f$
 *  will be determined by the solver. Updates for
 *  \f$\dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}}\f$
 *  and \f$\ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}}\f$ can be found by
 *  \f{eqnarray*}{
 *    \ddot{\mathbf{x}}_{(k)}^{n+1} & = &
 *       \frac{\mathbf{x}_{(k)}^{n+1-\alpha_{f}}
 *       -\mathbf{x}^{n}}{(1-\alpha_{f})\beta\left(\Delta t^{n}\right)^{2}}
 *       -\frac{\dot{\mathbf{x}}^{n}}{\beta\Delta t^{n}}
 *       +\left(1-\frac{1}{2\beta}\right)\ddot{\mathbf{x}}^{n}\\
 *    \dot{\mathbf{x}}_{(k)}^{n+1-\alpha_{f}} & = &
 *       \dot{\mathbf{x}}^{n}+(1-\alpha_{f})\,
 *       \Delta t^{n}\left[(1-\gamma)\mathbf{\ddot{x}}^{n}
 *                         +\gamma\mathbf{\ddot{x}}_{(k)}^{n+1}\right]\\
 *    \ddot{\mathbf{x}}_{(k)}^{n+1-\alpha_{m}} & = &
 *       \alpha_{m}\ddot{\mathbf{x}}^{n}
 *       +(1-\alpha_{m})\ddot{\mathbf{x}}_{(k)}^{n+1}
 *  \f}
 */
template <typename Scalar>
class StepperGeneralizedAlphaTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperGeneralizedAlphaTimeDerivative(
    Scalar                                         dt,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotDotOld,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDotScratch,
    Scalar rho_inf,
    Scalar alpha_m,
    Scalar alpha_f,
    Scalar beta,
    Scalar gamma)
    : dt_             (dt),
      xOld_           (xOld),
      xDotOld_        (xDotOld),
      xDotDotOld_     (xDotDotOld),
      xDotDotScratch_ (xDotDotScratch),
      rho_inf_        (rho_inf),
      alpha_m_        (alpha_m),
      alpha_f_        (alpha_f),
      beta_           (beta),
      gamma_          (gamma)
  {}

  /// Destructor
  virtual ~StepperGeneralizedAlphaTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    // If solving for x ...
    {
      // xDotDotScratch_
      Scalar c = 1.0/((1.0-alpha_f_)*beta_*dt_*dt_);
      Thyra::V_StVpStV(xDotDotScratch_.ptr(), c, *x, -c , *xOld_);
      c = -1.0/(beta_*dt_);
      Thyra::Vp_StV(xDotDotScratch_.ptr(), c, *xDotOld_);
      c = 1.0 - 1.0/(2.0*beta_);
      Thyra::Vp_StV(xDotDotScratch_.ptr(), c, *xDotDotOld_);

      // xDot
      const Scalar c1 = (1.0-alpha_f_)*dt_*(1.0-gamma_);
      const Scalar c2 = (1.0-alpha_f_)*dt_*gamma_;
      Thyra::V_StVpStV(xDot.ptr(), c1, *xDotDotOld_, c2, *xDotDotScratch_);
      Thyra::Vp_V(xDot.ptr(), *xDotOld_);

      // xDotDot
      Thyra::V_StVpStV(xDotDot.ptr(), alpha_m_, *xDotDotOld_,
                                (1.0-alpha_m_), *xDotDotScratch_);
    }
  }

private:

  Scalar                                         dt_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotOld_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotDotOld_;
  Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDotScratch_;

  Scalar rho_inf_;
  Scalar alpha_m_;
  Scalar alpha_f_;
  Scalar beta_;
  Scalar gamma_;

};

} // namespace Tempus

#endif // Tempus_StepperGeneralizedAlpha_decl_hpp
