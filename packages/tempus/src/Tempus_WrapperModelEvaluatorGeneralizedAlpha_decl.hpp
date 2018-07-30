// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapperModelEvaluatorGeneralizedAlpha_decl_hpp
#define Tempus_WrapperModelEvaluatorGeneralizedAlpha_decl_hpp

#include <functional>
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application application ModelEvaluator, \f$ f(\dot{x},x,t) \f$.
 *
 *  This class breaks the primary design principle for ModelEvaluators;
 *  it is not stateless!
 */
template <typename Scalar>
class WrapperModelEvaluatorGeneralizedAlpha
 : public Tempus::WrapperModelEvaluator<Scalar>
{
public:
  typedef Thyra::VectorBase<Scalar>  Vector;

  /// Constructor
  WrapperModelEvaluatorGeneralizedAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
    : appModel_(appModel)
  {}

  /// Set the underlying application ModelEvaluator
  void setAppModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me)
  { appModel_ = me; }

  /// Get the underlying application model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAppModel() const
  { return appModel_; }

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
    { return appModel_->create_W_op(); }

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    get_W_factory() const
    { return appModel_->get_W_factory(); }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const
    { return appModel_->get_f_space(); }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const
    { return appModel_->get_p_space(p); };

    Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const
    { return appModel_->get_p_names(p); }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const
    { return appModel_->get_x_space(); }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int i) const
    { return appModel_->get_g_space(i); }

    Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const
    {
      Thyra::ModelEvaluatorBase::InArgs<Scalar>
        nominalValues = appModel_->getNominalValues();
      return nominalValues;
    }

    /// Set InArgs the wrapper ModelEvalutor.
    virtual void setInArgs(Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs)
    { wrapperInArgs_.setArgs(inArgs); }

    /// Get InArgs the wrapper ModelEvalutor.
    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getInArgs()
    { return wrapperInArgs_; }

    /// Set OutArgs the wrapper ModelEvalutor.
    virtual void setOutArgs(Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs)
    { wrapperOutArgs_.setArgs(outArgs); }

    /// Get OutArgs the wrapper ModelEvalutor.
    virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar> getOutArgs()
    { return wrapperOutArgs_; }

    /// Set parameters for application implicit ModelEvaluator solve.
    virtual void setForSolve(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
      Thyra::ModelEvaluatorBase::InArgs<Scalar>   inArgs,
      Thyra::ModelEvaluatorBase::OutArgs<Scalar>  outArgs)
    {
      timeDer_ = timeDer;
      wrapperInArgs_.setArgs(inArgs);
      wrapperOutArgs_.setArgs(outArgs);
    }

    Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

    void evalModelImpl(
              const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;
  //@}

    enum SCHEME_TYPE {NEWMARK_IMPLICIT_AFORM, NEWMARK_IMPLICIT_DFORM};

private:

  /// Default constructor - not allowed
  WrapperModelEvaluatorGeneralizedAlpha();

private:

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;
  Teuchos::RCP<TimeDerivative<Scalar> >              timeDer_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar>          wrapperInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         wrapperOutArgs_;

};

} // namespace Tempus

#endif // Tempus_WrapperModelEvaluatorGeneralizedAlpha_hpp
