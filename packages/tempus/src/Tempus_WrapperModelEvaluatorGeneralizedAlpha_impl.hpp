// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapperModelEvaluatorGeneralizedAlpha_impl_hpp
#define Tempus_WrapperModelEvaluatorGeneralizedAlpha_impl_hpp

namespace Tempus {


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorGeneralizedAlpha<Scalar>::
createInArgs() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(appModel_->Np());
  inArgs.setSupports(MEB::IN_ARG_x);

  return inArgs;
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorGeneralizedAlpha<Scalar>::
createOutArgsImpl() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(appModel_->Np(),0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);

  return outArgs;
}


template <typename Scalar>
void
WrapperModelEvaluatorGeneralizedAlpha<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Thyra::VectorBase;

  auto x       = rcp_const_cast<VectorBase<Scalar>>(inArgs.get_x());
  auto xDot    = rcp_const_cast<VectorBase<Scalar>>(inArgs.get_x_dot());
  auto xDotDot = rcp_const_cast<VectorBase<Scalar>>(inArgs.get_x_dot_dot());

  timeDer_->compute(x, xDot, xDotDot);

  MEB::InArgs<Scalar>  appInArgs (wrapperInArgs_);
  appInArgs.set_x        (x);
  appInArgs.set_x_dot    (xDot);
  appInArgs.set_x_dot_dot(xDotDot);
  for (int i = 0; i < appModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      appInArgs.set_p(i, inArgs.get_p(i));
  }

  MEB::OutArgs<Scalar> appOutArgs(wrapperOutArgs_);
  appOutArgs.set_f(outArgs.get_f());
  appOutArgs.set_W_op(outArgs.get_W_op());
  if (outArgs.supports(MEB::OUT_ARG_W_prec)) {
    appOutArgs.set_W_prec(outArgs.get_W_prec());
  }

  // build residual and jacobian
  appModel_->evalModel(appInArgs, appOutArgs);
}


} // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorGeneralizedAlpha_impl_hpp
