// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKBT_hpp
#define Tempus_StepperRKBT_hpp

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include "Tempus_StepperERK.hpp"
#include "Tempus_RKButcherTableau.hpp"


namespace Tempus {


// ----------------------------------------------------------------------------
/** \brief Forward Euler Runge-Kutta Butcher Tableau
 *
 *  The tableau for Forward Euler (order=1) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} 0 & 0 \\ \hline
 *                       & 1 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERKNew_ForwardEuler :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_ForwardEuler()
  {
    this->setStepperType("RKNew Forward Euler");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_ForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Forward Euler");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c = [ 0 ]'\n"
                << "A = [ 0 ]\n"
                << "b = [ 1 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    A(0,0) = ST::zero();
    b(0) = ST::one();
    c(0) = ST::zero();
    int order = 1;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Runge-Kutta 4th order Butcher Tableau
 *
 *  The tableau for RK4 (order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                        1/2 &  0  & 1/2 &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 1/3 & 1/3 & 1/6 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERKNew_4Stage4thOrder :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_4Stage4thOrder()
  {
    this->setStepperType("RKNew Explicit 4 Stage");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_4Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 4 Stage");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "\"The\" Runge-Kutta Method (explicit):\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.2, pg 138\n"
                << "c = [  0  1/2 1/2  1  ]'\n"
                << "A = [  0              ] \n"
                << "    [ 1/2  0          ]\n"
                << "    [  0  1/2  0      ]\n"
                << "    [  0   0   1   0  ]\n"
                << "b = [ 1/6 1/3 1/3 1/6 ]'";
    return Description.str();
  }

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar onethird = one/(3*one);

    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =    zero; A(0,1) =    zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) = onehalf; A(1,1) =    zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) =    zero; A(2,1) = onehalf; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =    zero; A(3,1) =    zero; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = onethird; b(2) = onethird; b(3) = onesixth;

    // fill c:
    c(0) = zero; c(1) = onehalf; c(2) = onehalf; c(3) = one;

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RKNew Bogacki-Shampine Butcher Tableau
 *
 *  The tableau (order=3(2)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A   \\ \hline
 *      & b^T \\
 *      & b^{*T}
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  & 0    &     &     & \\
 *                        1/2 & 1/2  & 0   &     & \\
 *                        3/4 & 0    & 3/4 & 0   & \\
 *                         1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                            & 2/9  & 1/3 & 4/9 & 0 \\
 *                            & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *  \f]
 *  Reference:  P. Bogacki and L.F. Shampine.
 *              A 3(2) pair of Runge–Kutta formulas.
 *              Applied Mathematics Letters, 2(4):321 – 325, 1989.
 */
template<class Scalar>
class StepperERKNew_BogackiShampine32 :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_BogackiShampine32()
  {
    this->setStepperType("Bogacki-Shampine 3(2) Pair");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_BogackiShampine32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("Bogacki-Shampine 3(2) Pair");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "P. Bogacki and L.F. Shampine.\n"
                << "A 3(2) pair of Runge–Kutta formulas.\n"
                << "Applied Mathematics Letters, 2(4):321 – 325, 1989.\n"
                << "c =     [ 0     1/2  3/4   1  ]'\n"
                << "A =     [ 0                   ]\n"
                << "        [ 1/2    0            ]\n"
                << "        [  0    3/4   0       ]\n"
                << "        [ 2/9   1/3  4/9   0  ]\n"
                << "b     = [ 2/9   1/3  4/9   0  ]'\n"
                << "bstar = [ 7/24  1/4  1/3  1/8 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onethird = one/(3*one);
    const Scalar threefourths = (3*one)/(4*one);
    const Scalar twoninths = (2*one)/(9*one);
    const Scalar fourninths = (4*one)/(9*one);

    // Fill A:
    A(0,0) =     zero; A(0,1) =        zero; A(0,2) =      zero; A(0,3) = zero;
    A(1,0) =  onehalf; A(1,1) =        zero; A(1,2) =      zero; A(1,3) = zero;
    A(2,0) =     zero; A(2,1) =threefourths; A(2,2) =      zero; A(2,3) = zero;
    A(3,0) =twoninths; A(3,1) =    onethird; A(3,2) =fourninths; A(3,3) = zero;

    // Fill b:
    b(0) = A(3,0); b(1) = A(3,1); b(2) = A(3,2); b(3) = A(3,3);

    // Fill c:
    c(0) = zero; c(1) = onehalf; c(2) = threefourths; c(3) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(7*one/(24*one));
    bstar(1) = as<Scalar>(1*one/(4*one));
    bstar(2) = as<Scalar>(1*one/(3*one));
    bstar(3) = as<Scalar>(1*one/(8*one));
    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RKNew Merson Butcher Tableau
 *
 *  The tableau (order=4(5)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A   \\ \hline
 *      & b^T \\
 *      & b^{*T}
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}  0 & 0    &     &      &     & \\
 *                        1/3 & 1/3  & 0   &      &     & \\
 *                        1/3 & 1/6  & 1/6 & 0    &     & \\
 *                        1/2 & 1/8  & 0   & 3/8  &     & \\
 *                         1  & 1/2  & 0   & -3/2 & 2   & \\ \hline
 *                            & 1/6  & 0   & 0    & 2/3 & 1/6 \\
 *                            & 1/10 & 0   & 3/10 & 2/5 & 1/5 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 4.1, pg 167.
 *
 */
template<class Scalar>
class StepperERKNew_Merson45 :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_Merson45()
  {
    this->setStepperType("Merson 4(5) Pair");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_Merson45(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("Merson 4(5) Pair");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 4.1, pg 167\n"
                << "c =     [  0    1/3  1/3  1/2   1  ]'\n"
                << "A =     [  0                       ]\n"
                << "        [ 1/3    0                 ]\n"
                << "        [ 1/6   1/6   0            ]\n"
                << "        [ 1/8    0   3/8   0       ]\n"
                << "        [ 1/2    0  -3/2   2    0  ]\n"
                << "b     = [ 1/6    0    0   2/3  1/6 ]'\n"
                << "bstar = [ 1/10   0  3/10  2/5  1/5 ]'";
    return Description.str();
  }


protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages, true);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(1,0) = as<Scalar>(one/(3*one));;

    A(2,0) = as<Scalar>(one/(6*one));;
    A(2,1) = as<Scalar>(one/(6*one));;

    A(3,0) = as<Scalar>(one/(8*one));;
    A(3,2) = as<Scalar>(3*one/(8*one));;

    A(4,0) = as<Scalar>(one/(2*one));;
    A(4,2) = as<Scalar>(-3*one/(2*one));;
    A(4,3) = 2*one;

    // Fill b:
    b(0) = as<Scalar>(one/(6*one));
    b(3) = as<Scalar>(2*one/(3*one));
    b(4) = as<Scalar>(one/(6*one));

    // Fill c:
    c(0) = zero;
    c(1) = as<Scalar>(1*one/(3*one));
    c(2) = as<Scalar>(1*one/(3*one));
    c(3) = as<Scalar>(1*one/(2*one));
    c(4) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(1*one/(10*one));
    bstar(2) = as<Scalar>(3*one/(10*one));
    bstar(3) = as<Scalar>(2*one/(5*one));
    bstar(4) = as<Scalar>(1*one/(5*one));
    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RKNew 3/8th Rule Butcher Tableau
 *
 *  The tableau (order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/3 & 1/3 &  0  &     &    \\
 *                        2/3 &-1/3 &  1  &  0  &    \\
 *                         1  &  1  & -1  &  1  &  0 \\ \hline
 *                            & 1/8 & 3/8 & 3/8 & 1/8 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.2, pg 138.
 */
template<class Scalar>
class StepperERKNew_3_8Rule :
  virtual public StepperERK<Scalar>
{
public:

  StepperERKNew_3_8Rule()
  {
    this->setStepperType("RKNew Explicit 3/8 Rule");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_3_8Rule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 3/8 Rule");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.2, pg 138\n"
                << "c = [  0  1/3 2/3  1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/3  0          ]\n"
                << "    [-1/3  1   0      ]\n"
                << "    [  1  -1   1   0  ]\n"
                << "b = [ 1/8 3/8 3/8 1/8 ]'";
    return Description.str();
  }


protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onethird     = as<Scalar>(one/(3*one));
    const Scalar twothirds    = as<Scalar>(2*one/(3*one));
    const Scalar oneeighth    = as<Scalar>(one/(8*one));
    const Scalar threeeighths = as<Scalar>(3*one/(8*one));

    // Fill A:
    A(0,0) =      zero; A(0,1) = zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) =  onethird; A(1,1) = zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) = -onethird; A(2,1) =  one; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =       one; A(3,1) = -one; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) =oneeighth; b(1) =threeeighths; b(2) =threeeighths; b(3) =oneeighth;

    // Fill c:
    c(0) = zero; c(1) = onethird; c(2) = twothirds; c(3) = one;

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit 4 Stage 3rd order by Runge
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                         1  &  0  &  1  &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 2/3 &  0  & 1/6 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class StepperERKNew_4Stage3rdOrderRunge :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_4Stage3rdOrderRunge()
  {
    this->setStepperType("RKNew Explicit 4 Stage 3rd order by Runge");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_4Stage3rdOrderRunge(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 4 Stage 3rd order by Runge");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2  1   1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/2  0          ]\n"
                << "    [  0   1   0      ]\n"
                << "    [  0   0   1   0  ]\n"
                << "b = [ 1/6 2/3  0  1/6 ]'";
    return Description.str();
  }
protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) =    zero; A(0,1) = zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) = onehalf; A(1,1) = zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) =    zero; A(2,1) =  one; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =    zero; A(3,1) = zero; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = twothirds; b(2) = zero; b(3) = onesixth;

    // Fill c:
    c(0) = zero; c(1) = onehalf; c(2) = one; c(3) = one;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit 5 Stage 3rd order by Kinnmark and Gray
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}  0  &  0  &     &     &     &    \\
 *                         1/5 & 1/5 &  0  &     &     &    \\
 *                         1/5 &  0  & 1/5 &  0  &     &    \\
 *                         1/3 &  0  &  0  & 1/3 &  0  &    \\
 *                         2/3 &  0  &  0  &  0  & 2/3 &  0 \\ \hline
 *                             & 1/4 &  0  &  0  &  0  & 3/4 \end{array}
 *  \f]
 *  Reference:  Modified by P. Ullrich.  From the prim_advance_mod.F90
 *              routine in the HOMME atmosphere model code.
 */
template<class Scalar>
class StepperERKNew_5Stage3rdOrderKandG :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_5Stage3rdOrderKandG()
  {
    this->setStepperType("RKNew Explicit 5 Stage 3rd order by Kinnmark and Gray");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_5Stage3rdOrderKandG(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 5 Stage 3rd order by Kinnmark and Gray");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Kinnmark & Gray 5 stage, 3rd order scheme \n"
                << "Modified by P. Ullrich.  From the prim_advance_mod.F90 \n"
                << "routine in the HOMME atmosphere model code.\n"
                << "c = [  0  1/5  1/5  1/3  2/3  ]'\n"
                << "A = [  0                      ]\n"
                << "    [ 1/5  0                  ]\n"
                << "    [  0  1/5   0             ]\n"
                << "    [  0   0   1/3   0        ]\n"
                << "    [  0   0    0   2/3   0   ]\n"
                << "b = [ 1/4  0    0    0   3/4  ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar onefifth = one/(5*one);
    const Scalar onefourth = one/(4*one);
    const Scalar onethird = one/(3*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar threefourths = 3*one/(4*one);
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) =     zero; A(0,1) =     zero; A(0,2) =     zero; A(0,3) =      zero; A(0,4) = zero;
    A(1,0) = onefifth; A(1,1) =     zero; A(1,2) =     zero; A(1,3) =      zero; A(1,4) = zero;
    A(2,0) =     zero; A(2,1) = onefifth; A(2,2) =     zero; A(2,3) =      zero; A(2,4) = zero;
    A(3,0) =     zero; A(3,1) =     zero; A(3,2) = onethird; A(3,3) =      zero; A(3,4) = zero;
    A(4,0) =     zero; A(4,1) =     zero; A(4,2) =     zero; A(4,3) = twothirds; A(4,4) = zero;

    // Fill b:
    b(0) =onefourth; b(1) =zero; b(2) =zero; b(3) =zero; b(4) =threefourths;

    // Fill c:
    c(0) =zero; c(1) =onefifth; c(2) =onefifth; c(3) =onethird; c(4) =twothirds;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit 3 Stage 3rd order
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                       1/2 & 1/2 &  0  &     \\
 *                        1  & -1  &  2  &  0  \\ \hline
 *                           & 1/6 & 4/6 & 1/6  \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERKNew_3Stage3rdOrder :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_3Stage3rdOrder()
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_3Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c = [  0  1/2  1  ]'\n"
                << "A = [  0          ]\n"
                << "    [ 1/2  0      ]\n"
                << "    [ -1   2   0  ]\n"
                << "b = [ 1/6 4/6 1/6 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar two = Teuchos::as<Scalar>(2*one);
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar foursixth = 4*one/(6*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =    zero; A(0,1) = zero; A(0,2) = zero;
    A(1,0) = onehalf; A(1,1) = zero; A(1,2) = zero;
    A(2,0) =    -one; A(2,1) =  two; A(2,2) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = foursixth; b(2) = onesixth;

    // fill c:
    c(0) = zero; c(1) = onehalf; c(2) = one;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit 3 Stage 3rd order TVD
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                        1  &  1  &  0  &     \\
 *                       1/2 & 1/4 & 1/4 &  0  \\ \hline
 *                           & 1/6 & 1/6 & 4/6  \end{array}
 *  \f]
 *  Reference: Sigal Gottlieb and Chi-Wang Shu,
 *             'Total Variation Diminishing Runge-Kutta Schemes',
 *             Mathematics of Computation,
 *             Volume 67, Number 221, January 1998, pp. 73-85.
 *
 *  This is also written in the following set of updates.
    \verbatim
      u1 = u^n + dt L(u^n)
      u2 = 3 u^n/4 + u1/4 + dt L(u1)/4
      u^(n+1) = u^n/3 + 2 u2/2 + 2 dt L(u2)/3
    \endverbatim
 */
template<class Scalar>
class StepperERKNew_3Stage3rdOrderTVD :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_3Stage3rdOrderTVD()
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order TVD");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_3Stage3rdOrderTVD(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order TVD");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                  << "This Stepper is known as 'RKNew Explicit 3 Stage 3rd order TVD' or 'SSPERK33'.\n"
                  << "Sigal Gottlieb and Chi-Wang Shu\n"
                  << "`Total Variation Diminishing Runge-Kutta Schemes'\n"
                  << "Mathematics of Computation\n"
                  << "Volume 67, Number 221, January 1998, pp. 73-85\n"
                  << "c = [  0   1  1/2 ]'\n"
                  << "A = [  0          ]\n"
                  << "    [  1   0      ]\n"
                  << "    [ 1/4 1/4  0  ]\n"
                  << "b = [ 1/6 1/6 4/6 ]'\n"
                  << "This is also written in the following set of updates.\n"
                  << "u1 = u^n + dt L(u^n)\n"
                  << "u2 = 3 u^n/4 + u1/4 + dt L(u1)/4\n"
                  << "u^(n+1) = u^n/3 + 2 u2/2 + 2 dt L(u2)/3";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onefourth = one/(4*one);
    const Scalar onesixth = one/(6*one);
    const Scalar foursixth = 4*one/(6*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =      zero; A(0,1) =      zero; A(0,2) = zero;
    A(1,0) =       one; A(1,1) =      zero; A(1,2) = zero;
    A(2,0) = onefourth; A(2,1) = onefourth; A(2,2) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = onesixth; b(2) = foursixth;

    // fill c:
    c(0) = zero; c(1) = one; c(2) = onehalf;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit 3 Stage 3rd order by Heun
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                       1/3 & 1/3 &  0  &     \\
 *                       2/3 &  0  & 2/3 &  0  \\ \hline
 *                           & 1/4 &  0  & 3/4  \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class StepperERKNew_3Stage3rdOrderHeun :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_3Stage3rdOrderHeun()
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order by Heun");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_3Stage3rdOrderHeun(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit 3 Stage 3rd order by Heun");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/3 2/3 ]'\n"
                << "A = [  0          ] \n"
                << "    [ 1/3  0      ]\n"
                << "    [  0  2/3  0  ]\n"
                << "b = [ 1/4  0  3/4 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onethird = one/(3*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar onefourth = one/(4*one);
    const Scalar threefourths = 3*one/(4*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =     zero; A(0,1) =      zero; A(0,2) = zero;
    A(1,0) = onethird; A(1,1) =      zero; A(1,2) = zero;
    A(2,0) =     zero; A(2,1) = twothirds; A(2,2) = zero;

    // Fill b:
    b(0) = onefourth; b(1) = zero; b(2) = threefourths;

    // fill c:
    c(0) = zero; c(1) = onethird; c(2) = twothirds;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit Midpoint
 *
 *  The tableau (order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                      1/2 & 1/2 &  0  \\ \hline
 *                          &  0  &  1   \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class StepperERKNew_Midpoint :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_Midpoint()
  {
    this->setStepperType("RKNew Explicit Midpoint");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_Midpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit Midpoint");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2 ]'\n"
                << "A = [  0      ]\n"
                << "    [ 1/2  0  ]\n"
                << "b = [  0   1  ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =    zero; A(0,1) = zero;
    A(1,0) = onehalf; A(1,1) = zero;

    // Fill b:
    b(0) = zero; b(1) = one;

    // fill c:
    c(0) = zero; c(1) = onehalf;

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RKNew Explicit Trapezoidal
 *
 *  The tableau (order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                       1  &  1  &  0  \\ \hline
 *                          & 1/2 & 1/2  \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERKNew_Trapezoidal :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_Trapezoidal()
  {
    this->setStepperType("RKNew Explicit Trapezoidal");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_Trapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RKNew Explicit Trapezoidal");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "This Stepper is known as 'RKNew Explicit Trapezoidal' or 'Heuns Method' or 'SSPERK22'.\n"
                << "c = [  0   1  ]'\n"
                << "A = [  0      ]\n"
                << "    [  1   0  ]\n"
                << "b = [ 1/2 1/2 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
   typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) = zero; A(0,1) = zero;
    A(1,0) =  one; A(1,1) = zero;

    // Fill b:
    b(0) = onehalf; b(1) = onehalf;

    // fill c:
    c(0) = zero; c(1) = one;

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Explicit RKNew Butcher Tableau
 *
 *  The tableau (stage=5, order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\ \hline
 *      & \hat{b}^T
 *  \end{array}
 *
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 *
 *
 */
template<class Scalar>
class StepperERKNew_SSPERK54 :
  virtual public StepperERK<Scalar>
{
  public:
  StepperERKNew_SSPERK54()
  {
    this->setStepperType("SSPERK54");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_SSPERK54(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("SSPERK54");
    this->setupTableau();
    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Strong Stability Preserving Explicit RKNew (stage=5, order=4)\n"
                << std::endl;
    return Description.str();
  }

protected:

  void setupTableau()
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 5;
    const int order     = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) = A(0,1) =  A(0,2) = A(0,3) = A(0,4) = zero;

    A(1,0) = as<Scalar>(0.391752226571889);
    A(1,1) = A(1,2) = A(1,3) = A(0,4) = zero;

    A(2,0) = as<Scalar>(0.217669096261169);
    A(2,1) = as<Scalar>(0.368410593050372);
    A(2,2) = A(2,3) = A(2,4) = zero;

    A(3,0) = as<Scalar>(0.082692086657811);
    A(3,1) = as<Scalar>(0.139958502191896);
    A(3,2) = as<Scalar>(0.251891774271693);
    A(3,3) = A(2,4) = zero;

    A(4,0) = as<Scalar>(0.067966283637115);
    A(4,1) = as<Scalar>(0.115034698504632);
    A(4,2) = as<Scalar>(0.207034898597385);
    A(4,3) = as<Scalar>(0.544974750228520);
    A(4,4) = zero;

    // Fill b:
    b(0) = as<Scalar>(0.146811876084786);
    b(1) = as<Scalar>(0.248482909444976);
    b(2) = as<Scalar>(0.104258830331980);
    b(3) = as<Scalar>(0.274438900901350);
    b(4) = as<Scalar>(0.226007483236908);

    // fill c:
    c(0) = zero;
    c(1) = A(1,0);
    c(2) = A(2,0) + A(2,1);
    c(3) = A(3,0) + A(3,1) + A(3,1);
    c(4) = A(4,0) + A(4,1) + A(4,2) + A(4,3);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief General Explicit Runge-Kutta Butcher Tableau
 *
 *  The format of the Butcher Tableau parameter list is
    \verbatim
      <Parameter name="A" type="string" value="# # # ;
                                               # # # ;
                                               # # #">
      <Parameter name="b" type="string" value="# # #">
      <Parameter name="c" type="string" value="# # #">
    \endverbatim
 *  Note the number of stages is implicit in the number of entries.
 *  The number of stages must be consistent.
 *
 *  Default tableau is RK4 (order=4):
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                        1/2 &  0  & 1/2 &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 1/3 & 1/3 & 1/6 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERKNew_General :
  virtual public StepperERK<Scalar>
{
public:
  StepperERKNew_General()
  {
    this->setStepperType("General ERK");
    this->setupTableau();
    this->setupDefault();
  }

  StepperERKNew_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    const int order,
    const int orderMin,
    const int orderMax,
    const Teuchos::SerialDenseVector<int,Scalar>& bstar)
  {
    this->setStepperType("General ERK");
    this->setTableau(A,b,c,order,orderMin,orderMax,bstar);

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->tableau_->isImplicit() == true, std::logic_error,
      "Error - General ERKNew received an implicit Butcher Tableau!\n");

    this->setup(appModel, stepperRKAppAction, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }

  virtual std::string getDescription() const
  {
    std::stringstream Description;
    Description << this->getStepperType() << "\n"
      << "The format of the Butcher Tableau parameter list is\n"
      << "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
      << "                                           # # # ;\n"
      << "                                           # # #\"/>\n"
      << "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
      << "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
      << "Note the number of stages is implicit in the number of entries.\n"
      << "The number of stages must be consistent.\n"
      << "\n"
      << "Default tableau is RK4 (order=4):\n"
      << "c = [  0  1/2 1/2  1  ]'\n"
      << "A = [  0              ]\n"
      << "    [ 1/2  0          ]\n"
      << "    [  0  1/2  0      ]\n"
      << "    [  0   0   1   0  ]\n"
      << "b = [ 1/6 1/3 1/3 1/6 ]'";
    return Description.str();
  }

  void setupTableau()
  {
    if (this->tableau_ == Teuchos::null) {
      // Set tableau to the default if null, otherwise keep current tableau.
      auto stepper = Teuchos::rcp(new StepperERKNew_4Stage4thOrder<Scalar>());
      auto t = stepper->getTableau();
      this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
                                 this->getStepperType(),
                                 t->A(),t->b(),t->c(),
                                 t->order(),t->orderMin(),t->orderMax(),
                                 t->bstar()));
    }
  }

  void setTableau(const Teuchos::SerialDenseMatrix<int,Scalar>& A,
                  const Teuchos::SerialDenseVector<int,Scalar>& b,
                  const Teuchos::SerialDenseVector<int,Scalar>& c,
                  const int order,
                  const int orderMin,
                  const int orderMax,
                  const Teuchos::SerialDenseVector<int,Scalar>&
                    bstar = Teuchos::SerialDenseVector<int,Scalar>())
  {
    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,orderMin,orderMax,bstar));
  }

  virtual std::string getDefaultICConsistency() const { return "Consistent"; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicERK(pl);
    pl->set<std::string>("Initial Condition Consistency",
                         this->getDefaultICConsistency());

    // Tableau ParameterList
    Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
    tableauPL->set<std::string>("A",
     "0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.0; 0.0 0.5 0.0 0.0; 0.0 0.0 1.0 0.0");
    tableauPL->set<std::string>("b",
     "0.166666666666667 0.333333333333333 0.333333333333333 0.166666666666667");
    tableauPL->set<std::string>("c", "0.0 0.5 0.5 1.0");
    tableauPL->set<int>("order", 4);
    tableauPL->set<std::string>("bstar", "");
    pl->set("Tableau", *tableauPL);

    return pl;
  }
};


} // namespace Tempus


#endif // Tempus_StepperRKButcherTableau_hpp
