// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_06.cpp
    \brief Test line search.
*/

#define USE_HESSVEC 0

#include "ROL_TestObjectives.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Test body.

  try {

    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    // Define Status Test
    Teuchos::RCP<ROL::StatusTest<RealT> > status = Teuchos::rcp(new ROL::StatusTest<RealT>(*parlist));

    *outStream << "\n\n" << ROL::ETestObjectivesToString(ROL::TESTOBJECTIVES_ROSENBROCK) << "\n\n";

    // Initial Guess Vector 
    Teuchos::RCP<std::vector<RealT> > x0_rcp = Teuchos::rcp( new std::vector<RealT> );
    ROL::StdVector<RealT> x0(x0_rcp);

    // Exact Solution Vector
    Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp( new std::vector<RealT> );
    ROL::StdVector<RealT> z(z_rcp);

    // Get Objective Function
    Teuchos::RCP<ROL::Objective<RealT> > obj = Teuchos::null;
    ROL::getTestObjectives<RealT>(obj,x0,z,ROL::TESTOBJECTIVES_ROSENBROCK);

    // Get Dimension of Problem
    int dim = Teuchos::rcp_const_cast<std::vector<RealT> >(
                (Teuchos::dyn_cast<ROL::StdVector<RealT> >(x0)).getVector())->size();
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

    // Iteration Vector
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    x.set(x0);

    // Error Vector
    Teuchos::RCP<std::vector<RealT> > e_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> e(e_rcp);
    e.zero();

    for ( ROL::EDescent desc = ROL::DESCENT_STEEPEST; desc < ROL::DESCENT_LAST; desc++ ) {
      parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", ROL::EDescentToString(desc));
      *outStream << "\n\n" << ROL::EDescentToString(desc) << "\n\n";
      for (ROL::ELineSearch ls = ROL::LINESEARCH_BACKTRACKING; ls < ROL::LINESEARCH_USERDEFINED; ls++) {
        // Define Step
        parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type",ROL::ELineSearchToString(ls));
        Teuchos::RCP<ROL::LineSearchStep<RealT> > step = Teuchos::rcp(new ROL::LineSearchStep<RealT>(*parlist));
      
        // Define Algorithm
        ROL::Algorithm<RealT> algo(step,status,false);

        // Run Algorithm
        x.set(x0);
        algo.run(x, *obj, true, *outStream);

        // Compute Error
        e.set(x);
        e.axpy(-1.0,z);
        *outStream << "\nNorm of Error: " << e.norm() << "\n";
        //errorFlag += (int)(e.norm() < std::sqrt(ROL::ROL_EPSILON)); 
      }
      
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

