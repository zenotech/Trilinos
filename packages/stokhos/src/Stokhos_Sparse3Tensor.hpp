// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SPARSE3TENSOR_HPP
#define STOKHOS_SPARSE3TENSOR_HPP

#include <ostream>

#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*! 
   * \brief Data structure storing a sparse 3-tensor C(i,j,k) in a 
   * a compressed format.
   */
  template <typename ordinal_type, typename value_type>
  class Sparse3Tensor {
  public:
    
    //! Constructor
    Sparse3Tensor(ordinal_type sz);
    
    //! Destructor
    ~Sparse3Tensor();

    //! Return size
    ordinal_type size() const;

    //! Return number of n of non-zero's in C(:,:,k) for a given k
    ordinal_type num_values(ordinal_type k) const;
      
    //! Get value C(i,j,k) for 0 <= l < n for a given k
    void value(ordinal_type k, ordinal_type l, 
	       ordinal_type& i, ordinal_type& j, value_type& c) const;

    //! Add new term for given (i,j,k)
    void add_term(ordinal_type i, ordinal_type j, ordinal_type k,
		  const value_type& c);

    //! Print tensor
    void print(std::ostream& os) const;

  private:

    // Prohibit copying
    Sparse3Tensor(const Sparse3Tensor&);

    // Prohibit Assignment
    Sparse3Tensor& operator=(const Sparse3Tensor& b);

  protected:

    //! i-indices in Cijk for each k
    Teuchos::Array< Teuchos::Array<ordinal_type> > i_indices;

    //! j-indices in Cijk for each k
    Teuchos::Array< Teuchos::Array<ordinal_type> > j_indices;

    //! values in Cijk for each k
    Teuchos::Array< Teuchos::Array<value_type> > Cijk_values;

  }; // class Sparse3Tensor

  template <typename ordinal_type, typename value_type>
  std::ostream& 
  operator << (std::ostream& os, 
	       const Sparse3Tensor<ordinal_type, value_type>& Cijk) {
    Cijk.print(os);
    return os;
  }

} // namespace Stokhos

// Include template definitions
#include "Stokhos_Sparse3TensorImp.hpp"

#endif // STOKHOS_SPARSE3TENSOR_HPP
