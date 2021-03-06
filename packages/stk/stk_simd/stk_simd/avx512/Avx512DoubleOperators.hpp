// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// IWYU pragma: private, include <stk_simd/Simd.hpp>

namespace stk {
namespace simd {

namespace hidden {
static const simd::Bool TRUE_VAL(true);
}

inline simd::Double operator+ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_add_pd(a._data, b._data));
}

inline simd::Double operator- (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_sub_pd(a._data, b._data));
}

inline simd::Double operator* (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_mul_pd(a._data, b._data));
}
   
inline simd::Double operator/ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm512_div_pd(a._data, b._data));
}

inline simd::Double operator+ (const double a, const simd::Double& b) {
  return simd::Double(_mm512_add_pd(_mm512_set1_pd(a),b._data));
}

inline simd::Double operator+ (const simd::Double& a, const double b) {
  return simd::Double(_mm512_add_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator- (const double a, const simd::Double& b) {
  return simd::Double(_mm512_sub_pd(_mm512_set1_pd(a),b._data));
}
  
inline simd::Double operator- (const simd::Double& a, const double b) {
  return simd::Double(_mm512_sub_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator* (const double a, const simd::Double& b) {
  return simd::Double(_mm512_mul_pd(_mm512_set1_pd(a),b._data));
}

inline simd::Double operator* (const simd::Double& a, const double b) {
  return simd::Double(_mm512_mul_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Double operator/ (const double a, const simd::Double& b) {
  return simd::Double(_mm512_div_pd(_mm512_set1_pd(a),b._data));
}
  
inline simd::Double operator/ (const simd::Double& a, const double b) {
  return simd::Double(_mm512_div_pd(a._data, _mm512_set1_pd(b)));
}

inline simd::Bool operator< (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmplt_pd_mask(x._data, y._data));
}
    
inline simd::Bool operator<=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmple_pd_mask(x._data, y._data));
}

inline simd::Bool operator> (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmplt_pd_mask(y._data, x._data));
}

inline simd::Bool operator>=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmple_pd_mask(y._data, x._data));
}

inline simd::Bool operator==(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmpeq_pd_mask(x._data, y._data));
}
  
inline simd::Bool operator!=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm512_cmpneq_pd_mask(x._data, y._data));
}

inline simd::Bool operator&& (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool( (__mmask8)_mm512_kand(x._data, y._data));
}

inline simd::Bool operator|| (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool( (__mmask8)_mm512_kor(x._data, y._data));
}

inline simd::Bool operator! (const simd::Bool& x) {
  return simd::Bool( (__mmask8)_mm512_kxor(hidden::TRUE_VAL._data, x._data) );
}

} // namespace simd
} // namespace stk
