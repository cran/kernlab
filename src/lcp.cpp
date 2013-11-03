/* ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 2.0
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is the Suffix Array based String Kernel.
 *
 * The Initial Developer of the Original Code is
 * Statistical Machine Learning Program (SML), National ICT Australia (NICTA).
 * Portions created by the Initial Developer are Copyright (C) 2006
 * the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 *
 *   Choon Hui Teo <ChoonHui.Teo@rsise.anu.edu.au>
 *   S V N Vishwanathan <SVN.Vishwanathan@nicta.com.au>
 *
 * ***** END LICENSE BLOCK ***** */


// File    : sask/Code/LCP.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           11 Oct 2006


#ifndef LCP_CPP
#define LCP_CPP

#include "lcp.h"

// Threshold for compacting LCP[]
const Real THRESHOLD = 0.3;

LCP::LCP(const UInt32 &size): _p_array(0), 
															_idx_array(0), 
															_val_array(0),
															_size(size),
                              _is_compact(false),
                              _beg(0),
                              _end(0),
                              _cache(0),
                              _dist(0),
															array(new UInt32[size]){ }

LCP::~LCP()
{
  if(array) {delete [] array; array = 0;}
  if(_p_array) {delete [] _p_array; _p_array = 0;}
  if(_idx_array) {delete [] _idx_array; _idx_array = 0;}
  if(_val_array) {delete [] _val_array; _val_array = 0;}
}

/**
 *  Compact initial/original lcp array of n elements (i.e. 4n bytes)
 *  into a n byte array with 8 bytes of secondary storage. 
 *  
 */
ErrorCode  
LCP::compact(void){
  
  // Validate pre-conditions
	//assert(!array.empty() && array.size() == _size);
  assert(array);
  
  // Already compact. Nothing to do
  if (_is_compact) return NOERROR;
  
  // Count number of lcp-values >= 255.
//   UInt32 idx_len = std::count_if(array.begin(), array.end(), 
//                                  std::bind2nd(std::greater<int>(),254));


#ifdef _RWSTD_NO_CLASS_PARTIAL_SPEC
  UInt32 idx_len = 0;
  std::count_if(array, array + _size, std::bind2nd(std::greater<int>(),254),
		idx_len);
#else
  UInt32 idx_len = std::count_if(array, array + _size, 
                                 std::bind2nd(std::greater<int>(),254));
#endif

  // Compact iff  idx_len/|array| > THRESHOLD
 
	if((Real)idx_len/_size > THRESHOLD) {
		//std::cout<< "Not compacting " << std::endl;
    return NOERROR;
  }

	// std::cout<< "Compacting with : " << idx_len << std::endl;
	// We know how much space to use
//   _p_array.resize(_size);
//   _idx_array.resize(idx_len);
//   _val_array.resize(idx_len);

  _p_array   = new Byte1[_size];
  _idx_array = new UInt32[idx_len];
  _val_array = new UInt32[idx_len];

  // Hold pointers for later. Avoids function calls
//   _beg = _idx_array.begin();
//   _end = _idx_array.end();
//   _cache = _idx_array.begin();

  _beg   = _idx_array;
  _end   = _idx_array + idx_len;
  _cache = _idx_array;
 	_dist  = 0;
		

  for(UInt32 i=0, j=0; i<_size; i++) {
    if(array[i] < 255){
      _p_array[i] = array[i];
    }else {
      _p_array[i] = 255;
      _idx_array[j] = i;
      _val_array[j] = array[i];
      j++;
    }
  }
  //array.resize(0);
//   array.clear();
  delete [] array;
  array = 0;
	
	_is_compact = true;

  return NOERROR;
}

/**
 *  Retrieve lcp array values.
 *
 *  \param idx - (IN) Index of lcp array
 */
UInt32
LCP::operator [] (const UInt32 &idx) {
  
  // input is valid?
	// assert (idx >= 0 && idx < _size);
  
	if(!_is_compact){
    // LCP array has not been compacted yet!
		return array[idx];
  }
  
	if(_p_array[idx] < 255){
    // Found in primary index
    return (UInt32) _p_array[idx];
  }

  
  // svnvish: BUGBUG
  // Do some caching here.
	
	// 	// Now search in secondary index as last resort
	// 	std::pair< const_itr, const_itr > p = equal_range(_beg, _end, idx);
	// 	return _val_array[std::distance(_beg, p.first)];
  
	if (++_cache == _end){
		_cache = _beg;
		_dist = 0;
	}else{
		_dist++;
	}
	
	UInt32 c_idx = *(_cache);
	
	if (c_idx == idx){
		return _val_array[_dist];
	}


  //	_cache = std::equal_range(_beg, _end, idx).first;
	_cache = std::lower_bound(_beg, _end, idx);
#ifdef _RWSTD_NO_CLASS_PARTIAL_SPEC
	_dist = 0;
	std::distance(_beg, _cache, _dist);
#else
	_dist = std::distance(_beg, _cache);
#endif		
  //std::cout << "here" << std::endl;

	// _cache = equal_range(_beg, _end, idx).first;
	// _dist = std::distance(_beg, _cache);
	
	return _val_array[_dist];


// 	if (c_idx > idx){
// 		_cache = equal_range(_beg, _cache, idx).first;
// 	}else{
// 		_cache = equal_range(_cache, _end, idx).first;
// 	}

// 	//_cache = p.first;
// 	_dist = std::distance(_beg, _cache);	
// 	return _val_array[_dist];
	
}


/**
 *  Dump array elements to output stream.
 *  
 *  \param os  - (IN) Output stream
 *  \param lcp - (IN) LCP object.
 */
std::ostream& 
operator << (std::ostream& os, LCP& lcp){
  
  for( UInt32 i = 0; i < lcp._size; i++ ){
    os << "lcp[ " << i << "]: " << lcp[i] << std::endl;
  }
	return os;
}
#endif
