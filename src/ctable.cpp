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


// File    : sask/Code/ChildTable.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
#ifndef CTABLE_CPP
#define CTABLE_CPP

#include "ctable.h"
#include <cassert>

/**
 *  Return the value of idx-th "up" field of child table. 
 *   val = childtab[idx -1];
 *
 *  \param idx - (IN)  The index of child table.
 *  \param val - (OUT) The value of idx-th entry in child table's "up" field.
 */
ErrorCode 
ChildTable::up(const UInt32 &idx, UInt32 &val){

	if(idx == size()) {
		// Special case: To get the first 0-index
		val = (*this)[idx-1];
		return NOERROR;
	}

  // svnvish: BUGBUG
  // Do we need to this in production code?
	UInt32 lcp_idx = 0, lcp_prev_idx = 0;
  lcp_idx = _lcptab[idx];
	lcp_prev_idx = _lcptab[idx-1];

  assert(lcp_prev_idx > lcp_idx);
  val = (*this)[idx-1];

	return NOERROR;
}

/**
 *  Return the value of idx-th "down" field of child table.  Deprecated. 
 *    Instead use val = childtab[idx];
 *
 *  \param idx - (IN)  The index of child table.
 *  \param val - (OUT) The value of idx-th entry in child table's "down" field.
 */
ErrorCode 
ChildTable::down(const UInt32 &idx, UInt32 &val){
  
	// For a l-interval, l-[i..j], childtab[i].down == childtab[j+1].up
	// If l-[i..j] is last child-interval of its parent OR 0-[0..n], 
	//   childtab[i].nextlIndex == childtab[i].down

  // svnvish: BUGBUG
  // Do we need to this in production code?
// 	UInt32 lcp_idx = 0, lcp_nextidx = 0;
// 	lcp_nextidx = _lcptab[(*this)[idx]];
// 	lcp_idx = _lcptab[idx];
// 	assert(lcp_nextidx > lcp_idx);

	// childtab[i].down := childtab[i].nextlIndex
	val = (*this)[idx];
	
	return NOERROR;
}


/**
 *  Return the first l-index of a given l-[i..j] interval.
 *
 *  \param i   - (IN)  Left bound of l-[i..j]
 *  \param j   - (IN)  Right bound of l-[i..j]
 *  \param idx - (OUT) The first l-index.
 */

ErrorCode 
ChildTable::l_idx(const UInt32 &i, const UInt32 &j, UInt32 &idx){
  
	UInt32 up = (*this)[j]; 
  
	if(i < up && up <= j){
		idx = up;
  }else {
		idx = (*this)[i]; 
	}									
	return NOERROR;
}


/**
 *  Dump array elements to output stream
 *
 *  \param os - (IN) Output stream.
 *  \param ct - (IN) ChildTable object.
 */
std::ostream& 
operator << (std::ostream& os, const ChildTable& ct){
  
  for( UInt32 i = 0; i < ct.size(); i++ ){
    os << "ct[ " << i << "]: " << ct[i] << std::endl;
  }
	return os;
}

#endif
