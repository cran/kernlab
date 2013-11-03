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


// File    : sask/Code/KSpectrumWeight.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006

#ifndef KSPECTRUMWEIGHT_CPP
#define KSPECTRUMWEIGHT_CPP

#include "kspectrumweight.h"
#include <cassert>



/**
 *  K-spectrum weight function. Compute number of common (exactly) k character substring.
 *  
 *  \param floor_len - (IN) Length of floor interval of matched substring. (cf. gamma in VisSmo02).
 *  \param x_len     - (IN) Length of the matched substring. (cf. tau in VisSmo02).
 *  \param weight    - (OUT) The weight value.
 *
 */
ErrorCode
KSpectrumWeight::ComputeWeight(const UInt32 &floor_len, const UInt32 &x_len, Real &weight)
{
	//' Input validation
	assert(x_len >= floor_len);
		
	//' x_len == floor_len when the substring found ends on an interval.
	
 
	weight = 0.0;

	if(floor_len < k && x_len >= k) 
		weight = 1.0;
	
//   std::cout << "floor_len : " << floor_len
//             << "   x_len : " << x_len 
//             << "   weight : " << weight << std::endl;

	return NOERROR;
}

#endif


//' Question: Why return only 0 or 1?
//' Answer  : In k-spectrum method, any length of matched substring other than k 
//'           does not play a significant role in the string kernel. So, returning 1
//'           means that the substring weight equals to # of suffix in the current interval.
//'           When 0 is returned, it means that substring weight equals to the floor 
//'           interval entry in val[]. (See the definition of substring weight in 
//'           StringKernel.cpp)

//' Question: Why is the following a correct implementation of k-spectrum ?
//' Answer  : [Val precomputation phase] Every Interval with lcp < k has val := 0.
//'           For intervals with (lcp==k) or (lcp>k but floor_lcp<k), they have 
//'           same val := # of substring in common.
//'           [Kernel evaluation phase] When "if" statement is FALSE, the substring 
//'           weight is equal to its floor interval entry in val[] (see explanation above).
//'           When "if" statement is TRUE, it means that the x_len >= k but floor interval
//'           has val := 0 (floor_lcp < k). Hence, returning weight:=1 will make substring 
//'           weight equals to the size of the immediate ceil interval (# of substring in common).
