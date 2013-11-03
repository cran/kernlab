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


// File    : sask/Code/ConstantWeight.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006
//           12 Oct 2006


#ifndef CWEIGHT_CPP
#define CWEIGHT_CPP

#include "cweight.h"
#include <cassert>

/**
 *  Constant weight function. Computes number of common substrings. Every 
 *    matched substring is of same weight (i.e. 1)
 *  W(y,t) := tau - gamma
 *
 *  \param floor_len - (IN) Length of floor interval of matched substring.
 *                            (cf. gamma in VisSmo02).
 *  \param x_len     - (IN) Length of the matched substring.
 *                            (cf. tau in visSmo02).
 *  \param weight    - (OUT) The weight value.
 *
 */
ErrorCode
ConstantWeight::ComputeWeight(const UInt32 &floor_len, const UInt32 &x_len,	Real &weight)
{
	//' Input validation
	assert(x_len >= floor_len);
		
	//' x_len == floor_len when the substring found ends on an interval.

	weight = (x_len - floor_len);

//   std::cout << "floor_len : " << floor_len
//             << "   x_len : " << x_len
//             << "   weight : " << weight << std::endl;

	return NOERROR;
}

#endif
