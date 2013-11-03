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


// File    : sask/Code/ExpDecayWeight.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006

#ifndef EXPDECAYWEIGHT_CPP
#define EXPDECAYWEIGHT_CPP

#include <cmath>
#include <cassert>

#include "expdecayweight.h"

using namespace std;


/**
 *  Exponential Decay weight function.
 *  W(y,t) := (lambda^{-gamma} - lambda^{-tau}) / (lambda - 1)
 *
 *  \param floor_len - (IN) Length of floor interval of matched substring.
 *                            (cf. gamma in VisSmo02).
 *  \param x_len     - (IN) Length of the matched substring.
 *                            (cf. tau in visSmo02).
 *  \param weight    - (OUT) The weight value.
 *
 */

ErrorCode
ExpDecayWeight::ComputeWeight(const UInt32 &floor_len, const UInt32 &x_len,	Real &weight)

// ErrorCode
// ExpDecayWeight::ComputeWeight(const Real &floor_len, const Real &x_len,	Real &weight)
{
	//' Input validation
	assert(x_len >= floor_len);
	
	//' x_len == floor_len when the substring found ends on an interval.
	if(floor_len == x_len) {
		//' substring ended on an interval, so, get the val from val[]
		weight = 0.0;
	}
	else {
		//weight = (pow(-(floor_len-1), lambda) - pow(-x_len, lambda)) / (1-lambda);
    //weight = (pow(lambda,((Real)floor_len)) - pow(lambda, (Real)x_len+1)) / (1-lambda);
    //     double a=floor_len*-1.0;
    //     double b=x_len*-1.0;
    //     weight = (pow(lambda,a) - pow(lambda, b)) / (lambda-1);
    weight = (pow(lambda,Real(-1.0*floor_len)) - pow(lambda, Real(-1.0*x_len))) / (lambda-1);
  }
  
//   std::cout << "floor_len : " << floor_len
//             << "   x_len : " << x_len
//             << "   pow1 : " << pow(lambda,-((Real)floor_len))
//             << "   pow2 : " << pow(lambda,-(Real)x_len)
//             << "   weight : " << weight << std::endl;

	return NOERROR;
}

#endif
