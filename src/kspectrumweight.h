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


// File    : sask/Code/KSpectrumWeight.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006

#ifndef KSPECTRUMWEIGHT_H
#define KSPECTRUMWEIGHT_H

#include "datatype.h"
#include "errorcode.h"
#include "iweightfactory.h"
#include <iostream>

//' K-spectrum weight class
class KSpectrumWeight : public I_WeightFactory
{

   Real k;
	
public:

	/// Constructor
	KSpectrumWeight(const Real & k_=5.0):k(k_) {}

	/// Destructor
	virtual ~KSpectrumWeight(){}

	/// Compute weight
	ErrorCode ComputeWeight(const UInt32 &floor_len, const UInt32 &x_len, Real &weight);
};
#endif
