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


// File    : sask/Code/StringKernel.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006
//           10 Aug 2006


#ifndef STRINGKERNEL_H
#define STRINGKERNEL_H


#include "datatype.h"
#include "errorcode.h"
#include "esa.h"
#include "isafactory.h"
#include "ilcpfactory.h"
#include "iweightfactory.h"

//#include "W_msufsort.h"
#include "wkasailcp.h"

#include "cweight.h"
#include "expdecayweight.h"
#include "brweight.h"
#include "kspectrumweight.h"



//' Types of substring weighting functions
enum WeightFunction{CONSTANT, EXPDECAY, KSPECTRUM, BOUNDRANGE};

using namespace std;

class StringKernel {
 
  
 public:
	/// Variables
	ESA				      *esa;
	I_WeightFactory	*weigher;
	Real            *val;  //' val array. Storing precomputed val(t) values.
	Real			      *lvs;  //' leaves array. Storing weights for leaves.
	

	/// Constructors
	StringKernel();

	//' Given contructed suffix array
	StringKernel(ESA *esa_, int weightfn, Real param, int verb=INFO);

	//' Given text, build suffix array for it
	StringKernel(const UInt32 &size, SYMBOL *text, int weightfn, Real param, int verb=INFO);


	/// Destructor
	virtual ~StringKernel();

	//' Methods

	/// Precompute the contribution of each intervals (or internal nodes)
	void PrecomputeVal();
	
	/// Compute Kernel matrix
	void Compute_K(SYMBOL *xprime, const UInt32 &xprime_len, Real &value);

	/// Set leaves array, lvs[]
	void Set_Lvs(const Real *leafWeight, const UInt32 *len, const UInt32 &m);

	/// Set leaves array as lvs[i]=i for i=0 to esa->length
	void Set_Lvs();

 private:

  int _verb;

	/// An iterative auxiliary function used in PrecomputeVal()
	void IterativeCompute(const UInt32 &left, const UInt32 &right);

};
#endif
