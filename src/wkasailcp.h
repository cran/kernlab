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


// File    : sask/Code/W_kasai_lcp.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006


#ifndef W_KASAI_LCP_H
#define W_KASAI_LCP_H

#include "datatype.h"
#include "errorcode.h"
#include "ilcpfactory.h"
#include "lcp.h"


/**
 * Kasai et al's LCP array computation algorithm is
 * is slightly faster than Manzini's algorithm. However,
 * it needs inverse suffix array which costs extra memory.
 */
class W_kasai_lcp : public I_LCPFactory
{

 public:

	/// Constructor
	W_kasai_lcp(){}

	/// Desctructor
	virtual ~W_kasai_lcp(){}

	/// Compute LCP array.
	ErrorCode	ComputeLCP(const SYMBOL *text, const UInt32 &len, 
											 const UInt32 *sa, LCP& lcp);

};
#endif
