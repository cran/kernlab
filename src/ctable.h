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


// File    : sask/Code/ChildTable.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006


#ifndef CTABLE_H
#define CTABLE_H

#include <vector>
#include <iostream>

#include "datatype.h"
#include "errorcode.h"
#include "lcp.h"

// using namespace std;

/**
 *  ChildTable represents the parent-child relationship between
 *  the lcp-intervals of suffix array.
 *  Reference: AboKurOhl04
 */
class ChildTable : public std::vector<UInt32>
{
	
 private:
  // childtab needs lcptab to differentiate between up, down, and
	//  nextlIndex values.
	LCP& _lcptab;
  
 public:
	
	// Constructors
	ChildTable(const UInt32 &size, LCP& lcptab): std::vector<UInt32>(size), _lcptab(lcptab){}
  
	// Destructor
	virtual ~ChildTable() {}

  
	// Get first l-index of an l-[i..j] interval
	ErrorCode l_idx(const UInt32 &i, const UInt32 &j, UInt32 &idx);
  
	// .up field
	ErrorCode up(const UInt32 &idx, UInt32 &val);

	// .down field
	ErrorCode down(const UInt32 &idx, UInt32 &val);

	// .next field can be retrieved by accessing the array directly.

  friend std::ostream& operator << (std::ostream& os, const ChildTable& ct); 

};
#endif
