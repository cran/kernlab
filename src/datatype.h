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


// File    : sask/Code/DataType.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           11 Oct 2006


#ifndef DATATYPE_H
#define DATATYPE_H
	
//  #define UInt32  unsigned int
//  #define UInt64  unsigned long long
//  #define Byte1   unsigned char
//  #define Byte2   unsigned short
//  #define Real    double

typedef unsigned int UInt32;

// Seems that even using __extension__ g++ 4.6 will complain that
// ISO C++ 1998 does not support 'long long' ...
/*
#if defined __GNUC__ && __GNUC__ >= 2
__extension__ typedef unsigned long long UInt64;
#else
typedef unsigned long long UInt64;
#endif
*/

#include <stdint.h>
typedef uint64_t UInt64;

typedef unsigned char Byte1;
typedef unsigned short Byte2;
typedef double Real;

//  #define SENTINEL  '\n'
//  #define SENTINEL2 '\0'

const char SENTINEL  =  '\n';
const char SENTINEL2 =  '\0';

#ifndef UNICODE
// #  define SYMBOL  Byte1
  typedef Byte1 SYMBOL;
#else
// #  define SYMBOL  Byte2
  typedef Byte2 SYMBOL;
#endif

#endif
