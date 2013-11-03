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


// File    : sask/Code/StringKernel.cpp
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           12 Jul 2006
//           10 Aug 2006
//           11 Oct 2006


#ifndef STRINGKERNEL_CPP
#define STRINGKERNEL_CPP

#include <queue>
#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <memory>
#include <cstring>
#include <cstdio>

#include "stringkernel.h"

StringKernel::StringKernel(): esa(0), weigher(0), val(0), lvs(0)
{}


/**
 *  Construct string kernel given constructed enhanced suffix array.
 *
 *  \param esa_ - ESA instance.
 */
StringKernel::StringKernel(ESA *esa_, int weightfn, Real param, int verb): esa(esa_), val(new Real[esa_->size + 1]), lvs(0), _verb(verb)
{

    switch (weightfn)
    {
    case CONSTANT:
        weigher = new ConstantWeight();
        break;
    case EXPDECAY:
        weigher = new ExpDecayWeight(param);
        break;
    case KSPECTRUM:
        weigher = new KSpectrumWeight(param);
        break;
    case BOUNDRANGE:
        weigher = new BoundedRangeWeight(param);
        break;

    default:
      int nothing = 0;

    }
}

/**
 *  Construct string kernel when given only text and its length.
 *
 * \param text         - (IN) The text which SuffixArray and StringKernel correspond to.
 * \param text_length  - (IN) The length of #_text#.
 * \param verb         - (IN) Verbosity level.
 */
StringKernel::StringKernel(const UInt32 &size, SYMBOL *text, int weightfn, Real param, int verb): lvs(0), _verb(verb)
{

    // Build ESA.
    esa = new ESA(size, text, verb);

    // Allocate memory space for #val#
    val = new Real[esa->size + 1];

    // Instantiate weigher.
    switch (weightfn)
    {
    case CONSTANT:
        weigher = new ConstantWeight();
        break;
    case EXPDECAY:
        weigher = new ExpDecayWeight(param);
        break;
    case KSPECTRUM:
        weigher = new KSpectrumWeight(param);
        break;
    case BOUNDRANGE:
        weigher = new BoundedRangeWeight(param);
        break;

    default:
      int nothing = 0;
    }
}


/**
 *  StringKernel destructor.
 *
 */
StringKernel::~StringKernel()
{
    //' Delete objects and release allocated memory space.
    if (esa)
    {
        delete esa;
        esa = 0;
    }
    if (val)
    {
        delete [] val;
        val = 0;
    }
    if (lvs)
    {
        delete [] lvs;
        lvs = 0;
    }
    if (weigher)
    {
        delete weigher;
        weigher = 0;
    }
}




/**
 *  An Iterative auxiliary function used in PrecomputeVal().
 *
 *  Note: Every lcp-interval can be represented by its first l-index.
 *          Hence, 'val' is stored in val[] at the index := first l-index.
 *
 *  Pre: val[] is initialised to 0.
 *
 *  @param left Left bound of current interval
 *  @param right Right bound of current interval
 */
void
StringKernel::IterativeCompute(const UInt32 &left, const UInt32 &right)
{
    //std::cout << "In IterativeCompute() " << std::endl;

    //' Variables
    queue<pair<UInt32, UInt32> > q;
    vector<pair<UInt32, UInt32> > childlist;
    pair<UInt32, UInt32> p;
    UInt32 lb = 0;
    UInt32 rb = 0;
    UInt32 floor_len = 0;
    UInt32 x_len = 0;
    Real cur_val = 0.0;
    Real edge_weight = 0.0;


    //' Step 1: At root, 0-[0..size-1]. Store all non-single child-intervals onto #q#.
    lb = left;   //' Should be equal to 0.
    rb = right;  //' Should be equal to size-1.
    esa->GetChildIntervals(lb, rb, childlist);

    for (UInt32 jj = 0; jj < childlist.size(); jj++)
        q.push(childlist[jj]);


    //' Step 2: Do breadth-first traversal. For every interval, compute val and add
    //'           it to all its non-singleton child-intervals' val-entries in val[].
    //'         Start with child-interval [i..j] of 0-[0..size-1].
    //'         assert(j != size-1)
    while (!q.empty())
    {
        //' Step 2.1: Get an interval from queue, #q#.
        p = q.front();
        q.pop();

        //' step 2.2: Get the lcp of floor interval.
        UInt32 a = 0, b = 0;

        a = esa->lcptab[p.first];
        //svnvish: BUGBUG
        // Glorious hack. We have to remove it later.
        // This gives the lcp of parent interval
        if (p.second < esa->size - 1)
        {
            b = esa->lcptab[p.second + 1];
        }
        else
        {
            b = 0;
        }
        floor_len = (a > b) ? a : b;


        //' Step 2.3: Get the lcp of current interval.
        esa->GetLcp(p.first, p.second, x_len);


        //' Step 2.4: Compute val of current interval.
        weigher->ComputeWeight(floor_len, x_len, edge_weight);
        cur_val = edge_weight * (lvs[p.second + 1] - lvs[p.first]);


        //' Step 2.5: Add #cur_val# to val[].
        UInt32 firstlIndex1 = 0;
        esa->childtab.l_idx(p.first, p.second, firstlIndex1);
        val[firstlIndex1] += cur_val;

        //  std::cout << "p.first:"<<p.first
        //               << "   p.second:"<<p.second
        //               << "   firstlIndex1 : " << firstlIndex1
        //               << "   floor_len:"<<floor_len
        //               << "   x_len:"<<x_len<< std::endl;

        //' Step 2.6: Get all child-intervals of this interval.
        childlist.clear();
        esa->GetChildIntervals(p.first, p.second, childlist);


        //' Step 2.7: (a) Add #cur_val# to child-intervals' val-entries in val[].
        //'           (b) Push child-interval onto #q#.
        for (UInt32 kk = 0; kk < childlist.size(); kk++)
        {
            //' (a)
            UInt32 firstlIndex2 = 0;
            pair<UInt32, UInt32> tmp_p = childlist[kk];

            if (esa->text[esa->suftab[tmp_p.first]] == SENTINEL)
                continue;

            esa->childtab.l_idx(tmp_p.first, tmp_p.second, firstlIndex2);

            // assert( val[firstlIndex2] == 0 );
            val[firstlIndex2] = val[firstlIndex1]; // cur_val;

            //' (b)
            q.push(make_pair(tmp_p.first, tmp_p.second));
        }
    }

    //std::cout << "Out IterativeCompute() " << std::endl;
}



/**
 *  Precomputation of val(t) of string kernel. 
 *  Observation :Every internal node of a suffix tree can be represented by at
 *                 least one index of the corresponding lcp array. So, the val
 *                 of a node is stored in val[] at the index corresponding to that of 
 *                 the fist representative lcp value in lcp[].
 */
void
StringKernel::PrecomputeVal()
{
    //' Memory space requirement check.
    assert(val != 0);


    //' Initialise all val entries to zero!
    memset(val, 0, sizeof(Real)*esa->size + 1);


    //' Start iterative precomputation of val[]
    IterativeCompute(0, esa->size - 1);
}


/**
 *  Compute k(text,x) by performing Chang and Lawler's matching statistics collection 
 *    algorithm on the enhanced suffix array.
 *
 *  \param x     - (IN) The input string which is to be evaluated together with
 *                        the text in esa.
 *  \param x_len - (IN) The length of #x#.
 *  \param value - (IN) The value of k(x,x').
 */
void
StringKernel::Compute_K(SYMBOL *x, const UInt32 &x_len, Real &value)
{
    //' Variables
    UInt32 floor_i = 0;
    UInt32 floor_j = 0;
    UInt32 i = 0;
    UInt32 j = 0;
    UInt32 lb = 0;
    UInt32 rb = 0;
    UInt32 matched_len = 0;
    UInt32 offset = 0;
    UInt32 floor_len = 0;
    UInt32 firstlIndex = 0;
    Real edge_weight = 0.0;


    //' Initialisation
    value = 0.0;
    lb = 0;
    rb = esa->size - 1;


    //' for each suffix, xprime[k..xprime_len-1], find longest match in text
    for (UInt32 k = 0; k < x_len; k++)
    {

        //' Step 1: Matching
        esa->ExactSuffixMatch(lb, rb, offset, &x[k], x_len - k, i, j, matched_len,
                              floor_i, floor_j, floor_len);


        //' Step 2: Get suffix link for [floor_i..floor_j]
        esa->GetSuflink(floor_i, floor_j, lb, rb);
        assert((floor_j - floor_i) <= (rb - lb));  //' Range check


        //' Step 3: Compute contribution of this matched substring
        esa->childtab.l_idx(floor_i, floor_j, firstlIndex);
        assert(firstlIndex > floor_i && firstlIndex <= floor_j);
        assert(floor_len <= matched_len);


        weigher->ComputeWeight(floor_len, matched_len, edge_weight);
        value += val[firstlIndex] + edge_weight * (lvs[j + 1] - lvs[i]);

        //    std::cout << "i:"<<floor_i
        //               << "   j:"<<floor_j
        //               << "   firstlIndex : " << firstlIndex
        //               << "   val[] : " << val[firstlIndex]
        //               << "   edge_wt : " << edge_weight
        //               << "   value :" << value << std::endl;

        //' Step 4: Prepare for next iteration. This tells how many
        //'         char ExactSuffixMatch() can jump before matching
        //'         pattern to text.
        offset = (matched_len) ? matched_len - 1 : 0;
    }
}



/**
 *  Assign leaf weights. 
 *
 *  Assumption: SENTINEL is only used in delimiting the concatenated strings.
 *                No string should contain SENTINEL.
 *
 *  \param leafWeight - (IN) The leaf weights (or \alpha in SVM).
 *  \param len        - (IN) The array containing the length of each string in the Master String.
 *                             (Master String is the concatenation of SENTINEL-delimited strings.)
 *  \param m          - (IN) Size of the array #leafWeight# (equal to the number of datapoints).
 *
 */
void
StringKernel::Set_Lvs(const Real *leafWeight, const UInt32 *len, const UInt32 &m)
{
    //' Clean up previous lvs, if any.
    if (lvs)
    {
        delete lvs;
        lvs = 0;
    }

    //' Variables
    UInt32 pos = 0;


    //' Let n denotes the length of Master String, and
    //'     m denotes the number of strings in Master String.

    //' n := sum{|string_i|} + m, where
    //'   m is the number of delimiter (i.e. '\n') added.


    //' Create a cumulative array of #len.
    UInt32 *clen = new (nothrow) UInt32[m];


    //' clen[] is a cumulative array of len[]
    partial_sum(len, len + m, clen);
    assert(clen[m - 1] == esa->size);


    //' Allocate memory space for lvs[]
    lvs = new (nothrow) Real[esa->size + 1];
    assert(lvs);


    //' Assign leaf weight to lvs element according to its position in text.
    for (UInt32 j = 0; j < esa->size; j++)
    {
        pos = esa->suftab[j];
        UInt32 *p = upper_bound(clen, clen + m, pos);   //' O(log n)
        lvs[j + 1] = leafWeight[p - clen];
    }


    //' Compute cumulative lvs[]. To be used in matching statistics computation later.
    lvs[0] = 0.0;
    partial_sum(lvs, lvs + esa->size + 1, lvs);

    //chteo: [101006]
    delete [] clen;
    clen = 0;
}



/**
 *  Set lvs[i] = i, for i = 0 to esa->size
 *  Memory space for lvs[] will be allocated.
 */
void
StringKernel::Set_Lvs()
{
    //' Clean up previous lvs, if any.
    if (lvs)
    {
        delete lvs;
        lvs = 0;
    }

    //' Allocate memory space for lvs[]
    lvs = new (nothrow) Real[esa->size + 1];

    //' Check if memory correctly allocated.
    assert(lvs != 0);

    //' Range := [0..esa->size]
    UInt32 localsize = esa->size;
    for (UInt32 i = 0; i <= localsize; i++)
        lvs[i] = i;
}

#endif

extern "C" {

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

  SEXP stringtv(SEXP rtext, // text document 
		SEXP ltext, // list or vector of text documents to compute kvalues against 
		SEXP nltext, // number of text documents in ltext
		SEXP vnchar, // number of characters in text
		SEXP vnlchar, // characters per document in ltext
		SEXP stype, // type of kernel 
		SEXP param) // parameter for kernel
  {
    // R interface for text and list of text computation. Should return a vector of computed kernel values.
    // Construct ESASK
    UInt32 text_size = *INTEGER(vnchar);
    int number_ltext = *INTEGER(nltext);
    int *ltext_size =  (int *) malloc (sizeof(int) * number_ltext);
    memcpy(ltext_size, INTEGER(vnlchar), number_ltext*sizeof(int));
    int weightfn = *INTEGER(stype);
    const char *text = CHAR(STRING_ELT(rtext,0)); 
    Real kparam = *REAL(param); 
    double kVal;
    SEXP alpha;
    
    PROTECT(alpha = allocVector(REALSXP, number_ltext));
    
    // Check if stringlength reported from R is correct  
    if(strlen(text)!= text_size)
      text_size= strlen(text);

    StringKernel sk(text_size, (SYMBOL*)text, (weightfn - 1), kparam, 0);
    sk.Set_Lvs();
    sk.PrecomputeVal();
    
    for (int i=0; i<number_ltext; i++)
      {
	
	if(isList(ltext)){
	  const char *pattern = CHAR(VECTOR_ELT(ltext, i));
	  // Check if stringlength reported from R is correct
	  if(strlen(pattern)!=ltext_size[i])
	    ltext_size[i]= strlen(pattern);
	  sk.Compute_K((SYMBOL*)pattern, ltext_size[i], kVal);
	}
	else{
	  const char *pattern = CHAR(STRING_ELT(ltext, i));
	  // Check if stringlength reported from R is correct
	  if(strlen(pattern)!=ltext_size[i])
	    ltext_size[i]= strlen(pattern);
	  sk.Compute_K((SYMBOL*)pattern, ltext_size[i], kVal);
	}
	REAL(alpha)[i] = kVal;
      }

    free(ltext_size);
    UNPROTECT(1); 
    return alpha;

  }
}
