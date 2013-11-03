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
 * The Original Code is the MSufSort suffix sorting algorithm (Version 2.2).
 *
 * The Initial Developer of the Original Code is
 * Michael A. Maniscalco
 * Portions created by the Initial Developer are Copyright (C) 2006
 * the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 *
 *   Michael A. Maniscalco
 *
 * ***** END LICENSE BLOCK ***** */

#ifndef MSUFSORT_INDUCTION_SORTING_H
#define MSUFSORT_INDUCTION_SORTING_H

#include "introsort.h"


class InductionSortObject
{
public:
	InductionSortObject(unsigned int inductionPosition = 0, unsigned int inductionValue = 0, unsigned int suffixIndex = 0);

	bool operator <= (InductionSortObject & object);

	bool operator == (InductionSortObject & object);

	InductionSortObject& operator = (InductionSortObject & object);

	bool operator >= (InductionSortObject & object);

	bool operator > (InductionSortObject & object);

	bool operator < (InductionSortObject & object);

	unsigned int	m_sortValue[2];
};


inline bool InductionSortObject::operator <= (InductionSortObject & object)
{
	if (m_sortValue[0] < object.m_sortValue[0])
		return true;
	else
		if (m_sortValue[0] == object.m_sortValue[0])
			return (m_sortValue[1] <= object.m_sortValue[1]);
	return false;
}



inline bool InductionSortObject::operator == (InductionSortObject & object)
{
	return ((m_sortValue[0] == object.m_sortValue[0]) && (m_sortValue[1] == object.m_sortValue[1]));
}



inline bool InductionSortObject::operator >= (InductionSortObject & object)
{
	if (m_sortValue[0] > object.m_sortValue[0])
		return true;
	else
		if (m_sortValue[0] == object.m_sortValue[0])
			return (m_sortValue[1] >= object.m_sortValue[1]);
	return false;
}



inline InductionSortObject & InductionSortObject::operator = (InductionSortObject & object)
{
	m_sortValue[0] = object.m_sortValue[0];
	m_sortValue[1] = object.m_sortValue[1];
	return *this;
}




inline bool InductionSortObject::operator > (InductionSortObject & object)
{
	if (m_sortValue[0] > object.m_sortValue[0])
		return true;
	else
		if (m_sortValue[0] == object.m_sortValue[0])
			return (m_sortValue[1] > object.m_sortValue[1]);
	return false;
}



inline bool InductionSortObject::operator < (InductionSortObject & object)
{
	if (m_sortValue[0] < object.m_sortValue[0])
		return true;
	else
		if (m_sortValue[0] == object.m_sortValue[0])
			return (m_sortValue[1] < object.m_sortValue[1]);
	return false;
}




#endif
