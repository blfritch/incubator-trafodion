/* -*-C++-*-
****************************************************************************
*
* File:         HeapID.cpp
*
* Description:  Assign HeapID for a defined heap.
*
* Created:      3/1/99
* Language:     C++
*
*
// @@@ START COPYRIGHT @@@
//
// (C) Copyright 1999-2014 Hewlett-Packard Development Company, L.P.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
// @@@ END COPYRIGHT @@@
*
****************************************************************************
*/

#include "HeapLog.h"

#ifdef NA_DEBUG_HEAPLOG
// -----------------------------------------------------------------------
// Constructor and destructor.
// -----------------------------------------------------------------------
HeapID::HeapID()
{
  // Assign a unique ID.
  heapNum = HeapLogRoot::assignHeapNum();
}

HeapID::~HeapID()
{
  if (heapNum >= 0 &&
      HeapLogRoot::log != NULL)
    HeapLogRoot::deleteLogSegment(heapNum, TRUE);
}

#endif
