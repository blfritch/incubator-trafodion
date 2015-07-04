// **********************************************************************
// @@@ START COPYRIGHT @@@
//
// (C) Copyright 2015 Hewlett-Packard Development Company, L.P.
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
// **********************************************************************

#include "Platform.h"
#include "hs_faststats.h"
#include "hs_globals.h"
#include "BloomFilter.h"

/**
 * Estimates the number of distinct keys that will be stored in a Counting
 * Bloom Filter (CBF), based on a subsample of values. When the first rowset
 * has been read, the values of a given column in the rowset are hashed, and
 * the contents of the hash table are then used to determine frequency-of-
 * frequency and number of distinct values for the subsample. This information
 * is then used to estimate the number of distinct values for the full sample.
 *
 * The expected number of distinct values in an important parameter of the CBF,
 * and estimating it based on the first rowset (although this does not represent
 * a random subsample) allows us to set up the CBF without a separate preliminary
 * pass over the data in the sample table.
 *
 * @param numRows Number of rows available to look at.
 * @param ptr Pointer to the start of the data area for the column being processed.
 * @return The estimated number of distinct keys for the given column in the
 *         full sample file.
 */
template <class T, class E>
UInt32 FastStatsHist<T,E>::sizeCBF(Lng32 numRows, void* ptr)
{
  NAHashDictionary<T, UInt32> hashTable(fsHashFunc<T>, 7499, TRUE, STMTHEAP);
  T* dataPtr = (T*)ptr;
  UInt32* countPtr;
  UInt32 currValCount = 0;
  UInt32 maxCount = 0;
  UInt32 nulls = 0;
  short* nullIndic = group_->nullIndics;

  for (Lng32 i=0; i<numRows; i++, dataPtr++)
    {
      if (nullIndic && *nullIndic++ == -1)
        {
          nulls++;
          continue;
        }

      countPtr = hashTable.getFirstValue(dataPtr);

      // If the value was not found in the hash table, insert it with a count
      // of 1. Otherwise, increment the count of the found value.
      if (!countPtr )
        {
          hashTable.insert(dataPtr, new(STMTHEAP) UInt32(1));
          currValCount = 1;
        }
      else
        {
          (*countPtr)++;
          currValCount = *countPtr;
        }

      // currValCount now holds the updated number of occurrences of the value
      // used in this iteration of the loop. If it is greater than the previously
      // observed greatest count, store it as the new max.
      if (currValCount > maxCount )
        maxCount = currValCount;
    }

  // Iterate over keys in hash table, and for each key, use its frequency to
  // index into the frequency-of-frequencies array and increment the value at
  // that index. Also sum the frequencies.
  FrequencyCounts fc;
  NAHashDictionaryIterator<T, UInt32> hashIter(hashTable);
  T* valPtr;
  UInt32* freqPtr;
  Int64 totalFreq = 0;
  for (CollIndex i=0; i<hashIter.entries(); i++)
    {
      hashIter.getNext(valPtr, freqPtr);
      fc.increment(*freqPtr);
      totalFreq += *freqPtr;
      NADELETEBASIC(freqPtr, STMTHEAP);
    }

  // Now add frequency information for nulls.
  Float64 subSampleUec = hashIter.entries();
  if (nulls > 0)
    {
      subSampleUec++;      // add null as another distinct value
      fc.increment(nulls);
      totalFreq += nulls;
    }

  // Use the linear weighted combination estimator (combo of Jackknife and Schlosser
  // estimators) to get the estimated number of distinct values in the full sample.
  Float64 DshMax = CmpCommon::getDefaultNumeric(USTAT_DSHMAX);
  Float64 coeffOfVar;
  Float64 Duj;
  Float64 Dsh;
  UInt32 numKeys = lwcUecEstimate(subSampleUec,    // # distinct values in subsample
                                  numRows,         // size of subsample
                                  origSampleRows_, // size of full sample
                                  &fc, DshMax, coeffOfVar, Duj, Dsh);

  // Adjust estimate of number of distinct keys if high UEC.
  UInt32 geoNumKeys = (UInt32)sqrt(Duj * Dsh);
  if ((Float64)totalFreq / hashIter.entries() <= 1.5)
    numKeys = geoNumKeys;

  return numKeys;
}

/*
 * For the column this histogram represents, processes a rowset read from the
 * sample table, inserting each value into the appropriate CBF.
 *
 * @param numRows Number of rows in the rowset.
 */
template <class T, class E>
void FastStatsHist<T,E>::addRowset(Lng32 numRows)
{
  // Establish ptrs to the null indicator and data arrays bound to the column,
  // which are populated when the rowset is read.
  setUpValues<T>(group_, numRows);
  T* dataPtr = (T*)group_->data;
  Lng32 dataLen = group_->ISlength;
  short* nullIndic = group_->nullIndics;
  T localMin;
  T localMax;

  // If this is the initial rowset, use it as a subsample to determine n,
  // the number of distinct keys to be represented in the CBF.
  if (totalFreq_ == 0)
    {
      UInt32 numKeys = numRows;
      if (numKeys < origSampleRows_)
        numKeys = sizeCBF(numRows, dataPtr);

      // Now that we know the approximate number of distinct keys, the CBF can
      // be created. Note that the number of hash fns, probability of false positive,
      // and max non-overflow frequency are still specified by constants at this
      // point. False positive probability has a profound effect on the accuracy
      // of a CBF, so that parameter should probably come from a CQD instead.
      cbf_ = new(STMTHEAP) FastStatsCountingBloomFilter(STMTHEAP, 5, numKeys, .01, 255);
    }

  totalFreq_ += numRows;

  NABoolean firstNonNull = TRUE;
  for (Int64 i=0; i<numRows; i++)
    {
      if (nullIndic && *nullIndic++ == -1)
        nullCount_++;
      else
        {
          // Initial min and max for this rowset are set to the first non-null
          // value encountered. If this is the first rowset, also initialize
          // min_ and max_, the encoded overall min and max.
          if (firstNonNull)
            {
              localMin = *dataPtr;
              localMax = *dataPtr;
              if (totalFreq_ == 0)
                {
                  min_ = encode<T,E>(localMin);
                  max_ = encode<T,E>(localMax);
                }
              firstNonNull = FALSE;
            }
          else
            {
              // Update min or max if the current value is not within their range.
              if (localMin > *dataPtr)
                localMin = *dataPtr;
              else if (localMax < *dataPtr)
                localMax = *dataPtr;
            }

          cbf_->insert(getTrueDataPtr(dataPtr), getTrueDataLen(dataPtr, dataLen));
        }

      dataPtr++;
    }

  // Update the overall minimum and maximum encoded value if eclipsed by the local
  // min/max for this rowset.
  E encodedLocalMin = encode<T,E>(localMin);
  if (min_ > encodedLocalMin)
    min_ = encodedLocalMin;
  E encodedLocalMax = encode<T,E>(localMax);
  if (max_ < encodedLocalMax)
    max_ = encodedLocalMax;
}

/**
 * Generates the final histogram for a column, of the form utilized by the
 * optimizer.
 *
 * @param numEHIntervals Number of intervals in the faststats equi-height
 *                       histogram that we are using to generate the final
 *                       official histogram.
 * @param samplePercent The sampling rate used when the persistent sample table
 *                      was created, needed for the estimation of row count and
 *                      UEC in the underlying full table.
 */
template <class T, class E>
void FastStatsHist<T,E>::generateHSHistogram(Lng32 numEHIntervals, Float64 samplePercent)
{
  Int32 numEWIntervals = 4 * numEHIntervals;
  E width = (max_ - min_ + (numEWIntervals - 1)) / numEWIntervals; //@ZXbl
  if (width == 0) // range of values (e.g., _SALT_) may be less than # intervals
    width = 1;

  // Average height of intervals in equi-width histogram. This is passed to the
  // histogram ctor below when creating the equi-width histogram, and is used as
  // the initial number of elements in the NAArray underlying the interval (an
  // array of value/frequency pairs).
  Int32 keysPerEWInterval = cbf_->getAllKeys().entries() / numEWIntervals;

  FSHistogram<E> equiWidthHist(STMTHEAP, numEWIntervals, keysPerEWInterval);

  // Now compute the equi-width histogram.

  const NAList<simple_cbf_key>& keys = cbf_->getAllKeys();
  CollIndex intvlInx;
  UInt64 freq;
  E keyVal;

  // Iterate over each distinct value found and add it to the correct interval
  // of the equi-width histogram.
  for (CollIndex i=0; i<keys.entries(); i++ )
  {
     const simple_cbf_key& key = keys[i];

     // Look up the key in the CBF and find its frequency of occurrence.
     if (!cbf_->contain(key.getKey(), key.getKeyLen(), &freq))
       continue;  // why would the key not be found in CBF?

     recordKeyInfo<T>(key, freq, vcTotLen_);

     // compute the interval index for the key
     keyVal = encodeFromCbfKey<T,E>(key);
     intvlInx = (CollIndex)((keyVal - min_) / width);
     if (intvlInx == numEWIntervals)
       intvlInx--;

     if (intvlInx < 0)
       continue;  // shouldn't happen if min/max maintained correctly

     // Insert the encoded value and freq pair into the interval.
     KeyFreqPair<E> vf(keyVal, (UInt32)freq);
     equiWidthHist[intvlInx].append(vf);
  }

  // If a varchar, set the avg length.
  if (DFS2REC::isAnyVarChar(group_->ISdatatype))
    group_->avgVarCharSize = vcTotLen_ / (double)totalFreq_;

  // Now convert the equi-width histogram into equal height one.

  Float32 skRatio = 1.00;

  // Set the target interval height to the total frequency divided by the
  // desired number of intervals.
  Int32 height = totalFreq_ / numEHIntervals;

  // 3rd arg is initial size (# KeyFreqPairs) of each interval -- intervals of
  // an EH histogram do not consist of KeyFreqPairs.
  FSHistogram<E> equiHeightHist(STMTHEAP, numEHIntervals, 1);

  equiWidthHist.convertToEquiHeightHist(height, equiHeightHist);
  equiHeightHist.estimateRowsAndUecs(samplePercent, skRatio, group_);

  // Create the HSHistogram that is used to write info to the histogram tables.
  HSHistogram* hsHist = new (STMTHEAP) HSHistogram(equiHeightHist.entries());
  group_->groupHist = hsHist;

  HSDataBuffer boundary;
  T typedBoundary;

  // Set boundary for special interval 0, the minimum of all values encountered.
  // The boundary for all other intervals is the high value represented by that
  // interval.
  decode(min_, typedBoundary);
  setBufferValue(typedBoundary, group_, boundary);
  hsHist->setIntBoundary(0, boundary);
  freeDecodedValue(typedBoundary);

  // Set the boundary (max) value and most frequent value for each of the
  // regular intervals, which start at index 1 (the FSHistogram intervals start
  // at 0).
  for (CollIndex i=1; i<=equiHeightHist.entries(); i++)
    {
      FSInterval<E>& fsInterval = equiHeightHist[i-1];
      decode(fsInterval.getBoundary(), typedBoundary);
      setBufferValue(typedBoundary, group_, boundary);
      hsHist->setIntBoundary(i, boundary);
      freeDecodedValue(typedBoundary);

      // Pardon the reuse of boundary/typedBoundary variables for mfv.
      decode(fsInterval.getMfv(), typedBoundary);
      setBufferValue(typedBoundary, group_, boundary);
      hsHist->setIntMFVValue(i, boundary);
      freeDecodedValue(typedBoundary);

      hsHist->setIntRowCount(i, fsInterval.getRc());
      hsHist->setIntUec(i, fsInterval.getUec());
      hsHist->setIntMFVRowCount(i, fsInterval.getMfvc());
      hsHist->setIntMFV2RowCount(i, fsInterval.getMfvc2());
      hsHist->setIntSquareSum(i, fsInterval.getSquareCntSum());
      hsHist->setIntOrigUec(i, fsInterval.getSampleUec());
    }

  // Add a null interval at the end if there were any nulls observed. The function
  // called takes care of setting all needed values for the interval.
  if (nullCount_ > 0)
    hsHist->addNullInterval(nullCount_, 1);
}

/*
 * Scales up this interval's row count and UEC from the sampled values.
 * The estimated row count is derived by simple linear scaling, while the UEC
 * uses a linear weighted combination of the Jackknife and Schlosser estimators.
 *
 * @param sampleRate Percentage of rows included in the sample
 * @param skRatio Skew ratio.
 * @param[out] cov Coefficient of variation for this interval.
 */
template <class E>
void FSInterval<E>::estimateRowsAndUecs(Float64 sampleRate,
                                        Float32 skRatio,
                                        Float64& cov)
{
  if (sampleUec_ == 0)
    uec_ = rc_ = mfvc_ = mfvc2_ = 0;
  else if (sampleUec_ == 1)
    {
      uec_ = 1;
      rc_ = mfvc_ = freqCount_ / sampleRate;
      mfvc2_ = 0;
    }
  else if (sampleUec_ == freqCount_)
    {
      uec_ = rc_ = freqCount_ / sampleRate;
      mfvc_ = 1;
      mfvc2_ = (freqCount_ > 1 ? 1 : 0);
    }
  else if (sampleRate == 1.0)
    {
      uec_ = sampleUec_;
      rc_ = freqCount_;
      mfvc_ = sampleMfvc_;
      mfvc2_ = sampleMfvc2_;
    }
  else
    {
      // Before estimating rc/uec, remove skewed frequencies from consideration, and
      // reduce the sampled number of keys by the number having those frequencies.
      if (skRatio < 1.0)
        freqOfFreqs_.removeSkew(skRatio * freqCount_, sampleUec_);

      Float64 sampleUec = (Float64)sampleUec_;
      Float64 sampleRowCnt = (Float64)freqCount_;
      Float64 DshMax = CmpCommon::getDefaultNumeric(USTAT_DSHMAX);
      Float64 Duj = 0;
      Float64 Dsh = 0;
      Float64 estTotalRC = sampleRowCnt / sampleRate;
      Float64 uec = lwcUecEstimate(sampleUec, sampleRowCnt, estTotalRC, &freqOfFreqs_,
                                   DshMax, cov, Duj, Dsh);
      uec_ = (Int32)uec;
      rc_ = sampleRowCnt / sampleRate;
      mfvc_ = sampleMfvc_ / sampleRate;
      mfvc2_ = sampleMfvc2_ / sampleRate;
    }
}

/**
 * Converts this equi-width histogram to one with intervals of (approximately)
 * equal height, returning the equi-height histogram in \c equiHeightHistogram.
 *
 * @param targetHeight Target height of the intervals of the equi-height histogram.
 * @param equiHeightHistogram The equi-height histogram, empty on input.
 */
template <class E>
void FSHistogram<E>::convertToEquiHeightHist(Int32 targetHeight,
                                             FSHistogram& equiHeightHistogram)
{
  Int32 ehInx = 0; // index of current equi-height histogram interval being filled
  Int64 availableEHfreq = targetHeight;
  UInt64 remainingEWfreq = 0;

  for (CollIndex i=0; i<this->entries(); i++)  // for each EW interval
    {
      FSInterval<E>& ewInterval = (*this)[i];
      if (remainingEWfreq != 0)
        HS_ASSERT(remainingEWfreq == 0);
      remainingEWfreq = ewInterval.getFreqCount();
      CollIndex ewKeyInx = 0;

      while (ewKeyInx < ewInterval.entries())  // while there are key/freq pairs not
        {                                      //   yet assigned to an EH interval
          if (remainingEWfreq <= availableEHfreq)
            {
              // It fits; put all remaining keys for this EW interval into the
              // current interval of the equi-height (EH) histogram.
              FSInterval<E>& ehInterval = equiHeightHistogram[ehInx];
              availableEHfreq -= remainingEWfreq;
              remainingEWfreq = 0;
              while (ewKeyInx < ewInterval.entries())
                ehInterval.recordFrequencyInfo(ewInterval[ewKeyInx++]);
            }
          else
            {
              // If this is the first time key/frequency pairs from this EW
              // interval have been split between 2 intervals of the EH histogram,
              // sort the key/freq pairs by key, so we can assure that the
              // lower values go into the lower-numbered EH interval.
              if (!ewInterval.sorted())
                ewInterval.sort();

              while (availableEHfreq > 0 && ewKeyInx < ewInterval.entries())
                {
                  // The determination of whether to include a key in the current
                  // EH interval could be more sophisticated; currently it requires
                  // it not to exceed the target height at all, unless it is greater
                  // than the target height by itself, in which case it gets its
                  // own EH interval.
                  Int32 freq = ewInterval[ewKeyInx].freq;
                  if (freq <= availableEHfreq ||
                      (freq > targetHeight && availableEHfreq == targetHeight))
                    {
                      availableEHfreq -= freq;
                      remainingEWfreq -= freq;
                      equiHeightHistogram[ehInx].recordFrequencyInfo(ewInterval[ewKeyInx]);
                      ewKeyInx++;
                    }
                  else  // no room in ehInterval for this key's frequency
                    availableEHfreq = 0;
                }

              // If the current EH interval is full, move to the next one and
              // start assigning to it.
              if (availableEHfreq <= 0)
                {
                  ehInx++;
                  if (ehInx >= equiHeightHistogram.entries())
                    equiHeightHistogram.increaseBuckets(heap_, 10);
                  availableEHfreq = targetHeight;
                }
            } // else
        } //while
    } // for

  // Remove any unused intervals from the EH histogram, including the current
  // one if nothing was added to it.
  if (equiHeightHistogram[ehInx].getFreqCount() > 0)
    ehInx++;
  equiHeightHistogram.clearFrom(ehInx);
}

// Explicit instantiations of template member functions, so their definition
// can appear in this file instead of in .h file.
template UInt32 FastStatsHist<Int32,Int32>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<UInt32,UInt32>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<Int16,Int16>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<UInt16,UInt16>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<Int64,Int64>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<Float32,Float32>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<Float64,Float64>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<ISFixedChar,Float64>::sizeCBF(Lng32, void*);
template UInt32 FastStatsHist<ISVarChar,Float64>::sizeCBF(Lng32, void*);

template void FastStatsHist<Int32,Int32>::addRowset(Lng32 numRows);
template void FastStatsHist<UInt32,UInt32>::addRowset(Lng32 numRows);
template void FastStatsHist<Int16,Int16>::addRowset(Lng32 numRows);
template void FastStatsHist<UInt16,UInt16>::addRowset(Lng32 numRows);
template void FastStatsHist<Int64,Int64>::addRowset(Lng32 numRows);
template void FastStatsHist<Float32,Float32>::addRowset(Lng32 numRows);
template void FastStatsHist<Float64,Float64>::addRowset(Lng32 numRows);
template void FastStatsHist<ISFixedChar,Float64>::addRowset(Lng32 numRows);
template void FastStatsHist<ISVarChar,Float64>::addRowset(Lng32 numRows);

template void FastStatsHist<Int32,Int32>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<UInt32,UInt32>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<Int16,Int16>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<UInt16,UInt16>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<Int64,Int64>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<Float32,Float32>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<Float64,Float64>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<ISFixedChar,Float64>::generateHSHistogram(Lng32, Float64);
template void FastStatsHist<ISVarChar,Float64>::generateHSHistogram(Lng32, Float64);

template void FSInterval<Int32>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<UInt32>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<Int16>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<UInt16>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<Int64>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<Float32>::estimateRowsAndUecs(Float64, Float32, Float64&);
template void FSInterval<Float64>::estimateRowsAndUecs(Float64, Float32, Float64&);

