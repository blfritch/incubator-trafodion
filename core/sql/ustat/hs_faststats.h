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

#include "BloomFilter.h"
#include "dfs2rec.h"
#include "hs_globals.h"
#include "hs_const.h"
#include "hs_log.h"
#include "CharType.h"
#include <iostream>

struct HSColGroupStruct;

template <class E>
struct KeyFreqPair
{
  KeyFreqPair(E value = 0, UInt32 frequency = 0)
    : val(value), freq(frequency)
    {}

  ~KeyFreqPair() {};

  NABoolean operator<=(const KeyFreqPair& other)
  {
     if (val < other.val)
       return TRUE;

     if (val == other.val)
       return (freq <= other.freq);

     return FALSE;
  }

  NABoolean operator<(const KeyFreqPair& other)
  {
     if (val < other.val)
       return TRUE;

     if (val == other.val)
       return (freq < other.freq);

     return FALSE;
  }

  NABoolean operator==(const KeyFreqPair& other)
  {
    return (val  == other.val &&
            freq == other.getFreq());
  }

  KeyFreqPair& operator=(const KeyFreqPair& other)
  {
    val  = other.val;
    freq = other.freq;
    return *this;
  }

  E val;
  UInt32 freq;
};

/**
 * The purpose of this class is to allow a generic (independent of the type it
 * is instantiated with) pointer to FastStatsHist.
 */
class AbstractFastStatsHist : public NABasicObject
{
public:
    AbstractFastStatsHist() {}
  virtual ~AbstractFastStatsHist() {}

  virtual UInt32 sizeCBF(Lng32 numRows, void* dataPtr) = 0;
  virtual void addRowset(Lng32 numRows) = 0;
  virtual void generateHSHistogram(Lng32 numIntervals, Float64 samplePercent) = 0;
  virtual Int64 getTotalFrequency() const = 0;
};

template <class T, class E>
class FastStatsHist : public AbstractFastStatsHist
{
public:
  FastStatsHist(HSColGroupStruct* group, Int64 sampleRows)
    : group_(group),
      cbf_(NULL),
      nullCount_(0),
      totalFreq_(0),
      vcTotLen_(0),
      origSampleRows_(sampleRows)
    {}

  ~FastStatsHist()
  {
    delete cbf_;
  }

  virtual UInt32 sizeCBF(Lng32 numRows, void* ptr);
  virtual void addRowset(Lng32 numRows);
  virtual void generateHSHistogram(Lng32 numIntervals, Float64 samplePercent);
  virtual Int64 getTotalFrequency() const
  {
    return totalFreq_;
  }

private:
  HSColGroupStruct* group_;
  FastStatsCountingBloomFilter* cbf_;
  Int64 nullCount_;
  Int64 totalFreq_;
  E min_, max_;
  Int64 vcTotLen_;  // Total len of all varchars (if varchar col)
  Int64 origSampleRows_;
};

//@ZXbl: There are substantial differences between an EW interval and an EH interval.
//       They should probably be split into 2 subclasses of an abstract parent.
/**
 * Class representing an interval of a histogram. It is derived from an NAArray
 * of EncodeValueFreqPair objects, and so consists of some number of distinct
 * values, each paired with its frequency of occurrence.
 */
template <class E>
class FSInterval : public NAArray<KeyFreqPair<E>>
{
public:
  /**
   * Constructs an interval with the given row count and uec (0 by default).
   *
   * @param heap Heap to use in call to NAArray ctor.
   * @param count Initial number of elements in the NAArray.
   * @param rc Initial row count of the interval.
   * @param uec Initial UEC of the interval.
   */
  FSInterval(NAHeap* heap = NULL, Int32 count = 0, Int32 rc = 0, Int32 uec = 0)
   : NAArray<KeyFreqPair<E>>(heap, count),
     sampleUec_(0), freqCount_(0), squareCntSum_(0.0),
     mfvc_(0), mfvc2_(0), sampleMfvc_(0), sampleMfvc2_(0), mfv_(0),
     sorted_(FALSE), uec_(uec), rc_(rc), boundaryHasBeenSet_(FALSE)
   {}

  ~FSInterval()
  {}

  Int32 getFreqCount() const { return freqCount_; }

  // Return the scaled-up row count and # distinct values
  Int32 getRc() const { return rc_; }
  Int32 getUec() const { return uec_; }

  // Return the row count and # distinct values from the sample.
  Int32 getSampleRc() const { return freqCount_; }
  Int32 getSampleUec() const { return sampleUec_; }

  Float64 getSquareCntSum() const { return squareCntSum_; }
  NABoolean sorted() const { return sorted_; }
  E getBoundary() const { return boundary_; }
  E getMfv() const { return mfv_; }
  UInt32 getMfvc() const { return mfvc_; }
  UInt32 getMfvc2() const { return mfvc2_; }

  // To sort KeyFreqPairs, we operate directly on the array underlying the
  // NAArray object, which can be done because we never remove any elements
  // (which would produce "holes" in the array).
  void sort()
  {
    KeyFreqPair<E>* arr = this->getArrIfIntact();
    HS_ASSERT(arr != NULL);
    qsort(arr, this->entries(), sizeof(KeyFreqPair<E>), compareOp);
    sorted_ = TRUE;
  }

  // Comparison function for qsort of KeyFreqPairs. Since all keys are guaranteed
  // to be distinct, we can skip equality comparison.
  static int compareOp(const void* v1, const void* v2)
  {
    return *(KeyFreqPair<E>*)v1 < *(KeyFreqPair<E>*)v2 ? -1 : 1;
  }

  void append(KeyFreqPair<E>& x)
  {
    insertAt(this->entries(), x);
    freqCount_ += x.freq;
  }

  void recordFrequencyInfo(KeyFreqPair<E>& kfp)
  {
    UInt32 freq = kfp.freq;
    sampleUec_++;
    freqCount_ += freq;
    squareCntSum_ += freq * freq;
    freqOfFreqs_.increment(freq);
    if (freq > sampleMfvc_)
      {
        sampleMfvc2_ = sampleMfvc_;
        sampleMfvc_ = freq;
        mfv_ = kfp.val;
      }
    else if (freq > sampleMfvc2_)
      sampleMfvc2_ = freq;

    // Boundary is highest value in the interval.
    if (kfp.val > boundary_ || !boundaryHasBeenSet_)
      {
        boundary_ = kfp.val;
        boundaryHasBeenSet_ = TRUE;
      }
  }

  void estimateRowsAndUecs(Float64 sampleRate, Float32 skRatio, Float64& cov);

private:
  void swap(UInt32 left, UInt32 right)
  {
    KeyFreqPair<E> tmp = this->at(left);
    this->at(left) = this->at(right);
    this->at(right) = tmp;
  }

  FrequencyCounts freqOfFreqs_;
  Int32 sampleUec_;      // # distinct values found
  Int32 freqCount_;      // total frequency of all keys of this interval
  Float64 squareCntSum_; // sum of squared frequencies (to compute std dev)
  UInt32 sampleMfvc_;    // Most frequent value count from sampled data
  UInt32 sampleMfvc2_;   // 2nd most frequent value count from sampled data
  E mfv_;                // Most frequently occurring value
  NABoolean sorted_;     // true if KeyFreqPairs have been sorted

  Int32 rc_;             // estimated number of rows for the interval
  Int32 uec_;            // estimated number of distinct values
  UInt32 mfvc_;          // Estimated most frequent value count
  UInt32 mfvc2_;         // Estimated 2nd most frequent value count
  NABoolean boundaryHasBeenSet_;
  E boundary_;           // upper boundary value for the interval
};


/**
 * A histogram, consisting of a list of FSInterval objects.
 */
template <class E>
class FSHistogram : public NAList<FSInterval<E>>
{
public:
  /**
   * Constructs a histogram with the given number of intervals, all of the
   * same height.
   *
   * @param heap Heap to use in NAList ctor.
   * @param buckets Number of intervals to create for the new histogram.
   * @param height Height of each interval. Passed to each interval ctor as
   *               the initial number of array elements, but the interval's
   *               initial rowcount and frequency are left at 0.
   */
  FSHistogram(NAHeap* heap, Int32 buckets, Int32 height)
     : NAList<FSInterval<E>>(heap, buckets),
       heap_(heap), buckets_(buckets), height_(height)

  {
    for (CollIndex i=0; i<buckets; i++)
      this->insert(FSInterval<E>(heap, height));
  };

  ~FSHistogram()
  {};

  Int32 keyCountAt(Int32 bucketIdx);

  void increaseBuckets(NAHeap* heap, Int32 more)
  {
    Int32 n = this->entries();
    for (CollIndex i=0; i<more; i++)
      this->insertAt(n+i, FSInterval<E>(heap, 10));
  }

  void convertToEquiHeightHist(Int32 height, FSHistogram& equiHeightHistogram);

  void estimateRowsAndUecs(Float64 sampleRate, Float32 skRatio, HSColGroupStruct* group)
  {
    Float64 intervalCov;  // Coefficient of variation for a single interval
    Float64 totalCov = 0; // Sum of cov for intervals estimated
    NABoolean estimated = sampleRate < 1.0;

    for (CollIndex i=0; i<this->entries(); i++)
      {
        this->at(i).estimateRowsAndUecs(sampleRate, skRatio, intervalCov);
        if (estimated)
          totalCov += intervalCov;
      }

    if (estimated)
      group->coeffOfVar = totalCov / this->entries();
  }

  void display(ostream& out, NABoolean isEquiHeight)
  {
    if (isEquiHeight)
      out << "Equi-height Histogram:" << endl;
    else
      out << "Equi-width Histogram:" << endl;
    out << "Total intervals=" << this->entries() << endl;

    Int64 totalKeys = 0;
    Int64 totalFreq = 0;

    for (CollIndex i=0; i<this->entries(); i++)
      {
        FSInterval<E>& intv = this->at(i);

        totalKeys += (isEquiHeight ? intv.getSampleUec() : intv.entries());
        totalFreq += intv.getFreqCount();

        out << " intv[" << i << "]: keys="
            << (isEquiHeight ? intv.getSampleUec() : intv.entries())
            << ", totalFreq="
            << intv.getFreqCount()
            << '.';
        if (isEquiHeight)
          out << " ScaledUp: rc="
              << intv.getRc()
              << "; uec="
              << intv.getUec();
        out << endl;
      }

    cout << "total keys=" << totalKeys << endl;
    cout << "total frequency=" << totalFreq << endl;
  }

protected:
  Int32 buckets_;
  Int32 height_;
  NAHeap* heap_;
};

//Nonmember function defined in optimizer/EncodedValue.cpp.
extern Float64 EncVal_encodeString(const char * str, Lng32 strLen, CharType *cType);

// Templates for getting data ptr and length. For non-char types, the passed
// ptr and len are the true ones; for char types, the ptr is used to access
// the actual string and length.
template <class T>
inline char* getTrueDataPtr(T* ptr)
{ return (char*)ptr; }

template <class T>
inline UInt32 getTrueDataLen(T* ptr, Lng32 defaultLen)
{ return defaultLen; }

template <>
inline char* getTrueDataPtr(ISFixedChar* ptr)
{ return ptr->getContent(); }

template <>
inline UInt32 getTrueDataLen(ISFixedChar* ptr, Lng32)
{ return ptr->getLength(); }


template <>
inline char* getTrueDataPtr(ISVarChar* ptr)
//{ return ptr->getContent(); }
{ return ptr->getContent() + sizeof(UInt16); }

template <>
inline UInt32 getTrueDataLen(ISVarChar* ptr, Lng32)
//{ return *(UInt16*)(ptr->getContent()) + sizeof UInt16; }
{ return ptr->getLength(); }


// Templates for encoding data values. Non-char types have identity encoding
// (values encode as themselves). Char types are encoded as doubles.
template <class T, class E>
inline E encode(T& val)
{ return val; }

template <>
inline Float64 encode(ISFixedChar& fc)
{
  //printf("*** Encoding %.20s (length %d).\n", fc.getContent(), fc.getLength());
  return EncVal_encodeString(fc.getContent(), fc.getLength(), fc.getCharType());
}

template <>
inline Float64 encode(ISVarChar& vc)
{
  return EncVal_encodeString(vc.getContent() + sizeof(UInt16),
                             *(UInt16*)(vc.getContent()),
                             vc.getCharType());
}

// Templates for producing encoded value given a CBF key. For non-char types,
// just cast as the encoding type (which is the same as the actual type) and
// dereference.
template <class T, class E>
inline E encodeFromCbfKey(const simple_cbf_key& key)
{ return *((E*)key.getKey()); }

template <>
inline Float64 encodeFromCbfKey<ISFixedChar, Float64>(const simple_cbf_key& key)
{
  return EncVal_encodeString(key.getKey(), key.getKeyLen(), ISFixedChar::getCharType());
}

// varchar cbf key value is the string portion of ISVarChar content.
template <>
inline Float64 encodeFromCbfKey<ISVarChar, Float64>(const simple_cbf_key& key)
{
  return EncVal_encodeString(key.getKey(), key.getKeyLen(), ISVarChar::getCharType());
}

// Templates for preprocessing a rowset for a given column before adding the
// values to the CBF. Numeric columns require nothing, but for chars we have
// to set ptr and length fields in the objects representing the char values.
template <class T>
inline void setUpValues(HSColGroupStruct* group, Lng32 numRows)
{}

template <>
inline void setUpValues<ISFixedChar>(HSColGroupStruct* group, Lng32 numRows)
{
  // Set up elements of data array, which are pointers to char values.
  ISFixedChar* fixedCharPtr = (ISFixedChar*)group->data;
  char* strDataPtr = (char*)group->strData;
  for (Int32 i=0; i<numRows; i++)
    {
      fixedCharPtr->setContent(strDataPtr);
      strDataPtr += group->ISlength;
      fixedCharPtr++;
    }
}

template <>
inline void setUpValues<ISVarChar>(HSColGroupStruct* group, Lng32 numRows)
{
  // Set up elements of data array, which are pointers to char values.
  ISVarChar* varCharPtr = (ISVarChar*)group->data;
  char* strDataPtr = (char*)group->strData;
  for (Int32 i=0; i<numRows; i++)
    {
      varCharPtr->setContent(strDataPtr);
      strDataPtr += (group->ISlength + sizeof(UInt16) + (group->ISlength % 2));
      //strDataPtr += (*(UInt16*)strDataPtr + sizeof(UInt16));
      varCharPtr++;
    }
}


// Hash function used to get subsample frequency information from first rowset,
// to calculate estimated uec used as n parameter for CBF. The hash key for
// varchar columns is the entire content of the buffer representing an ISVarChar
// 2-byte len field and following char string.
template <class T>
inline ULng32 fsHashFunc(const T& key)
{
  return ExHDPHash::hash((char*)&key, ExHDPHash::NO_FLAGS, sizeof(T));
}

template <>
inline ULng32 fsHashFunc(const ISFixedChar& key)
{
  return ExHDPHash::hash(key.getContent(), ExHDPHash::NO_FLAGS, key.getLength());
}

template <>
inline ULng32 fsHashFunc(const ISVarChar& key)
{
  return ExHDPHash::hash(key.getContent(), ExHDPHash::NO_FLAGS, key.getLength() + sizeof(UInt16));
}

template <class T>
inline void recordKeyInfo(const simple_cbf_key&, Int64, Int64&)
{}

template <>
inline void recordKeyInfo<ISVarChar>(const simple_cbf_key& key, Int64 count, Int64& sum)
{
  //sum += *(UInt16*)(key.getKey());
  sum += (key.getKeyLen() * count);
}

template <class T, class E>
inline void decode(E encVal, T& decVal)
{
  decVal = encVal;
}

template <>
inline void decode(Float64 encVal, ISFixedChar& fchar)
{
  NAString* decodedValue = fchar.getCharType()->convertToString(encVal, STMTHEAP);
  char* copy = new(STMTHEAP) char[fchar.getLength() + 1];

  // Remove quotes added by decoder, add blank padding.
  memset(copy, ' ', fchar.getLength());
  strncpy(copy, decodedValue->data() + 1, decodedValue->length() - 2);
  copy[fchar.getLength()] = '\0';
  fchar.setContent(copy);
  delete decodedValue;
}

template <>
inline void decode(Float64 encVal, ISVarChar& vchar)
{
  NAString* decodedValue = vchar.getCharType()->convertToString(encVal, STMTHEAP);
  char* copy = new(STMTHEAP) char[decodedValue->length()
                                  + sizeof(UInt16)   // plus length field
                                  -2                 // minus trimmed quotes
                                  + 1];              // plus null terminator

  // Set length field, accounting for surrounding quotes that will be removed.
  *(UInt16*)copy = decodedValue->length() - 2;

  // Remove quotes added by decoder.
  strncpy(copy + sizeof(UInt16), decodedValue->data() + 1, decodedValue->length() - 2);
  copy[*(UInt16*)copy + sizeof(UInt16)] = '\0';
  vchar.setContent(copy);
  delete decodedValue;
}

template <class T>
inline void freeDecodedValue(T& val)
{}

template <>
inline void freeDecodedValue(ISFixedChar& fchar)
{
  NADELETEBASIC(fchar.getContent(), STMTHEAP);
  fchar.setContent(NULL);
}

template <>
inline void freeDecodedValue(ISVarChar& vchar)
{
  NADELETEBASIC(vchar.getContent(), STMTHEAP);
  vchar.setContent(NULL);
}

template <class T>
inline void fsDisplay(T& val)
{
  printf("%f", (Float64)val);
}

template <>
inline void fsDisplay<ISFixedChar>(ISFixedChar& val)
{
  printf("%.*s", val.getLength(), val.getContent());
}

template <>
inline void fsDisplay<ISVarChar>(ISVarChar& val)
{
  printf("%.*s", *(UInt16*)val.getContent(), val.getContent() + sizeof(UInt16));
}
