// swigbindings.i: SWIG wrapper for the FMDIndex class.

// Name the module something Java-y that isn't going to clash with the FMD
// class.
%module FMDUtil

// Set up STL compatibility
%include "std_string.i"
%include "std_vector.i"

// Set up int types (int64_t)
%include "stdint.i"

// Set up exception not-killing-the-process-ness. See
// <http://www.swig.org/Doc1.3/Library.html#Library_nn17>
%include "exception.i"

// Java can't handle these operator names.
%rename(operatorLeftShift) operator<<;
%rename(operatorEquals) operator==;
%rename(operatorNotEquals) operator!=;
%rename(operatorIncrement) operator++;

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Note that build() on an FMDIndexBuilder produces a new object, which we
// should delete if we don't need it.
%newobject FMDIndexBuilder::build;

// Add some renames for the bit vector types.
%rename(RangeVector) NibbleVector;
%rename(RangeEncoder) NibbleEncoder;


%include "BitVector.hpp"
%include "NibbleVector.hpp"

typedef NibbleVector RangeVector;
typedef NibbleEncoder RangeEncoder;

%include "Mapping.hpp"
%include "TextPosition.hpp"
%{
  #include "FMDIndex.hpp"
%}
%include "FMDIndex.hpp"
%include "FMDIndexIterator.hpp"
%{
  #include "FMDIndexBuilder.hpp"
%}
%include "FMDIndexBuilder.hpp"





// We need to use the inner vector iterator classes to look at vectors. Give a
// partial definition under a new name.
class RangeVectorIterator
{
public:
  explicit RangeVectorIterator(const RangeVector& par);
  ~RangeVectorIterator();

  size_t rank(size_t value, bool at_least = false);

  size_t select(size_t index);
  size_t selectNext();

  pair_type valueBefore(size_t value);
  pair_type valueAfter(size_t value);
  pair_type nextValue();

  pair_type selectRun(size_t index, size_t max_length);
  pair_type selectNextRun(size_t max_length);

  bool isSet(size_t value);

  size_t countRuns();
};

// Since we will need to load and save range vectors to files, we need to expose
// a minimal C FILE API.
FILE* fopen(char* filename, char* mode);
void fclose(FILE* file);

// Java needs to work with vectors of mappings coming back from the map method.
%template(MappingVector) std::vector<Mapping>; 

// Java also needs to work with vectors of int64_ts coming back from the map
// method when working on ranges.
%template(IntVector) std::vector<int64_t>; 

// Whenever any of the JNI classes loads, load the native library.
%pragma(java) jniclasscode=%{
  static {
    FMDNativeLoader.load();
  }
%}

// We already worked around the inner classes thing.
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS
