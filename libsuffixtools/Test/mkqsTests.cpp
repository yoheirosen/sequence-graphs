// Test the changes I made to quicksort.

#include <vector>

#include "Util/ReadTable.h"
#include "SuffixArray.h"
#include "SuffixCompare.h"
#include "Util/mkqs.h"

#include "mkqsTests.h"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( MKQSTests );

// Define constants
const std::string MKQSTests::filename = "Test/haplotypes.fa";

void MKQSTests::setUp() {

    // Set up the basic suffix array
    readTable = new ReadTable(filename);
    suffixArray = new SuffixArray(readTable, 1);

}


void MKQSTests::tearDown() {
    delete suffixArray;
    delete readTable;
}

/**
 * Test building a sampled suffix array.
 */
void MKQSTests::testSort() {
    
    // Make the comparison objects to use in the sort.
    SuffixCompareRadix radix_compare(readTable, 6);
    SuffixCompareIndex index_compare;
    
    // Get our own copy of the SA data since it's private
    SAElemVector elems;
    for(size_t i = 0; i < suffixArray->getSize(); i++) {
        elems.push_back(suffixArray->get(i));
    }
    
    // Sort the whole suffix array
    mkqs2(&elems[0], elems.size(), 0, radix_compare,
        index_compare);
    
    
    for(int i = 0; i < (int)elems.size() - 1; i++) {
        // Scan it and assert order.
        
        SAElem a = elems[i];
        SAElem b = elems[i + 1];
        
        // We need a < b
        
        // First compare the strings
        std::string stringA(radix_compare.getChrPtr(a));
        std::string stringB(radix_compare.getChrPtr(b));
        
        // a is smaller, or they're the same
        CPPUNIT_ASSERT(stringA <= stringB);
        
        if(stringA == stringB) {
        
            // If they're the same, index_compare must sort them right.
            CPPUNIT_ASSERT(index_compare(a, b));
        
        }
        
        
        
    }

}
