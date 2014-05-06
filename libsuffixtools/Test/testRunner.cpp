// test.cpp: Main test runner for testing the libsuffixtools library. Uses
// CppUnit. Mostly borrowed from the CPPUnit documentation at <http://cppunit.so
// urceforge.net/doc/lastest/money_example.html#sec_running_test>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
#include <string>

/**
 * Main function: run all registered tests.
 */
int main(int argc, char** argv) {
    // Dump our arguments to avoid unused warnings, and so we can see them.
    std::cout << "Command line: ";
    for(int i = 0; i < argc; i++) {
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Running unit tests..." << std::endl;

    // Get the top level suite from the registry
    CppUnit::Test *suite = 
        CppUnit::TestFactoryRegistry::getRegistry().makeTest();

    // Adds the test to the list of test to run
    CppUnit::TextUi::TestRunner runner;
    runner.addTest( suite );

    // Run the tests.
    bool wasSucessful = runner.run();

    // Return error code 1 if the one of test failed.
    return wasSucessful ? 0 : 1;
}

#include <cppunit/extensions/HelperMacros.h>

/**
 * Test fixture.
 */
class BWTTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE( BWTTest );
    CPPUNIT_TEST( testConstruction );
    CPPUNIT_TEST_SUITE_END();
    
    // Keep a string saying where to get the haplotypes to test with.
    static const std::string filename;

public:
    void setUp();
    void tearDown();

    void testConstruction();
};

// Test fixture definition.

#include "Util/ReadTable.h"
#include "SuffixArray.cpp"

// Register the fixture to be run.
CPPUNIT_TEST_SUITE_REGISTRATION( BWTTest );

// Define constants
const std::string BWTTest::filename = "Test/haplotypes.fa";

void BWTTest::setUp() {
}


void BWTTest::tearDown() {
}

/**
 * Test building a BWT.
 */
void BWTTest::testConstruction() {
    
    // Make a read table from the headers in the file
	ReadTable* readTable = new ReadTable(filename);
    CPPUNIT_ASSERT(readTable != NULL);
    
    // Make a suffix array in 1 thread.
    SuffixArray* suffixArray = new SuffixArray(readTable, 1);
    CPPUNIT_ASSERT(suffixArray != NULL);
    
    // Validate it
    suffixArray->validate(readTable);
    
    // Try writing it out to nowhere.
    std::string devnull("/dev/null");
    suffixArray->writeBWT(devnull, readTable);
    suffixArray->writeIndex(devnull);
    
    // Clean up
    delete readTable;
    delete suffixArray;

}

