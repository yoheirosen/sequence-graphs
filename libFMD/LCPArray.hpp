#ifndef LCPARRAY_HPP
#define LCPARRAY_HPP

#include <vector>
#include <string>

// Depend on the libsuffixtools stuff.
#include <ReadTable.h>
#include <SuffixArray.h>

/**
 * Defines an array suitable for holding Longest Common Prefix information
 * between successive suffixes in a suffix array or FMD-index. Allows efficient
 * queries for the position of the next or previous smaller value from any
 * position.
 *
 * This particular implementation stores a full index, and so will use an
 * inordinately large amount of space.
 */
class LCPArray {

public:
    /**
     * Create a new LCPArray by building it from the given suffix array, using
     * the given collection of strings to get the actual characters.
     *
     * Neither the suffix array nor the read table are needed after the
     * constructor returns.
     */
    LCPArray(const SuffixArray& suffixArray, const ReadTable& strings);
    
    /**
     * Load an LCPArray from the given file. Uses platform-dependent byte
     * order and size_t size.
     */
    LCPArray(const std::string& filename);
    
    /**
     * Save an LCPArray to the given file. Uses platform-dependent byte order
     * and size_t size.
     */
    void save(const std::string& filename) const;
    
    /**
     * Get the longest common prefix value at a given index.
     */
    inline size_t operator[](size_t index) const {
        return values[index];
    }
    
    /**
     * Get the index of the previous smaller value before the given index.
     */
    inline size_t getPSV(size_t index) const {
        return psvs[index];
    }
    
    /**
     * Get the index of the next smaller value after the given index.
     */
    inline size_t getNSV(size_t index) const {
        return nsvs[index];
    }

protected:
    // TODO: Should these helpers go on the increasingly inaccurately named
    // ReadTable?

    /**
     * Get the character at the given offset into the given suffix in the given
     * collection of strings.
     */
    static inline char getFromSuffix(const SAElem& suffix, size_t offset,
        const ReadTable& strings) {
    
        // Split out the text and offset from the suffix, and offset the offset.
        return strings.getChar(suffix.getID(), suffix.getPos() + offset);    
        
    }
        
    /**
     * Get the length of a suffix form the collection of strings.
     */
    static inline size_t getSuffixLength(const SAElem& suffix,
        const ReadTable& strings) {
        
        // Get the length of the "read", minus the offset.
        return strings.getReadLength(suffix.getID()) - suffix.getPos();
        
    }
    
    // Store all the LCP array entries.
    std::vector<size_t> values;
    
    // Store the index of the previous smaller value for each position.
    std::vector<size_t> psvs;
    
    // Store the index of the next smaller value for each position.
    std::vector<size_t> nsvs;
    
private:
    
    /**
     * Helper method for the constructor. Visit the suffix tree node at the
     * given interval (start inclusive, end exclusive), find all its internal
     * nodes, fill in the LCP array and PSV/NSV indices for the boumdaries of
     * those internal nodes, and then recurse into each.
     */
    void exploreSuffixTreeNode(const SuffixArray& suffixArray, 
        const ReadTable& strings, size_t depth, size_t start, size_t end);
    
    /**
     * Binary-search the suffix array at the given depth for the given
     * character. Returns the first index of that character between start
     * (inclusive) and end (exclusive), or the position of the next larger
     * character (or end) if it is not present.
     */
    size_t binarySearch(const SuffixArray& suffixArray, 
        const ReadTable& strings, size_t depth, size_t start, size_t end,
        char value) const; 
    
};

#endif
