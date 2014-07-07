#ifndef MAPATTEMPTRESULT_HPP
#define MAPATTEMPTRESULT_HPP

/**
 * A triple to hold the return values from FMD::mapPosition() or
 * FMD::partialMap(). Holds a flag for whether the mapping succeeded or not,
 * a flag for whether the failure was due to multimapping and FMDPosition
 * corresponding either to where the character mapped or the longest
 * search starting at the character that did actually return results, and the
 * number of characters in the FMDPosition's search pattern.
 */
struct MapAttemptResult
{
    bool is_mapped;
    FMDPosition position;
    size_t characters;
    
    //Flag to indicate that we've multimapped
    bool is_multimap;

};

#endif
