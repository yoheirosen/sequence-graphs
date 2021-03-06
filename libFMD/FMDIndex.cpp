#include <iostream>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <Util.h>

#include "FMDIndex.hpp"
#include "util.hpp"
#include "Log.hpp"

FMDIndex::FMDIndex(std::string basename, SuffixArray* fullSuffixArray): 
    names(), starts(), lengths(), cumulativeLengths(), genomeAssignments(),
    endIndices(), genomeRanges(), genomeMasks(), bwt(basename + ".bwt"), 
    suffixArray(basename + ".ssa"), fullSuffixArray(fullSuffixArray) {
    
    // TODO: Too many initializers

    Log::info() << "Loading " << basename << std::endl;

    // We already loaded the index itself in the initializer. Go load the
    // length/order metadata.
    
    // Open the contig name/start/length/genome file for reading.
    std::ifstream contigFile((basename + ".contigs").c_str());
    
    // Keep a cumulative length sum
    size_t lengthSum = 0;
    
    // Have a string to hold each line in turn.
    std::string line;
    while(std::getline(contigFile, line)) {
        // For each <contig>\t<start>\t<length>\t<genome> line...
        
        // Make a stringstream we can read out of.
        std::stringstream lineData(line);
        
        // Read in the name of the contig
        std::string contigName;
        lineData >> contigName;
        
        // Add it to the vector of names in number order.
        names.push_back(contigName);
        
        // Read in the contig start position on its scaffold
        size_t startNumber;
        lineData >> startNumber;
        
        // Add it to the vector of starts in number order.
        starts.push_back(startNumber);
        
        // Read in the contig's length
        size_t lengthNumber;
        lineData >> lengthNumber;
        
        // Add it to the vector of sizes in number order
        lengths.push_back(lengthNumber);
        
        // And to the vector of cumulative lengths
        cumulativeLengths.push_back(lengthSum);
        lengthSum += lengthNumber;
        
        // Read in the number of the genome that the contig belongs to.
        size_t genomeNumber;
        lineData >> genomeNumber;
        
        // Add it to the vector of genome assignments in number order
        genomeAssignments.push_back(genomeNumber);
    }
        
        
    // Close up the contig file. We read our contig metadata.
    contigFile.close();
    
    // Now read the genome bit masks.
    
    // What file are they in? Make sure to hold onto it while we construct the
    // stream with its c_str pointer.
    std::string genomeMaskFile = basename + ".msk";
    
    // Open the file where they live.
    std::ifstream genomeMaskStream(genomeMaskFile.c_str(), std::ios::binary);
    
    while(genomeMaskStream.peek() != EOF && !genomeMaskStream.eof()) {
        // As long as there is data left to read
        
        // Read a new BitVector from the stream and put it in our list.
        genomeMasks.push_back(new BitVector(genomeMaskStream));
        
        // This lets us autodetect how many genomes there are.
    }
    
    // Now invert the contig-to-genome index to make the genome-to-contig-range
    // index.
    
    // How many genomes are there?
    size_t numGenomes = genomeMasks.size();
    
    // Make the genome range vector big enough. Fill it with empty ranges for
    // genomes that somehow have no contigs.
    genomeRanges.resize(numGenomes, std::make_pair(0, 0));
    
    // Start out the first range.
    std::pair<size_t, size_t> currentRange = std::make_pair(0, 0);
    
    // Which genome are we on?
    size_t currentGenome;
    
    for(std::vector<size_t>::iterator i = genomeAssignments.begin(); 
        i != genomeAssignments.end(); ++i) {
    
        if(i == genomeAssignments.begin()) {
            // We're starting. Grab this genome as our current genome (since the
            // file may not start with genome 0).
            // TODO: Maybe just require that?
            currentGenome = *i;
        }
        
        if(*i == currentGenome) {
            // This is the same genome as the last one. Extend the range.
            currentRange.second++;
        } else {
            // This is now a different genome. Save the range.
            genomeRanges[currentGenome] = currentRange;
            
            // Make a new range that comes after it, with length 1, since we
            // include this contig.
            currentRange = std::make_pair(currentRange.second,
                currentRange.second + 1);
            
            // Remember what new genome we're looking at.
            currentGenome = *i;
            
            if(currentGenome >= numGenomes) {
                // Complain if we have a genome number higher than the number of
                // masks we loaded.
                throw std::runtime_error(
                    "Got a contig for a genome with no mask!");
            }
        }
    
    }
    
    // Save the last range
    genomeRanges[currentGenome] = currentRange;
    
    // First make sure the vector is big enough for them.
    endIndices.resize(getNumberOfContigs());
    
    for(int64_t i = 0; i < getNumberOfContigs() * 2; i++) {
        // The first #-of-texts rows in the BWT table have a '$' in the F
        // column, so the L column (what our BWT string actually is) will have
        // the last real character in some text.
        
        // Locate it to a text and offset
        TextPosition position = locate(i);
        
        if(position.getText() % 2 == 0) {
            // This is a forward strand. Save the index of the last real
            // character in the forward strand of the contig.
            endIndices[position.getText() / 2] = i;
        }
        
    }
    
    Log::info() << "Loaded " << names.size() << " contigs in " << numGenomes <<
        " genomes" << std::endl;
}

FMDIndex::~FMDIndex() {
    if(fullSuffixArray != NULL) {
        // If we were holding a full SuffixArray, throw it out.
        delete fullSuffixArray;
    }
    
    for(std::vector<BitVector*>::iterator i = genomeMasks.begin(); 
        i != genomeMasks.end(); ++i) {
        
        // Also delete all the genome masks we loaded.
        delete (*i);
    }
}

size_t FMDIndex::getContigNumber(TextPosition base) const {
    // What contig corresponds to that text? Contigs all have both strands.
    return base.getText() / 2;
}

bool FMDIndex::getStrand(TextPosition base) const {
    // What strand corresponds to that text? Strands are forward, then reverse.
    return base.getText() % 2 == 1;
}

size_t FMDIndex::getOffset(TextPosition base) const {
    // What base offset, 1-based, from the left, corresponds to this pair_type,
    // which may be on either strand.
    if(getStrand(base) == 0) {
        // We're on the forward strand. Make offset 1-based.
        return base.getOffset() + 1;
    } else {
        // We're on the reverse strand, so we measured from the end. Make it
        // 1-based.
        return getContigLength(getContigNumber(base)) - base.getOffset();
    }
}

std::string FMDIndex::getName(TextPosition base) const {

    // Unpack the coordinate parts.
    size_t contig = getContigNumber(base);
    size_t offset = getOffset(base);
    
    // Work out what to name the position.
    std::stringstream nameStream;
    // Leave off the strand.
    nameStream << "N" << contig << "B" << offset;
    return nameStream.str(); 
    
}

size_t FMDIndex::getBaseID(TextPosition base) const {
    // Get the cumulative total of bases by the start of the given contig
    size_t total = cumulativeLengths[getContigNumber(base)];
    
    // Add in the offset of this base from the start of its contig, convert back
    // to 0-based, and return.
    return total + getOffset(base) - 1;
}

size_t FMDIndex::getNumberOfContigs() const {
    // How many contigs do we know about?
    return names.size();
}
    
const std::string& FMDIndex::getContigName(size_t index) const {
    // Get the name of that contig.
    return names[index];
}

size_t FMDIndex::getContigStart(size_t index) const {
    // Get the start of that contig.
    return starts[index];
}

size_t FMDIndex::getContigLength(size_t index) const {
    // Get the length of that contig.
    return lengths[index];
}

size_t FMDIndex::getContigGenome(size_t index) const {
    // Get the genome that that contig belongs to.
    return genomeAssignments[index];
}
    
size_t FMDIndex::getNumberOfGenomes() const {
    // Get the number of genomes that are in the index.
    // TODO: Several things are this length. There should only be one.
    return genomeMasks.size();
}
    
std::pair<size_t, size_t> FMDIndex::getGenomeContigs(size_t genome) const {
    // Get the range of contigs belonging to the given genome.
    return genomeRanges[genome];
}

bool FMDIndex::isInGenome(int64_t bwtIndex, size_t genome) const {
    return BitVectorIterator(*genomeMasks[genome]).isSet(bwtIndex);
}

const BitVector& FMDIndex::getGenomeMask(size_t genome) const {
    return *genomeMasks[genome];
}

int64_t FMDIndex::getTotalLength() const {
    // Sum all the contig lengths and double (to make it be for both strands).
    // See <http://stackoverflow.com/a/3221813/402891>
    return std::accumulate(lengths.begin(), lengths.end(), 0) * 2;
}

int64_t FMDIndex::getBWTLength() const {
    return bwt.getBWLen();
}

FMDPosition FMDIndex::getCoveringPosition() const {
    // We want an FMDPosition that covers the entire BWT.
    
    // Just construct one. TODO: Will this extend properly since it goes to 0?
    return FMDPosition(0, 0, getBWTLength() - 1);
}
    
   
FMDPosition FMDIndex::getCharPosition(char c) const {
    // Starting a search with this character
    
    // See BWTAlgorithms::initInterval
    
    // Start the forward string with this character.
    int64_t forwardStart = bwt.getPC(c);
    
    // Start the reverse string with its complement.
    int64_t reverseStart = bwt.getPC(complement(c));
    
    // Get the offset to the end of the first interval (as well as the second).
    int64_t offset = bwt.getOcc(c, bwt.getBWLen() - 1) - 1;

    // Make the FMDPosition.
    return FMDPosition(forwardStart, reverseStart, offset);
    
}
     
void FMDIndex::extendFast(FMDPosition& range, char c, bool backward) const {
    // Extend the search with this character in an optimized way. We work on our
    // argument so we don't need to bother with copies.
    
    // Skip any sort of argument validation.
    
    if(!backward) {
        // Flip the arguments around so we can work on the reverse strand.
        c = complement(c);
        range.flipInPlace();
    }
    
    // Read occurrences of everything from the BWT
    
    // What rank among occurrences is the first instance of every character in
    // the BWT range?
    AlphaCount64 startRanks = bwt.getFullOcc(range.getForwardStart() - 1);
    
    // And the last? If endOffset() is 0, this will be 1 character later than
    // the call for startRanks, which is what we want.
    AlphaCount64 endRanks = bwt.getFullOcc(range.getForwardStart() + 
        range.getEndOffset());
        
    // Get the number of suffixes that had '$' (end of text) next. TODO: should
    // this be '\0' instead?
    range.setReverseStart(range.getReverseStart() + 
        (endRanks.get('$') - startRanks.get('$')));
        
    for(size_t base = 0; base < NUM_BASES; base++) {
        // For each base in alphabetical order by reverse complement
        
        // Work out the length of the interval this base gets.
        size_t intervalLength = endRanks.get(BASES[base]) - 
            startRanks.get(BASES[base]);
        
        if(BASES[base] == c) {
            // This is the base we're looking for. Finish up and break out of
            // the loop.
            
            // Range reverse start is already set.
            
            // Set the range forward start.
            range.setForwardStart(bwt.getPC(c) + startRanks.get(c));
            
            // Set the range length.
            range.setEndOffset((int64_t)intervalLength - 1);
            
            // Now we've put together the range.
            break;
            
        } else {
            // This is not the base we're looking for. Budge the reverse strand
            // interval over by the length of the interval, to account for the
            // bit this base took up.
            range.setReverseStart(range.getReverseStart() + intervalLength);
        }
    }
    
    if(!backward) {
        // Flip the result since we were working on the opposite strand.
        range.flipInPlace();
    }
    
}
   
FMDPosition FMDIndex::extend(FMDPosition range, char c, bool backward) const {
    // Extend the search with this character.
    
    // More or less directly implemented off of algorithms 2 and 3 in "Exploring
    // single-sample SNP and INDEL calling with whole-genome de novo assembly"
    // (Li, 2012). However, our character indices are one less, since we don't
    // allow search patterns to include the end-of-text symbol. We also use
    // alphabetical ordering instead of the paper's N-last ordering in the FM-
    // index, and consequently need to assign reverse ranges in alphabetical
    // order by reverse complement.
  
    if(!backward) {
        // We only really want to implement backwards search. Flip the interval,
        // do backwards search with the complement of the base, and then flip
        // back.
        return extend(range.flip(), complement(c), true).flip();        
    }

    if(c == '\0') {
        throw std::runtime_error("Can't extend with null byte!");
    }

    if(!isBase(c)) {
        std::string errorMessage = std::string("Character #");
        errorMessage.push_back(c);
        errorMessage += std::string(" is not a DNA base.");
        throw std::runtime_error(errorMessage);
    }

    Log::trace() << "Extending " << range << " backwards with " << c <<
        std::endl;

    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];

    for(size_t base = 0; base < NUM_BASES; base++) {
        // Go through the bases in arbitrary order.

        Log::trace() << "\tThinking about base " << base << "(" << 
            BASES[base] << ")" << std::endl;

        // Count up the number of characters < this base.
        int64_t start = bwt.getPC(c);

        Log::trace() << "\t\tstart = " << start << std::endl;

        // Get the rank among occurrences of the first instance of this base in
        // this slice.
        int64_t forwardStartRank = bwt.getOcc(BASES[base], 
            range.getForwardStart() - 1);
        
        // Get the same rank for the last instance. TODO: Is the -1 right here?
        int64_t forwardEndRank = bwt.getOcc(BASES[base], 
            range.getForwardStart() + range.getEndOffset()) - 1;

        // Fill in the forward-strand start position and range end offset for
        // this base's answer.
        answers[base].setForwardStart(start + forwardStartRank);
        answers[base].setEndOffset(forwardEndRank - forwardStartRank);

        Log::trace() << "\t\tWould go to: " << answers[base] << std::endl;
    }

    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character (i.e. how much of the current range
    // is devoted to things where an endOfText comes next) implicitly: it's
    // whatever part of the length of the range is unaccounted-for by the other
    // characters. We need to use the length accessor because ranges with one
    // thing have the .end_offset set to 0.
    int64_t endOfTextLength = range.getLength();

    for(size_t base = 0; base < NUM_BASES; base++)
    {
        // Go through the bases in order and account for their lengths.
        endOfTextLength -= answers[base].getLength();
    }


    Log::trace() << "\tendOfTextLength = " << endOfTextLength << std::endl;

    // The endOfText character is the very first character we need to account
    // for when subdividing the reverse range and picking which subdivision to
    // take.
    Log::trace() << "\tendOfText reverse_start would be " << 
        range.getReverseStart() << std::endl;

    // Next, allocate the range for the base that comes first in alphabetical
    // order by reverse complement.
    answers[0].setReverseStart(range.getReverseStart() + endOfTextLength);
    Log::trace() << "\t" << BASES[0] << " reverse_start is " << 
        answers[0].getReverseStart() << std::endl;

    for(size_t base = 1; base < NUM_BASES; base++)
    {
        // For each subsequent base in alphabetical order by reverse complement
        // (as stored in BASES), allocate it the next part of the reverse range.

        answers[base].setReverseStart(answers[base - 1].getReverseStart() + 
            answers[base - 1].getLength());
        Log::trace() << "\t" << BASES[base] << " reverse_start is " << 
        answers[base].getReverseStart() << std::endl;
    }

    // Now all the per-base answers are filled in.

    for(size_t base = 0; base < NUM_BASES; base++)
    {
        // For each base in arbitrary order
        if(BASES[base] == c)
        {
            // This is the base we're actually supposed to be extending with. Return
            // its answer.
            
            Log::trace() << "Moving " << range << " to " << answers[base] << 
                " on " << BASES[base] << std::endl;
            
            return answers[base];
        }
    }

    // If we get here, they gave us something not in BASES somehow, even though
    // we checked already.
    throw std::runtime_error("Unrecognized base");
}

FMDPosition FMDIndex::count(std::string pattern) const {
    if(pattern.size() == 0) {
        // We match everything! Say the whole range of the BWT.
        return getCoveringPosition();
    }
    

    // Start at the end and select the first character.
    FMDPosition position = getCharPosition(pattern[pattern.size() - 1]);
    
    for(int i = pattern.size() - 2; !position.isEmpty() && i >= 0; i--) {
        // Extend backwards with each character
        position = extend(position, pattern[i], true);
    }
    
    // We either ran out of matching locations or finished the pattern.
    return position;

}

TextPosition FMDIndex::locate(int64_t index) const {
    // Wrap up locate functionality.
    
    // We're going to fill in an SAElem: a composite thing where the high bits
    // encode the text number, and the low bits encode the offset.
    SAElem bitfield;

    if(fullSuffixArray != NULL) {
        // We can just look at the full suffix array cheat sheet.
        bitfield = fullSuffixArray->get(index);
        
    } else {
        // We need to use the sampled suffix array.
        
        // Run the libsuffixtools locate. 
        bitfield = suffixArray.calcSA(index, &bwt); 
        
    }
    
    // Unpack it and convert to our own format.
    return TextPosition(bitfield.getID(), bitfield.getPos());
}

int64_t FMDIndex::getContigEndIndex(size_t contig) const {
    // Looks a bit like the metadata functions from earlier. Actually pulls info
    // from the same file.
    return endIndices[contig];
}

char FMDIndex::display(int64_t index) const {
    // Just pull straight from the BWT string.
    return bwt.getChar(index);
}

char FMDIndex::displayFirst(int64_t index) const {
    // Our BWT supports this natively.
    return bwt.getF(index);
}

std::string FMDIndex::displayContig(size_t index) const {
    // We can't efficiently un-locate, so we just use a vector of the last BWT
    // index in every contig. This works since there are no 0-length contigs.
    int64_t bwtIndex = getContigEndIndex(index);
    
    // Make a string to hold all the bases.
    std::string bases;

    for(size_t i = 0; i < getContigLength(index); i++) {
        // Until we have the right number of bases...
        
        // Grab this base
        bases.push_back(display(bwtIndex));
        
        // LF-map to the previous position in the contig (or off the front end).
        bwtIndex = getLF(bwtIndex);
    }
    
    // Flip the string around so it's front to front.
    std::reverse(bases.begin(), bases.end());
    
    // Give it back.
    return bases;
    
}

int64_t FMDIndex::getLF(int64_t index) const {
    // Find the character we're looking at
    char toFind = display(index);
    
    // Find the start of that character in the first column. It's just the
    // number of characters less than it, counting text stops.
    int64_t charBlockStart = bwt.getPC(toFind);
    
    // Find the rank of that instance of that character among instances of the
    // same character in the last column. Subtract 1 from occurrences since the
    // first copy should be rank 0.
    int64_t instanceRank = bwt.getOcc(toFind, index) - 1;
    
    // Add that to the start position to produce the LF mapping.
    return charBlockStart + instanceRank;
}

std::vector<Mapping> FMDIndex::map(const std::string& query,
    const BitVector* mask, int minContext, int start, int length) const {
	
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }

    // Make an itarator for the mask, if needed, so we can query it.
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
        new BitVectorIterator(*mask);
        
    if(maskIterator == NULL) {
        Log::debug() << "Mapping " << length << " bases to all genomes." <<
            std::endl;
    } else {
        Log::debug() << "Mapping " << length << " bases to one genome only." <<
            std::endl;
    }
    
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
        std::endl;

    // We need a vector to return.
    std::vector<Mapping> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to track our FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base.
    // Other fields get overwritten.
    location.position = EMPTY_FMD_POSITION;

    for(size_t i = start; i < start + length; i++)
    {
        if(location.position.isEmpty(maskIterator))
        {
            Log::debug() << "Starting over by mapping position " << i <<
                std::endl;
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->mapPosition(query, i, maskIterator);
        } else {
            Log::debug() << "Extending with position " << i << std::endl;
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the right (not
            // backwards) with the next base.
            location.position = this->extend(location.position, query[i],
                false);
            location.characters++;
        }

        if(location.is_mapped && location.characters >= minContext &&
            location.position.getLength(maskIterator) == 1) {
            
            // It mapped. We didn't do a re-start and fail, we have enough
            // context to be confident, and there's exactly one thing in our
            // interval.

            // Take the first (only) thing in the bi-interval's forward strand
            // side, not accounting for the mask.
            int64_t start = location.position.getForwardStart();
            
            if(maskIterator != NULL) {
                // Account for the mask. The start position of the interval may
                // be masked out. Get the first 1 after (or at) the start,
                // instead of the start itself. Since the interval is nonempty
                // under the mask, we know this will exist.
                start = maskIterator->valueAfter(start).first;
            }

            // Locate it, and then report position as a (text, offset) pair.
            // This will give us the position of the first base in the pattern,
            // which lets us infer the position of the last base in the pattern.
            TextPosition textPosition = locate(start);

            Log::debug() << "Mapped " << location.characters << "/" << 
                minContext << " context to text " << textPosition.getText() << 
                " position " << textPosition.getOffset() << std::endl;

            // Correct to the position of the last base in the pattern, by
            // offsetting by the length of the pattern that was used. A
            // 2-character pattern means we need to go 1 further right in the
            // string it maps to to find where its rightmost character maps.
            textPosition.setOffset(textPosition.getOffset() + 
                (location.characters - 1));

            // Add a Mapping for this mapped base.
            mappings.push_back(Mapping(textPosition));

            // We definitely have a non-empty FMDPosition to continue from

        } else {

            Log::debug() << "Failed (" << 
                location.position.getLength(maskIterator) << " options for " <<
                location.characters << " context)." << std::endl;

            if(location.is_mapped && location.position.isEmpty(maskIterator)) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                Log::debug() << "Restarting from here..." << std::endl;

                // Move the loop index back
                i--;

                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base.

            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little left context to 
                //   be unique
                // - It was an initial mapping with a nonexistent left context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Add an empty/unmapped Mapping.
                mappings.push_back(Mapping());

                // Mark that the next iteration will be an extension (if we had
                // any results this iteration; if not it will just restart)
                location.is_mapped = true;

            }
        }
    }

    // We've gone through and attempted the whole string. Give back our answers.
    return mappings;

}

std::vector<Mapping> FMDIndex::map(const std::string& query, int64_t genome, 
    int minContext, int start, int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return map(query, genome == -1 ? NULL : genomeMasks[genome], minContext, 
        start, length);    
}

char op_increase (char i) { return ++i; }

std::vector<Mapping> FMDIndex::mapBoth(const std::string& query, int64_t genome, 
    int minContext, int start, int length) const {
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        
        // We need to do this ourselves since we go clipping out that bit of the
        // string to reverse complement.
        length = query.length() - start;
    }
    
    // Map it forward
    std::vector<Mapping> forward = map(query, genome, minContext, start,
        length);
    
    // Make a reversed copy of the appropriate region of the query string.
    
    // Where does our selected region end (as a reverse iterator)?
    std::string::const_reverse_iterator reverseStart = query.rbegin() + 
        (query.size() - length);
    
    // And what's one before it started?
    std::string::const_reverse_iterator reverseEnd = query.rend() - start;
    
    // Copy the string over, complementing bases as we go with libsuffixtools's
    // complement function that goes char -> char. See
    // <http://stackoverflow.com/a/7531885/402891>. Also remember to make a
    // back_inserter so we actually can put in new characters.
    std::string reverseComplemented;
    std::transform(reverseStart, reverseEnd, 
        std::back_inserter(reverseComplemented), (char(*)(char))complement);
        
    // Map it backward
    std::vector<Mapping> reverse = map(reverseComplemented, genome, minContext);
    
    if(forward.size() != reverse.size()) {
        throw std::runtime_error("Forward and reverse region size mismatch!");
    }
    
    for(size_t i = 0; i < forward.size(); i++) {
        // Go through and disambiguate in place to resolve multi-mappings and
        // such. Make sure to read reverse backwards.
        
        forward[i] = disambiguate(forward[i], reverse[reverse.size() - i - 1]);
    }
    
    // Give back the disambiguated vector.
    return forward;
    
}

std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> FMDIndex::Cmap(const BitVector& ranges,
    const std::string& query, const BitVector* mask, int minContext, int start,
    int length) const {
    
    // Map to a range.
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }
    
    Log::debug() << "Mapping with (two-sided) minimum " << minContext << " context." <<
        std::endl;

    // Make an iterator for ranges, so we can query it.
    BitVectorIterator rangeIterator(ranges);
    
    // And one for the mask, if needed
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
        new BitVectorIterator(*mask);

    // We need a vector to return.
    std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to trackour FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    creditMapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base
    location.position = EMPTY_FMD_POSITION;

    for(int i = start + length - 1; i >= start; i--) {
        // Go from the end of our selected region to the beginning.
	
        Log::trace() << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;
	    
	// Need to prevent an extension from overrunning the query contig, if
	// this will happen, also trigger a restart
	    	    
        if(location.position.isEmpty() || i < location.characters) {
            Log::debug() << "Starting over by mapping position " << i << std::endl;
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->CmapPosition(rangeIterator, query, i, maskIterator);
        } else {
            Log::debug() << "Extending with position " << i << " with characters = " << location.characters << std::endl;
	    // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the left
            // (backwards) with the next base.
            location.position = this->extend(location.position, query[i - location.characters + 1], true);
	    location.position = this->extend(location.position, query[i - location.characters], true);
            location.characters++;
	    if(location.characters > location.maxCharacters) {
		location.maxCharacters++;
	    }
        }

        // What range index does our current left-side position (the one we just
        // moved) correspond to, if any?
        int64_t range = location.position.range(rangeIterator, maskIterator);
	
	if(location.characters < minContext && location.maxCharacters >=minContext) {
	    location.characters = minContext;
	}

        if(location.is_mapped && location.characters >= minContext && 
            !location.position.isEmpty(maskIterator) && range != -1) {
            
            // It mapped. We didn't do a re-start and fail, we have sufficient
            // context to be confident, and our interval is nonempty and
            // subsumed by a range.

            Log::debug() << i << " Mapped " << location.characters << 
                " context to " << location.position << " in range #" << range <<
                std::endl;

            // Remember that this base mapped to this range
            mappings.push_back(std::make_pair(range,std::make_pair(location.characters,location.maxCharacters)));
            
            // We definitely have a non-empty FMDPosition to continue from

        } else {

            Log::debug() << "Failed at " << i << " " << location.position << " (" << 
                location.position.ranges(rangeIterator, maskIterator) <<
                " options for " << location.characters << " context)." << 
                std::endl;
                
            if(location.is_mapped && location.position.isEmpty(maskIterator)) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                Log::debug() << "Restarting from here..." << std::endl;

                // Move the loop index towards the end we started from (right)
                i++;
		
                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base.

            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little right context to 
                //   be unique to a range.
                // - It was an initial mapping with a nonexistent right context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Say it corresponds to no range.
                mappings.push_back(std::make_pair(-1,std::make_pair(0,0)));

                // Mark that the next iteration will be an extension (if we had
                // any results this iteration; if not it will just restart)
                location.is_mapped = true;

            }
        }
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Get rid of the mask iterator if needed
    if(maskIterator != NULL) {
        delete maskIterator;
    }

    // Give back our answers.
    return mappings;
    
    
}

std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> FMDIndex::Cmap(const BitVector& ranges, 
    const std::string& query, int64_t genome, int minContext, int start,
    int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return Cmap(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext, start, length);    
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::map(const BitVector& ranges,
    const std::string& query, const BitVector* mask, int minContext, int start,
    int length) const {
    
    // RIGHT-map to a range.
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }
    
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
        std::endl;

    // Make an iterator for ranges, so we can query it.
    BitVectorIterator rangeIterator(ranges);
    
    // And one for the mask, if needed
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
        new BitVectorIterator(*mask);

    // We need a vector to return.
    std::vector<std::pair<int64_t,size_t>> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to trackour FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MapAttemptResult location;
    // Make sure the scratch position is empty so we re-start on the first base
    location.position = EMPTY_FMD_POSITION;

    for(int i = start + length - 1; i >= start; i--) {
        // Go from the end of our selected region to the beginning.

        Log::trace() << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;

        if(location.position.isEmpty()) {
            Log::debug() << "Starting over by mapping position " << i <<
                std::endl;
            // We do not currently have a non-empty FMDPosition to extend. Start
            // over by mapping this character by itself.
            location = this->mapPosition(rangeIterator, query, i, maskIterator);
        } else {
            Log::debug() << "Extending with position " << i << std::endl;
            // The last base either mapped successfully or failed due to multi-
            // mapping. Try to extend the FMDPosition we have to the left
            // (backwards) with the next base.
            location.position = this->extend(location.position, query[i], true);
            location.characters++;
        }

        // What range index does our current left-side position (the one we just
        // moved) correspond to, if any?
        int64_t range = location.position.range(rangeIterator, maskIterator);

        if(location.is_mapped && location.characters >= minContext && 
            !location.position.isEmpty(maskIterator) && range != -1) {
            
            // It mapped. We didn't do a re-start and fail, we have sufficient
            // context to be confident, and our interval is nonempty and
            // subsumed by a range.

            Log::debug() << "Mapped " << location.characters << 
                " context to " << location.position << " in range #" << range <<
                std::endl;

            // Remember that this base mapped to this range
            mappings.push_back(std::make_pair(range,location.characters - 1));
	    
            
            // We definitely have a non-empty FMDPosition to continue from

        } else {

            Log::debug() << "Failed at " << location.position << " (" << 
                location.position.ranges(rangeIterator, maskIterator) <<
                " options for " << location.characters << " context)." << 
                std::endl;
                
            if(location.is_mapped && location.position.isEmpty(maskIterator)) {
                // We extended right until we got no results. We need to try
                // this base again, in case we tried with a too-long left
                // context.

                Log::debug() << "Restarting from here..." << std::endl;

                // Move the loop index towards the end we started from (right)
                i++;

                // Since the FMDPosition is empty, on the next iteration we will
                // retry this base.

            } else {
                // It didn't map for some other reason:
                // - It was an initial mapping with too little right context to 
                //   be unique to a range.
                // - It was an initial mapping with a nonexistent right context
                // - It was an extension that was multimapped and still is

                // In none of these cases will re-starting from this base help
                // at all. If we just restarted here, we don't want to do it
                // again. If it was multimapped before, it had as much left
                // context as it could take without running out of string or
                // getting no results.

                // It didn't map. Say it corresponds to no range.
                mappings.push_back(std::make_pair(-1,0));

                // Mark that the next iteration will be an extension (if we had
                // any results this iteration; if not it will just restart)
                location.is_mapped = true;

            }
        }
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Get rid of the mask iterator if needed
    if(maskIterator != NULL) {
        delete maskIterator;
    }

    // Give back our answers.
    return mappings;
    
    
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::map(const BitVector& ranges, 
    const std::string& query, int64_t genome, int minContext, int start,
    int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return map(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext, start, length);    
}

FMDIndex::iterator FMDIndex::begin(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that automatically searches out the first
    // suffix of the right length.
    return FMDIndex::iterator(*this, depth, false, reportDeadEnds);
}
     
FMDIndex::iterator FMDIndex::end(size_t depth, bool reportDeadEnds) const {
    // Make a new suffix tree iterator that is just a 1-past-the-end sentinel.
    return FMDIndex::iterator(*this, depth, true, reportDeadEnds);
}

MapAttemptResult FMDIndex::mapPosition(const std::string& pattern,
    size_t index, BitVectorIterator* mask) const {

    Log::debug() << "Mapping " << index << " in " << pattern << std::endl;
  
    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MapAttemptResult result;

    // Do a backward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);
    result.characters = 1;
    if(result.position.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if(result.position.getLength(mask) == 1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    if(index == 0) {
        // The rest of the function deals with characters to the left of the one
        // we start at. If we start at position 0 there can be none.
        return result;
    }

    Log::trace() << "Starting with " << result.position << std::endl;

    do {
        // Now consider the next character to the left.
        index--;

        // Grab the character to extend with.
        char character = pattern[index];

        Log::trace() << "Index " << index << " in " << pattern << " is " << 
            character << "(" << character << ")" << std::endl;

        // Backwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position, character,
            true);

        Log::trace() << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;
        if(next_position.isEmpty(mask)) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        } else if(next_position.getLength(mask) == 1) {
            // We have successfully mapped to exactly one place. Update our
            // result to reflect the additional extension and our success, and
            // return it.
            result.position = next_position;
            result.characters++;
            result.is_mapped = true;
            return result;      
        }

        // Otherwise, we still map to a plurality of places. Record the
        // extension and loop again.
        result.position = next_position;
        result.characters++;

    } while(index > 0);
    // Continue considering characters to the left until we hit the start of the
    // string.

    // If we get here, we ran out of upstream context and still map to multiple
    // places. Just give our multi-mapping FMDPosition and unmapped result.
    return result;
}

creditMapAttemptResult FMDIndex::CmapPosition(BitVectorIterator& ranges, 
    const std::string& pattern, size_t index, BitVectorIterator* mask) const {
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    creditMapAttemptResult result;

    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);    
    result.characters = 1;
    result.maxCharacters = 1;
    if(result.position.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if (result.position.range(ranges, mask) != -1) {
        // We've already mapped.
        result.is_mapped = true;
    }

    Log::trace() << "Starting with " << result.position << std::endl;

    FMDPosition found_position;
    FMDPosition next_position;
    
    for(size_t i = 1; index + i < pattern.size() && 1 + index > i; i++) {
		
        // Dual extend with subsequent characters.
        next_position = this->extend(result.position,
            pattern[index + i], false);
	next_position = this->extend(next_position,pattern[index - i], true);
	
        Log::debug() << "Now at " << next_position << " after " << pattern[i] << std::endl;
        if(next_position.isEmpty(mask)) {
            // The next place we would go is empty, so return the result holding
            // the last position.
	    Log::debug() << "Couldn't find more context" << std::endl;
	    result.characters = result.maxCharacters;
            return result;
        }

        if(!result.is_mapped && next_position.range(ranges, mask) != -1) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success, and
            // return it.
    
    	    result.position = next_position;
	    result.maxCharacters++;
	    result.characters = result.maxCharacters;
	    Log::debug() << "Extended " << i << " times" << std::endl;
	    result.is_mapped = true;
	    found_position = result.position;
	    
        } else if(result.is_mapped && next_position.range(ranges, mask) != -1) {
	    result.position = next_position;
	    result.maxCharacters++;
	    Log::debug() << "Restart continue " << i << std::endl;

	} else {
	    // Otherwise, we still map to a plurality of ranges. Record the
	    // extension and loop again.
	
	    result.position = next_position;
	    result.maxCharacters++;
	    result.characters = result.maxCharacters;
	    
	}
    }
    
    if(result.is_mapped) {
	result.position = found_position;
    }

    // If we get here, we ran out of downstream context and still map to
    // multiple ranges. Just give our multi-mapping FMDPosition and unmapped
    // result.
    return result;


}

MapAttemptResult FMDIndex::mapPosition(BitVectorIterator& ranges, 
    const std::string& pattern, size_t index, BitVectorIterator* mask) const {
    
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MapAttemptResult result;

    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.position = this->getCharPosition(pattern[index]);
    result.characters = 1;
    if(result.position.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if (result.position.range(ranges, mask) != -1) {
        // We've already mapped.
        result.is_mapped = true;
        return result;
    }

    Log::trace() << "Starting with " << result.position << std::endl;
    
    FMDPosition found_position;

    for(index++; index < pattern.size(); index++) {
        // Forwards extend with subsequent characters.
        FMDPosition next_position = this->extend(result.position,
            pattern[index], false);

        Log::trace() << "Now at " << next_position << " after " << 
            pattern[index] << std::endl;
        if(next_position.isEmpty(mask)) {
            // The next place we would go is empty, so return the result holding
            // the last position.
            return result;
        }

        if(next_position.range(ranges, mask) != -1) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success, 
            // but continue the search
	  
	    result.position = next_position;
	    result.characters++;
	    result.is_mapped = true;
	    found_position = result.position;
	    
        } else {
	    // Otherwise, we still map to a plurality of ranges. Record the
	    // extension and loop again.
	
	    result.position = next_position;
	    result.characters++;
	
	}
    }
    
    if(result.is_mapped) {
	result.position = found_position;
    }
    
    return result;

}

Mapping FMDIndex::disambiguate(const Mapping& left, 
    const Mapping& right) const {

    if(right.is_mapped) {
        
        // Either left is mapped, and we need to compare right to it, or it
        // isn't, and we need to return a left-ified version of right. Either
        // way we need to work out what right would look like on the left.
        
        // That requires the length of the contig it is on.
        size_t contigLength = getContigLength(getContigNumber(right.location));
            
        // What's the new text we map to? Just toggle the low bit (since text 2n
        // and 2n + 1 go together).
        size_t flippedText = right.location.getText() ^ 1;
        
        // And what's our offset on that text? -1 because we want a 0-based
        // answer still.
        size_t flippedOffset = contigLength - right.location.getOffset() - 1;
    
        // TODO: as an optimization, if left did map, try only computing the
        // above as they are needed for the comparison.
    
        // TODO: Add a flip method to TextPositions (taking contig length) and
        // use that.
    
        if(left.is_mapped) {
            // Check if they are mapped to the same base
            
            if(left.location.getText() == flippedText && 
                left.location.getOffset() == flippedOffset) {
                
                // We successfully 2-sided-mapped. Return either one.
                return left;
            } else {
                // We didn't map at all. Construct an unmapped answer.
                return Mapping(TextPosition(0, 0), false);
            }
            
            
        } else {
            // Only right is mapped. Return a flipped version of it.
            return Mapping(TextPosition(flippedText, flippedOffset));
        }
    } else {
        // Right isn't mapped, so whether left is mapped or not, we can return
        // it and be correct.
        return left;
    }
}

MisMatchAttemptResults FMDIndex::misMatchExtend(MisMatchAttemptResults& prevMisMatches,
	char c, bool backward, size_t z_max, BitVectorIterator* mask, bool startExtension, bool finishExtension) const {
    MisMatchAttemptResults nextMisMatches;
    nextMisMatches.is_mapped = prevMisMatches.is_mapped;
    nextMisMatches.characters = prevMisMatches.characters;
    
    // Note that we do not flip parameters when !backward since
    // FMDIndex::misMatchExtend uses FMDIndex::extend which performs
    // this step itself
    
    if(prevMisMatches.positions.size() == 0) {
	throw std::runtime_error("Tried to extend zero length mismatch vector");
    }
    
    if(prevMisMatches.positions.front().first.isEmpty(mask)) {
	throw std::runtime_error("Can't extend an empty position");
    }

    if(c == '\0') {
        throw std::runtime_error("Can't extend with null byte!");
    }

    if(!isBase(c)) {
        std::string errorMessage = std::string("Character #");
        errorMessage.push_back(c);
        errorMessage += std::string(" is not a DNA base.");
        throw std::runtime_error(errorMessage);
    }
    
    std::pair<FMDPosition,size_t> m_position;
    std::pair<FMDPosition,size_t> m_position2;
        
    for(std::vector<std::pair<FMDPosition,size_t>>::iterator it =
	  prevMisMatches.positions.begin(); it != prevMisMatches.positions.end(); ++it) {
	
	// Store each successive element as an m_position

	m_position.first = it->first;
	m_position.second = it->second;
		
	// extend m_position by correct base. Do not do this if the
	// finishExtension flag is true--in this case it's already been
	// done
    
	if(startExtension) {
		
	    m_position2.first = extend(m_position.first, c, backward);
	    m_position2.second = m_position.second;
	
	    if(m_position2.first.getLength(mask) > 0) {
		nextMisMatches.positions.push_back(m_position2);	
	    }
	    
	} else if(finishExtension) {	    
	    // Extend by all mismatched positions
	    if(m_position.second < z_max) {
		for(size_t base = 0; base < NUM_BASES; base++) {
		    if(BASES[base] != c) {
			m_position2.first = extend(m_position.first, BASES[base], backward);
			m_position2.second = m_position.second;
			m_position2.second++;
		    
			// If the position exists at all in the FMDIndex, place
			// it in the results vector
		    
			if(m_position2.first.getLength(mask) > 0) {
			    nextMisMatches.positions.push_back(m_position2);
			}
		    }
		}
	    }
	} else {
	  
	    m_position2.first = extend(m_position.first, c, backward);
	    m_position2.second = m_position.second;
	
	    if(m_position2.first.getLength(mask) > 0) {
		nextMisMatches.positions.push_back(m_position2);
		
	    }
	    
	    if(m_position.second < z_max) {
		for(size_t base = 0; base < NUM_BASES; base++) {
		    if(BASES[base] != c) {
			m_position2.first = extend(m_position.first, BASES[base], backward);
			m_position2.second = m_position.second;
			m_position2.second++;
		    
			// If the position exists at all in the FMDIndex, place
			// it in the results vector
		    
			if(m_position2.first.getLength(mask) > 0) {
				nextMisMatches.positions.push_back(m_position2);
		    
			} 
		    }
		}
	    }
	}		    
    }
    
    // If no results are found, place an empty FMDPosition in the
    // output vector
        
    if(nextMisMatches.positions.size() == 0) {
	nextMisMatches.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
    }
    
    // Return all matches
    
    // Or if there are matches, but not unique matches of at least
    // minimum context length, return the entire vector of positions
    // generated in this run, to use as starting material for the next
            
    return nextMisMatches;
    
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::misMatchMap(const BitVector& ranges,
    const std::string& query, const BitVector* mask, int minContext, size_t z_max,
    int start, int length) const {
    
    if(length == -1) {
	// Fix up the length parameter if it is -1: that means the whole rest of
	// the string.
	length = query.length() - start;
    }
	
    Log::debug() << "Mapping with minimum " << minContext << " context." <<
	std::endl;
	
    // Make an iterator for ranges, so we can query it.
    BitVectorIterator rangeIterator(ranges);

    // And one for the mask, if needed
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
    new BitVectorIterator(*mask);
    
    // We need a vector to return.
    std::vector<std::pair<int64_t,size_t>> mappings;
    
    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to trackour FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.
    MisMatchAttemptResults search;
    // Make sure the scratch position is empty so we re-start on the first base
    search.positions.push_back(
	std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
    search.characters = 0;
    search.maxCharacters = 0;

    MisMatchAttemptResults searchExtend;
    // Make sure the scratch position is empty so we re-start on the first base
    searchExtend.positions.push_back(
	std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
    search.characters = 0;
    searchExtend.maxCharacters = 0;
    
    int64_t range;
    
    for(int i = start + length - 1; i >= start; i--) {
	// Go from the end of our selected region to the beginning.

	Log::info() << "On position " << i << " from " <<
	    start + length - 1 << " to " << start << std::endl;

	if(search.positions.size() == 1 && search.positions.front().first.isEmpty()) {
	    Log::info() << "Starting over by mapping position " << i << std::endl;
	    // We do not currently have a non-empty FMDPosition to extend. Start
	    // over by mapping this character by itself.
	    search = this->misMatchMapPosition(rangeIterator, query, i, minContext,
		z_max, maskIterator);

	    if(search.is_mapped && search.characters >= minContext && 
		!search.positions.front().first.isEmpty(maskIterator) && range != -1
		&& search.positions.size() == 1) {

		// It mapped. We didn't do a re-start and fail, we have sufficient
		// context to be confident, and our interval is nonempty and
		// subsumed by a range.
		
		range = search.positions.front().first.range(rangeIterator, maskIterator);
		
		Log::info() << "Mapped " << search.characters << 
		" context to " << search.positions.front().first << " in range #" << range <<
		std::endl;
	    
		// Remember that this base mapped to this range
		mappings.push_back(std::make_pair(range,searchExtend.maxCharacters - 1));
	    
		// We definitely have a non-empty FMDPosition to continue from
	    } else {
		 // It didn't map. Say it corresponds to no range.
		 mappings.push_back(std::make_pair(-1,0));
		    
		 // Mark that the next iteration will be an extension (if we had
		 // any results this iteration; if not it will just restart)
		 search.is_mapped = true;
	    }
	} else {
	    
	    // The last base either mapped successfully or failed due to multi-
	    // mapping. Try to extend the FMDPosition we have to the left
	    // (backwards) with the next base.
	    
	    // Extend by *only* mismatched bases. Do not extend by the correct base yet.
	    searchExtend = this->misMatchExtend(search, query[i], true, z_max, maskIterator, false, true);
	    
	    // Check if mismatch extension gives you any results. If so, restart. See discussion
	    // of mis-identifying mapped positions in the email thread
	    if(searchExtend.positions.size() > 1 || !searchExtend.positions.front().first.isEmpty()) {
		i++;
		search.positions.clear();
		search.positions.push_back(
		    std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
		search.characters = 0;
		search.maxCharacters = 0;
		
	    } else {
		
		// If no mismatch extension results exist, we can safely extend by the correct base
		// and be assured we are passing forward a complete set of search results
		
		Log::info() << "Extending with position " << i << std::endl;
		
		search = this->misMatchExtend(search, query[i], true, z_max, maskIterator, true, false);
		search.characters++;
		
		// What range index does our current left-side position (the one we just
		// moved) correspond to, if any?
		range = search.positions.front().first.range(rangeIterator, maskIterator);
		
		if(search.is_mapped && search.characters >= minContext && 
		    !search.positions.front().first.isEmpty(maskIterator) && range != -1
		    && search.positions.size() == 1) {
		    
		    // It mapped. We didn't do a re-start and fail, we have sufficient
		    // context to be confident, and our interval is nonempty and
		    // subsumed by a range.
		    
		    Log::info() << "Mapped " << search.characters << 
		    " context to " << search.positions.front().first << " in range #" << range <<
		    std::endl;
		
		
		    // Remember that this base mapped to this range
		    mappings.push_back(std::make_pair(range,searchExtend.characters - 1));
		
		    // We definitely have a non-empty FMDPosition to continue from
		   
		} else {
		
		    if(search.is_mapped && search.positions.front().first.isEmpty(maskIterator)
			&& searchExtend.positions.size() == 1) {
		    
			Log::info() << "Failed at " << searchExtend.positions.front().first << " (" << 
			searchExtend.positions.size() << " mismatch search results for " <<
			searchExtend.characters << " context)." << std::endl;
			// We extended right until we got no results. We need to try
			// this base again, in case we tried with a too-long left
			// context.
		
			Log::info() << "Restarting from here..." << std::endl;
		
			search = searchExtend;
		
			// Move the loop index towards the end we started from (right)
			i++;
		
			// Since the FMDPosition is empty, on the next iteration we will
			// retry this base.

		    } else {
			    
			Log::info() << "Failed at " << search.positions.front().first << " (" << 
			search.positions.size() <<
			" mismatch search results for " << search.characters << " context)." << 
			std::endl;
			
			// It didn't map for some other reason:
			// - It was an initial mapping with too little right context to 
			//   be unique to a range.
			// - It was an initial mapping with a nonexistent right context
			// - It was an extension that was multimapped and still is
			
			// In none of these cases will re-starting from this base help
			// at all. If we just restarted here, we don't want to do it
			// again. If it was multimapped before, it had as much left
			// context as it could take without running out of string or
			// getting no results.
			
			// It didn't map. Say it corresponds to no range.
			mappings.push_back(std::make_pair(-1,0));
			
			// Mark that the next iteration will be an extension (if we had
			// any results this iteration; if not it will just restart)
			search.is_mapped = true;
		    }
		}
	    }
	}
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());
    
    // Get rid of the mask iterator if needed
    if(maskIterator != NULL) {
	delete maskIterator;
    }

    // Give back our answers.
    return mappings;
}

std::vector<std::pair<int64_t,size_t>> FMDIndex::misMatchMap(const BitVector& ranges, 
    const std::string& query, int64_t genome, int minContext, size_t z_max, int start, int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return misMatchMap(ranges, query, genome == -1 ? NULL : genomeMasks[genome], 
        minContext, z_max, start, length);    
}

MisMatchAttemptResults FMDIndex::misMatchMapPosition(BitVectorIterator& ranges, 
    const std::string& pattern, size_t index, size_t minContext, size_t z_max,
    BitVectorIterator* mask) const {
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MisMatchAttemptResults result;
            
    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.positions.push_back(std::pair<FMDPosition,size_t>(this->getCharPosition(pattern[index]),0));
    result.characters = 1;
    result.maxCharacters = 1;
    if(result.positions.front().first.isEmpty(mask)) {

        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        result.is_mapped = true;    
        return result;
    } else if (result.positions.front().first.range(ranges, mask) != -1) {
        // We've already mapped.

        result.is_mapped = true;
        return result;
    }
    
    std::vector<std::pair<FMDPosition,size_t>> found_positions;
                
    for(index++; index < pattern.size(); index++) {
	      
        // Forwards extend with subsequent characters.
      
	// Necessary here to create new result set since we're editing every
	// position of the last one. Is there a way around this? Don't think so...
            
	MisMatchAttemptResults new_result = this->misMatchExtend(result, pattern[index], false, z_max, mask, false, false);
			
	if(new_result.positions.front().first.isEmpty(mask)) {
	    if(result.positions.size() == 1 && result.characters >= minContext) {
		result.is_mapped = true;
		result.characters = result.maxCharacters;
		return result;
	    } else {
		// Last position multimapped but this extension now
		// maps nowhere. Return an empty set of positions
	      
		// Or we don't map anywhere with at least the minimum
		// context
	      
		result.positions.clear();
		result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
		result.is_mapped = false;
		result.characters = 1;
		return result;
	    }
        }
        
        
        if(!result.is_mapped && new_result.positions.front().first.range(ranges, mask) != -1 &&
	  new_result.positions.size() == 1 && new_result.characters >= minContext) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success, and
            // return it.
	    
	    result.positions = new_result.positions;
	    result.characters++;
	    result.maxCharacters++;
	    result.is_mapped = true;
            found_positions = result.positions;      
        } else if(result.is_mapped && new_result.positions.front().first.range(ranges, mask) != -1) {
	
	    result.positions = new_result.positions;
	    result.maxCharacters++;
	    
	} else {
	  
	    // Otherwise, we still map to a plurality of ranges. Record the
	    // extension and loop again.
	    
	    result.positions = new_result.positions;
	    result.characters++;
	    result.maxCharacters++;
	}	  
   }

   if (result.is_mapped) {
      result.positions = found_positions;
   } else {
	
	// If we get here, we ran out of downstream context and still map to
	// multiple ranges. Just give our multi-mapping FMDPosition and unmapped
	// result.
    
	result.positions.clear();
	result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
	result.is_mapped = false;
    }
	
    return result;

}

std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> FMDIndex::CmisMap(const BitVector& ranges,
    const std::string& query, const BitVector* mask, int minContext, size_t z_max, int start, int length) const {
    
    // Map to a range.
    
    if(length == -1) {
        // Fix up the length parameter if it is -1: that means the whole rest of
        // the string.
        length = query.length() - start;
    }
    
    Log::debug() << "Mapping with (two-sided) minimum " << minContext << " context." <<
        std::endl;

    // Make an iterator for ranges, so we can query it.
    BitVectorIterator rangeIterator(ranges);
    
    // And one for the mask, if needed
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
        new BitVectorIterator(*mask);

    // We need a vector to return.
    std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> mappings;

    // Keep around the result that we get from the single-character mapping
    // function. We use it as our working state to track our FMDPosition and how
    // many characters we've extended by. We use the is_mapped flag to indicate
    // whether the current iteration is an extension or a restart.

    MisMatchAttemptResults location;
    
    // Make sure the scratch position is empty so we re-start on the first base
    
    location.positions.clear();
    location.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
    location.characters = 1;
    

    for(int i = start + length - 1; i >= start; i--) {
        // Go from the end of our selected region to the beginning.
	
        Log::debug() << "On position " << i << " from " <<
            start + length - 1 << " to " << start << std::endl;
	    
        location = this->CmisMatchMapPosition(rangeIterator, query, i, minContext, z_max, maskIterator);

        // What range index does our current left-side position (the one we just
        // moved) correspond to, if any?
        int64_t range = location.positions.front().first.range(rangeIterator, maskIterator);

        if(location.is_mapped) {
            
            // It mapped. We didn't do a re-start and fail, we have sufficient
            // context to be confident, and our interval is nonempty and
            // subsumed by a range.

            // Remember that this base mapped to this range
            mappings.push_back(std::make_pair(range,std::make_pair(location.characters,location.maxCharacters)));
            
            // We definitely have a non-empty FMDPosition to continue from

        } else {

                mappings.push_back(std::make_pair(-1,std::make_pair(0,0)));

        }
    }

    // We've gone through and attempted the whole string. Put our results in the
    // same order as the string, instead of the backwards order we got them in.
    // See <http://www.cplusplus.com/reference/algorithm/reverse/>
    std::reverse(mappings.begin(), mappings.end());

    // Get rid of the mask iterator if needed
    if(maskIterator != NULL) {
        delete maskIterator;
    }

    // Give back our answers.
    return mappings;
    
}

std::vector<std::pair<int64_t,std::pair<size_t,size_t>>> FMDIndex::CmisMap(const BitVector& ranges, 
    const std::string& query, int64_t genome, int minContext, size_t z_max, int start, int length) const {
    
    // Get the appropriate mask, or NULL if given the special all-genomes value.
    return CmisMap(ranges, query, genome == -1 ? NULL : genomeMasks[genome], z_max, 
        minContext, start, length);    
}

MisMatchAttemptResults FMDIndex::CmisMatchMapPosition(BitVectorIterator& ranges, 
	const std::string& pattern, size_t index, size_t z_max, size_t minContext, BitVectorIterator* mask) const {
    
    // We're going to right-map so ranges match up with the things we can map to
    // (downstream contexts)

    // Initialize the struct we will use to return our somewhat complex result.
    // Contains the FMDPosition (which we work in), an is_mapped flag, and a
    // variable counting the number of extensions made to the FMDPosition.
    MisMatchAttemptResults result;

    // Do a forward search.
    // Start at the given index, and get the starting range for that character.
    result.is_mapped = false;
    result.positions.push_back(std::pair<FMDPosition,size_t>(this->getCharPosition(pattern[index]),0));    
    result.characters = 1;
    result.maxCharacters = 1;
    if(result.positions.front().first.isEmpty(mask)) {
        // This character isn't even in it. Just return the result with an empty
        // FMDPosition; the next character we want to map is going to have to
        // deal with having some never-before-seen character right upstream of
        // it.
        return result;
    } else if (result.positions.front().first.range(ranges, mask) != -1) {
        // We've already mapped.
        result.is_mapped = true;
    }
    
    std::vector<std::pair<FMDPosition,size_t>> found_positions;
    MisMatchAttemptResults new_result;
    MisMatchAttemptResults new_result2;
    MisMatchAttemptResults rightFirstResults;
    MisMatchAttemptResults leftFirstResults;
    
    for(size_t i = 1; index + i < pattern.size() && 1 + index > i; i++) {
		
        // Dual extend with subsequent characters.
      
	new_result2 = this->misMatchExtend(result, pattern[index + i], false, z_max, mask, false);
	
	if(new_result2.positions.front().first.isEmpty(mask)) {
	    if(result.positions.size() == 1 && result.maxCharacters >= minContext) {
		result.is_mapped = true;
		result.characters = result.maxCharacters;
		return result;
	    } else {
		result.positions.clear();
		result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
		result.is_mapped = false;
		result.characters = 1;
		result.maxCharacters = 1;
		return result;
	    }
        }
	
	new_result = this->misMatchExtend(new_result2, pattern[index - i], true, z_max, mask);
	
        if(new_result.positions.front().first.isEmpty(mask)) {
	    if(result.positions.size() == 1 && result.maxCharacters >= minContext) {
		result.is_mapped = true;
		result.characters = result.maxCharacters;
		return result;
	    } else {
		result.positions.clear();
		result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
		result.is_mapped = false;
		result.characters = 1;
		result.maxCharacters = 1;
		return result;
	    }
        }

        if(!result.is_mapped && new_result.positions.size() == 1 &&
	    new_result.positions.front().first.range(ranges, mask) != -1
	    && result.maxCharacters >= minContext) {
            // We have successfully mapped to exactly one range. Update our
            // result to reflect the additional extension and our success
    
    	    result.positions = new_result.positions;
	    result.maxCharacters++;
	    result.characters = result.maxCharacters;
	    result.is_mapped = true;
	    found_positions = result.positions;
	    
        } else if(result.is_mapped && new_result.positions.front().first.range(ranges, mask) != -1) {
	    result.positions = new_result.positions;
	    result.maxCharacters++;

	} else {
	    // Otherwise, we still map to a plurality of ranges. Record the
	    // extension and loop again.
	
	    result.positions = new_result.positions;
	    result.maxCharacters++;
	    result.characters = result.maxCharacters;
	    
	}
    }
    
    if (result.is_mapped) {
      result.positions = found_positions;
   } else {
	// If we get here, we ran out of downstream context and still map to
	// multiple ranges. Just give our multi-mapping FMDPosition and unmapped
	// result.
    
	result.positions.clear();
	result.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
	result.is_mapped = false;
	result.characters = 1;
	result.maxCharacters = 1;
    }

    // If we get here, we ran out of downstream context and still map to
    // multiple ranges. Just give our multi-mapping FMDPosition and unmapped
    // result.
    return result;
}

// In case anyone wants it later... the following two functions implement a
// mismatch extend which returns results sorted by number of mismatches

MisMatchAttemptResults FMDIndex::sortedMisMatchExtend(MisMatchAttemptResults& prevMisMatches,
	char c, bool backward, size_t z_max, BitVectorIterator* mask) const {
    MisMatchAttemptResults nextMisMatches;
    nextMisMatches.is_mapped = false;
    nextMisMatches.characters = prevMisMatches.characters;
    
    // Note that we do not flip parameters when !backward since
    // FMDIndex::misMatchExtend uses FMDIndex::extend which performs
    // this step itself
    
    if(prevMisMatches.positions.size() == 0) {
	throw std::runtime_error("Tried to extend zero length mismatch vector");
    }
    
    if(prevMisMatches.positions.front().first.isEmpty(mask)) {
	throw std::runtime_error("Can't extend an empty position");
    }

    if(c == '\0') {
        throw std::runtime_error("Can't extend with null byte!");
    }

    if(!isBase(c)) {
        std::string errorMessage = std::string("Character #");
        errorMessage.push_back(c);
        errorMessage += std::string(" is not a DNA base.");
        throw std::runtime_error(errorMessage);
    }
    
    // z tracks how many mismatches each range has in order to
    // sort our range-holding data structure. We initialize z with
    // the minimum number of mismatches for which a context existed
    // the previous extension
        
    size_t z = prevMisMatches.positions.front().second;
    
    // To store our FMDPositions under consideration in order of
    // number z of mismatches so that we place them on the queue
    // in the correct order
    
    std::vector<std::pair<FMDPosition,size_t>> waitingMatches;
    std::vector<std::pair<FMDPosition,size_t>> waitingMisMatches;
    
    std::pair<FMDPosition,size_t> m_position;
    std::pair<FMDPosition,size_t> m_position2;
        
    // Output has flag to mark whether we have found a unique
    // hit at the most favourable z-level, which has our
    // minimum context length, in which case we want to terminate
    // search. If not flagged, we want to pass the entire queue
    // as working material for the next extension since we don't
    // know which "mismatch path" will find us our end result
    
        
    for(std::vector<std::pair<FMDPosition,size_t>>::iterator it =
	prevMisMatches.positions.begin(); it != prevMisMatches.positions.end(); ++it) {
	// Store each successive element as an m_position

	m_position.first = it->first;
	m_position.second = it->second;
	
	// Check if we exhausted the search over all sequences in the
	// queue with z mismatches. If so, search all exact-base
	// extensions of level-z sequences to see if there is a single unique
	// hit. Else add all such sequences to the queue, and next search
	// all mismatched extensions for a unique hit.
	
	if(m_position.second != z) {
	  
	    // Built-in check to make sure our data structure holding all extended
	    // context ranges has actually been arranged in mismatch number order
	    // by the previous iteration of processMisMatchPositions
	  
	    if(m_position.second < z) {
		throw std::runtime_error("Generated misordered mismatch list");
	    }
	    
	    // Evaluate and search all extensions 
	  	  	    
	    processMisMatchPositions(nextMisMatches, waitingMatches, waitingMisMatches, mask);
	    
	    // Else we need to see what we get at subsequent z-levels
	    z++;
	    
	}

	// extend m_position by correct base
		
	m_position2.first = extend(m_position.first, c, backward);
	m_position2.second = z;
	waitingMatches.push_back(m_position2);
	
	if(z < z_max) {
	    for(size_t base = 0; base < NUM_BASES; base++) {
		if(BASES[base] != c) {
		    m_position2.first = extend(m_position.first, BASES[base], backward);
		    m_position2.second = z;
		    m_position2.second++;
		    waitingMisMatches.push_back(m_position2);
		
		}
	    }
	}
    }
    
    
    // We have extended and searched the entire queue of matches from the
    // last level
        
    processMisMatchPositions(nextMisMatches, waitingMatches, waitingMisMatches, mask);

    
    // If there are matches, but not unique matches of at least
    // minimum context length, return the entire "heap"
    // generated in this run, to use as starting material for the next
    
    // If there are unique matches, these will also get passed. Don't
    // need the flag 
        
    if(nextMisMatches.positions.size() == 0) {
	nextMisMatches.positions.push_back(std::pair<FMDPosition,size_t>(EMPTY_FMD_POSITION,0));
	nextMisMatches.is_mapped = 0;
    }
            
    return nextMisMatches;
}

void FMDIndex::processMisMatchPositions(
		MisMatchAttemptResults& nextMisMatches,
		std::vector<std::pair<FMDPosition,size_t>>& waitingMatches,
		std::vector<std::pair<FMDPosition,size_t>>& waitingMisMatches,
		BitVectorIterator* mask) const {
		  
    while(!waitingMatches.empty()) {
	if(waitingMatches.back().first.getLength(mask) > 0) {
	    nextMisMatches.positions.push_back(waitingMatches.back());
	}
      
	waitingMatches.pop_back();
    }
  
    while(!waitingMisMatches.empty()) {
	if(waitingMisMatches.back().first.getLength(mask) > 0) {
	    nextMisMatches.positions.push_back(waitingMisMatches.back());
	}
      
	waitingMisMatches.pop_back();
  }
  
  return;
  }