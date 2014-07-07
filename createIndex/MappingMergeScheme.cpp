#include <stdexcept>

#include <Log.hpp>
#include <Util.h> // From libsuffixtools, for reverse_complement

#include "MappingMergeScheme.hpp"



MappingMergeScheme::MappingMergeScheme(const FMDIndex& index, 
    const BitVector& rangeVector, 
    const std::vector<std::pair<std::pair<size_t, size_t>, bool> >& rangeBases, 
    const BitVector& includedPositions, size_t genome, size_t minContext) : 
    MergeScheme(index), threads(), queue(NULL), rangeVector(rangeVector), 
    rangeBases(rangeBases), includedPositions(includedPositions), 
    genome(genome), minContext(minContext) {
    
    // Nothing to do
    
}

MappingMergeScheme::~MappingMergeScheme() {
    
    join(); // Join our threads

    if(queue != NULL) {
        // Get rid of the queue if we made one.
        delete queue;
        queue = NULL;
    }
    
}

ConcurrentQueue<Merge>& MappingMergeScheme::run() {
    
    if(queue != NULL) {
        // Don't let people call this twice.
        throw std::runtime_error(
            "Called run() twice on a MappingMergeScheme.");
    }
    
    // Grab the limits of the contig range belonging to this genome.
    auto genomeContigs = index.getGenomeContigs(genome);
    
    // Figure out how many threads to run. Let's do one per contig in this
    // genome.
    size_t threadCount = genomeContigs.second - genomeContigs.first;
    
    Log::info() << "Running Mapping merge on " << threadCount << " threads" <<
        std::endl;
    
    // Make the queue    
    queue = new ConcurrentQueue<Merge>(threadCount);
    
    for(size_t contig = genomeContigs.first; contig < genomeContigs.second; 
        contig++) {
        
        // Start up a thread for every contig in this genome.
        threads.push_back(std::thread(&MappingMergeScheme::generateMerges,
                this, contig));
    }

    // Return a reference to the queue.
    return *queue;
    
}

void MappingMergeScheme::join() {

    for(std::vector<std::thread>::iterator i = threads.begin(); 
        i != threads.end(); ++i) {
    
        // Join each thread.
        (*i).join();    
        
    }
    
    // Probably not a good idea to join the same threads twice.
    threads.clear();
    
}

void MappingMergeScheme::generateMerge(size_t queryContig, size_t queryBase, 
    size_t referenceContig, size_t referenceBase, bool orientation) const {
    
    if(referenceBase > index.getContigLength(referenceContig) || 
        referenceBase == 0) {
        
        // Complain that we're trying to talk about an off-the-contig position.
        throw std::runtime_error("Reference base " + 
            std::to_string(referenceBase) + 
            " on contig " + std::to_string(referenceContig) + 
            " is beyond 1-based bounds of length " + 
            std::to_string(index.getContigLength(referenceContig)));
    }
    
    if(queryBase > index.getContigLength(queryContig) || 
        queryBase == 0) {
        
        // Complain that we're trying to talk about an off-the-contig position.
        throw std::runtime_error("Query base " + std::to_string(queryBase) + 
            " on contig " + std::to_string(queryContig) + 
            " is beyond 1-based bounds of length " + 
            std::to_string(index.getContigLength(queryContig)));
    }
    
    // Where in the query do we want to come from? Always on the forward strand.
    // Correct offset to 0-based.
    TextPosition queryPos(queryContig * 2, queryBase - 1);
    
    // Where in the reference do we want to go to? May b e on the forward or
    // reverse strand. Correct text and offset for orientation, and correct
    // offset to 0-based (at the same time).
    TextPosition referencePos(referenceContig * 2 + orientation, 
        orientation ? (index.getContigLength(referenceContig) - referenceBase) :
        (referenceBase - 1));
    
    // Make a Merge between these positions.
    Merge merge(queryPos, referencePos);
        
    // Send that merge to the queue.
    // Lock the queue.
    auto lock = queue->lock();
    // Spend our lock to add something to it.
    queue->enqueue(merge, lock);
    
}

void MappingMergeScheme::generateMerges(size_t queryContig) const {
    
    // What's our thread name?
    std::string threadName = "T" + std::to_string(genome) + "." + 
        std::to_string(queryContig);
        
    // How many bases have we mapped or not mapped
    size_t mappedBases = 0;
    size_t unmappedBases = 0;
    size_t nonReversibleMappings = 0;
    
    // Grab the contig as a string
    std::string contig = index.displayContig(queryContig);
    
    // How many positions are available to map to?
    Log::info() << threadName << " mapping " << contig.size() << 
        " bases via " << BitVectorIterator(includedPositions).rank(
        includedPositions.getSize()) << " bottom-level positions" << std::endl;
	
    // Make an iterator for ranges, so we can query it.
    BitVectorIterator rangeIterator(rangeVector);
    
    // And one for the mask, if needed
    BitVectorIterator* maskIterator = (mask == NULL) ? NULL : 
        new BitVectorIterator(&includedPositions);
    
    // Map it on the right
    std::vector<std::pair<int64_t,std::pair<bool,bool>>> rightMappings = index.map(rangeVector, contig, 
        genome, &includedPositions, minContext);
    
    // Map it on the left
    std::vector<std::pair<int64_t,std::pair<bool,bool>>> leftMappings = index.map(rangeVector, 
        reverseComplement(contig), genome, &includedPositions, minContext);

    
    // Flip the left mappings back into the original order. They should stay as
    // other-side ranges.
    std::reverse(leftMappings.begin(), leftMappings.end());
    
    for(size_t i = 0; i < leftMappings.size(); i++) {
        // For each position, look at the mappings.
	
	if(leftMappings[i].second.first) {
            // We have a left mapping. Grab its base.
            auto leftBase = rangeBases[leftMappings[i].first];
            
            if(rightMappings[i].second.first) {
                // We have a right mapping. Grab its base too.
                auto rightBase = rangeBases[rightMappings[i].first];
                
                // Compare the position (contig, base) pairs (first) and the
                // orientation flags (second)
                if(leftBase.first == rightBase.first && 
                    leftBase.second != rightBase.second) {
		    // These are opposite faces of the same base.
		  
		    if(leftMappings[i].second.second || rightMappings[i].second.second) {
		      // Nothing else maps to our target base
		      generateMerge(queryContig, i + 1, leftBase.first.first, 
                        leftBase.first.second, !leftBase.second);
		      
		    } else {
			// Other positions in the query genome map to our target base
			// Invoke method to handle multiple query bases
		      
			//ie we won't map back at all

		    }
		  
		} else {
		  
		    // We have a multimap with different target sides. We
		    // want to pass this to the multimap handling procedure
		    // below
		  
		    if(leftMappings[i].second.second || rightMappings[i].second.second) {
			//TODO segergate i
			//TODO segregate left and right targets
		    } else {
			unmappedBases++;
				      
		    }
		  
		}
		
	    } else {
		    
		    //We have a left map only
	      
		    if(leftMappings[i].second.second) {
			// Nothing else maps to our target base
			generateMerge(queryContig, i + 1, leftBase.first.first, 
			    leftBase.first.second, !leftBase.second);
		      
		    } else {
			// Other positions in the query genome map to our target base
			// Invoke method to handle multiple query bases
		      
		    }	      
	    }
	    
        } else if(rightMappings[i].second.first) {
	    // We have a right mapping only. Grab its base.
	    auto rightBase = rangeBases[rightMappings[i].first];
	    
		 if(leftMappings[i].second.second) {
		    // Nothing else maps to our target base
		    generateMerge(queryContig, i + 1, leftBase.first.first, 
                    leftBase.first.second, !leftBase.second);
		      
		 } else {
		    // Other positions in the query genome map to our target base
		    // Invoke method to handle multiple query bases
		      
		 }
	  
	} else if(leftMappings[i].second.second || rightMappings[i].second.second) {
	    // We have a multimapping to some set of bases where these
	    // same bases map back only to position i. Non-unique query
	    // mappings will occur in the other direction
	  
	    //TODO: segregate i
	  
	    if(leftMappings[i].second.second) {
		//TODO: segregate all left targets of i
	      
	    }
	      
	    if(rightMappings[i].second.second) {
		//TODO: segergate all right targets of i
	    }
	  
	    //TODO: do we need to handle sorting multimapping types at this level?
	    //This is a math question!
	  
	  
	} else {
	  
            // Didn't map this one
            unmappedBases++;
        }
        
    }
    // Close the queue to say we're done.
    auto lock = queue->lock();
    queue->close(lock);
    
    // Report that we're done.
    Log::info() << threadName << " finished (" << mappedBases << "|" << 
        unmappedBases << ")" << std::endl;
	
}















