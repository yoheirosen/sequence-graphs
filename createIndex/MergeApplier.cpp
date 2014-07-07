#include <iostream>

#include <Log.hpp>

#include "MergeApplier.hpp"

MergeApplier::MergeApplier(const FMDIndex& index,
    ConcurrentQueue<Merge>& source, stPinchThreadSet* target): index(index), 
    source(source), target(target), thread(&MergeApplier::run, this) {
    
    // Already started running. See <http://stackoverflow.com/a/10673671/402891>
    
}

void MergeApplier::join() {
    thread.join();
}

void MergeApplier::run() {
    // OK, do the actual merging.
  
    if(flagNonReversible) {
	std::unordered_set<TextPosition> targets;
    }
    
    while(true) {
        // Lock the queue when either there's something in it or all its writers
        // have closed it. If neither happens we just wait here forever.
        auto lock = source.waitForNonemptyOrEnd();
        
        if(source.isEmpty(lock)) {
            // We've finished our job.
            return;
        }
        
	// We are merging all mappings without regard to reversibility
        
	// If we get here, there's actual work to do. Trade our lock for a value
	// to work on.
	Merge merge = source.dequeue(lock);
        
        if(flagNonReversible &&  merge.second) {
	  
	  // We want to segregate bases whose mapping is reversible (and merge
	  // them, from those whose mapping is non-reversible. If the user
	  // specified that we should not merge non-reversible bases
	  
	  // TODO use bitvector to flag non-reversible (vs. just unmapped
	  // positions
        
	} else {
        
	  // Now actually apply the merge.
        
	  // Unpack the first TextPosition to be merged
	  size_t firstContigNumber = index.getContigNumber(merge.first.first);
	  size_t firstStrand = index.getStrand(merge.first.first);
	  size_t firstOffset = index.getOffset(merge.first.first);
        
	  // And the second
	  size_t secondContigNumber = index.getContigNumber(merge.first.second);
	  size_t secondStrand = index.getStrand(merge.first.second);
	  size_t secondOffset = index.getOffset(merge.first.second);
	  
	  if(flagNonReversible){
	      // We pay special attention to the case that a base in the reference
	      // genome is mapped to by distinct left-side unique and right-side
	      // unique bases, since this gives us a multimapping in the reverse
	      // direction which is not caught by our mapping-level checks.
	    
	      if(targets.count(merge.first.second) = 1) {
		  // We found the exceptional case above, you want to flag both
		  // bases, and reverse the previous pinching operation.
		
		  // TODO: make sure you've got the pinching order correct... in
		  // order to not disturb reference genome pinches you should
		  // only be popping the last-added element off the block
		
		  // TODO: implement pop operation for blocks. Right now you
		  // have to destruct then reconstruct the block... Hopefully
		  // these collisions will be rare enough that this
		  // shouldn't be too huge a problem
		  
		  stPinchSegment* unPinchSegment = stPinchThreadSet_getSegment(target, 
			  secondContigNumber);
		  
		  stPinchBlock* unPinchBlock = stPinchSegment_getBlock(unPinchSegment);
		  
		  stPinchBlockIt unPinchBlockIt = stPinchBlock_getSegmentIterator(unPinchBlock);
		  
		  std::vector<stPinchSegment*> rebuildPinchBlock;
		  std::vector<bool> rebuildOrientations;
		  
		  rebuildPinchBlock.push_back(stPinchSegment_getFirst(unPinchBlock));
		  rebuildOrientations.push_back(stPinchSegment_getBlockOrientation(stPinchSegment_getFirst(unPinchBlock)));
		  
		  size_t unPinchBlockCounter =  stPinchBlock_getDegree(unPinchBlock);
		  
		  for(i=1; i < unPinchBlockCounter - 1; i++) {
		      rebuildPinchBlock.push_back(stPinchBlockIt_getNext(unPinchBlockIt));		    
		  }
		  
		  stPinchBlock_destruct(unPinchBlock);
		  
		  for(i=1; i < unPinchBlockCounter; i++) {
		      stPinchThread_pinch(rebuildPinchBlock[0],rebuildPinchBlock[i],rebuildOrientations[0] != rebuildOrientations[i]);
		      
		  }

		  //We have unpinched the whole block, then rebuild without the most recent element
	      }
	    
	  }
        
	  // What orientation should we use for the second strand, given that we
	  // are pinching against the first strand in orientation 1 (reverse)?
	  bool orientation = firstStrand == secondStrand;
	  
	  // Grab the first pinch thread
	  stPinchThread* firstThread = stPinchThreadSet_getThread(target,
	      firstContigNumber);
            
	  // And the second
	  stPinchThread* secondThread = stPinchThreadSet_getThread(target,
	      secondContigNumber);
            
	  // Log the pinch if applicable.
	  Log::trace() << "\tPinching #" << firstContigNumber << ":" <<
	      firstOffset << " strand " << firstStrand << " and #" << 
	      secondContigNumber << ":" << secondOffset << " strand " << 
	      secondStrand << " (orientation: " << orientation << ")" <<
	      std::endl;
        
	  // Perform the pinch
	  stPinchThread_pinch(firstThread, secondThread, firstOffset,
	      secondOffset, 1, orientation);
        
	}
    }
    
}


