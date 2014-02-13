package edu.ucsc.genome

import scala.util.parsing.combinator.RegexParsers

/**
 * Represents an interval in a linear contig. Uses 1-based inclusive indexing
 * and parses/outputs contig:start-end notation. Semantically, runs from the
 * left Side of the start base to the right Side of the end base. Start must be
 * less than or equal to end. Contig may not contain ":".
 */
class BaseRange(val contig: String, val start: Long, val end: Long) {
    
    if(start > end) {
        // We can't have backwards bounds
        throw new Exception("Can't create a negative-length range %s:%d-%d"
            .format(contig, start, end))
    }
    
    if(contig.contains(":")) {
        // We can't have pathological contig names.
        throw new Exception(
            "Cannot make a range on %s because contig name contains \":\""
            .format(contig))
    }
    
    def this(string: String) {
        // If only I could define some vals here this would be much simpler. And
        // we could have error checking. Everything would be great.
        this(string.split(":")(0), string.split(":")(1).split("-")(0).toLong, 
            string.split(":")(1).split("-")(1).toLong)
    }
    
    /**
     * Test if the given Position is within this range.
     */
    def contains(position: Position): Boolean = {
        // It's in range if it's on the correct contig, and it's between start
        // and end inclusive. Face doesn't matter since we implicitly represent
        // the outer faces of the bases we end at.
        position.contig == contig && position.base >= start && 
            position.base <= end
    }
    
    /**
     * Output in contig:start-end format.
     */
    override def toString = "%s:%d-%d".format(contig, start, end)    
    
}