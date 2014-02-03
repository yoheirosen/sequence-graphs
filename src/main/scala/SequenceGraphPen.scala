package edu.ucsc.genome

/**
 * Provides a "pen" that can be used to "draw" Sequence Graph elements (Sides,
 * AlleleGroups, Adjacencies, and Anchors).
 *
 * Exists mostly as a workaround for the fact that the Sequence Graph object
 * constructors as generated by Avro are awkward: if you want to do anything
 * other than specify all the fields, you need to instantiate a new Builder just
 * to build your object, and everything that needs a unique ID demands that you
 * provide it with one.
 *
 * We want to replace this:
 * {{{
 * // Make the Anchor
 * val anchor = new Anchor(
 *     new Edge(nextID, startSide.id, endSide.id), 
 *     new PloidyBounds(1, null, null), 
 *     sample)
 * nextID += 1
 * }}}
 * with this:
 * {{{
 * pen.drawAnchor(startSide, endSide)
 * }}}
 *
 * We have some implicit conversions in the package (like Side to its Long ID)
 * that make this a bit easier.
 *
 * The pen, when constructed, takes the base ID, the number of IDs it is allowed
 * to generate (before infringing on some other parallel task's ID space), and
 * the genome name it is drawing.
 */
class SequenceGraphPen(idStart: Long, idCount: Long, genome: String) {
    // Keep track of the next ID available to allocate
    var nextID = idStart
    
    /**
     * Allocate and return an ID. If we have run out of IDs, raise an exception.
     */
    protected def id: Long = {
        // Grab the next ID value
        val toReturn = nextID
        // Say the next next ID is 1 after that one
        nextID += 1
        
        if(toReturn >= idStart + idCount) {
            // Complain we ran out of IDs
            throw new Exception(
                "Pen can't allocate %d from ID block of %d at %d".format(
                toReturn, idCount, idStart))
        }
        
        // Return the ID we generated
        toReturn
    }
    
    /**
     * Produce a new Side on the given contig, with the given base and face.
     * Optionally can flag the Side as reference; by default it is not.
     */
    def drawSide(contig: String, base: Long, face: Face,
        isReference: Boolean = true): Side = {
        
        new Side(id, new Position(contig, base, face), isReference)
    }
    
    /**
     * Make a new Edge between the Sides with the given IDs. This is not the
     * same as adjacencies; drawAdjacency draws those.
     */
    def drawEdge(left: Long, right: Long): Edge = new Edge(id, left, right)
    
    /**
     * Make a new AlleleGroup between the given pair of Side IDs, carrying the
     * given Allele and using the given PloidyBounds (which defaults to a
     * constant ploidy of 1).
     */
    def drawAlleleGroup(left: Long, right: Long, allele: Allele, 
        ploidy: PloidyBounds = 1): AlleleGroup = {
        
        new AlleleGroup(drawEdge(left, right), allele, ploidy, genome)
    }
    
    /**
     * Make a new Adjacency between the given Sides, with the given ploidy
     * bounds (defaulting to a constant 1).
     */
    def drawAdjacency(left: Long, right: Long, 
        ploidy: PloidyBounds = 1): Adjacency = {
        
        new Adjacency(drawEdge(left, right), ploidy, genome)
    }
    
    /**
     * Makes a new Anchor between the given Sides, with the given ploidy bounds
     * (defaulting to a constant 1).
     */
    def drawAnchor(left: Long, right: Long, 
        ploidy: PloidyBounds = 1): Anchor = {
        
        new Anchor(drawEdge(left, right), ploidy, genome)
    }
    
    
}

