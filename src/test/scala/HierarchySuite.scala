package edu.ucsc.genome

import org.scalatest._
import java.io._
import java.nio.file._
import org.apache.commons.io._
import scala.sys.process._

/**
 * A test utility class that provides a createIndex-generated hierarchy fixture
 * on disk.
 */
abstract class HierarchySuite extends FunSuite with BeforeAndAfterAll {

    // Holds the temporary directory
    var scratch: Path = null
    
    // Holds the directory of the index
    var indexName: String = null
    
    /**
     * This can be overridden to specify the sequences to index.
     */
    def sequences = Seq("ACTAGT")
    
    /**
     * This can be overridden to specify the context length to use for merging
     * positions.
     */
    def contextLength = 3
    
    override def beforeAll = {
        // Make a new temp directory
        scratch = Files.createTempDirectory("test")
        
        // Decide on the index directory
        indexName = scratch.resolve("index").toString
    
        val files = for((sequence, index) <- sequences.zipWithIndex) yield {
    
            // Pick a name
            val fasta = scratch.resolve("fasta%d".format(index)).toString
            
            // Write a FASTA
            val fastaWriter = new FileWriter(fasta)
            
            // Write the sequence
            fastaWriter.write(">seq%d\n".format(index))
            fastaWriter.write("%s\n".format(sequence))
            
            fastaWriter.close
            
            // Return the file name to be indexed.
            fasta
        }
        
        // Invoke the createIndex tool on all the FASTAs, with the given context
        // length.
        val returnCode = (Seq("./createIndex.sh", "--context", 
            contextLength.toString, indexName) ++ files)!

        if(returnCode != 0) {
            // Don't continue on with the test if we can't construct the
            // fixture.
            throw new Exception("Invoking createIndex produced return code %d"
                .format(returnCode))
        }
        
    }
    
    override def afterAll = {
        // Clean it up
        FileUtils.deleteDirectory(new File(scratch.toString))
    }
    
}