import scala.collection.SeqLike
import perfplot._
import HWCounters.IvyBridge;
import java.io._
import java.util.regex.Pattern
import java.util.regex.Matcher

object AdvancedRoofline {
	def main( args: Array[String] ) {
		def extractName( raw_name:String ) = {
			val pattern = Pattern.compile("bst_\\d+");
			val m = pattern.matcher(raw_name)
			if( m.find() ) m.group().replace("_","") 
			else "bst_000"
		}

		val src = scala.io.Source.fromFile("names_sizes.txt");
		val lines = src.mkString.split("\n");
		src.close
		val raw_impl_map = 
		(for( i <- 0 to lines.size-1 ) yield (i, ((x:Array[String]) => (extractName(x(0)), x(1).toLong)) (lines(i).split(" ")))).groupBy( _._2._1 )


		def getSizeFromTo( l: Seq[(Int, (String, Long))] ) = {
			val sortedLines = l.sortWith(_._1 > _._1)
			val sizes = sortedLines map (_._2._2)
			val max = sortedLines.head._1
			val min = sortedLines.last._1
			(min, max, sizes)
		}

		val impl_map = raw_impl_map map ( (x) => (x._1, getSizeFromTo(x._2))  )

		val counters = CommandService.Counters.apply(new File("."))

		impl_map foreach ( (t) => { val k = t._1; val v = t._2; run_kernel( counters, v._3, k, v._1, v._2 ) } )

		print( (for( i <- 0 to 10 ) yield counters.getFlops(i)) mkString "," )
	}

  def run_kernel (counters: CommandService.Counters, sizes: Seq[Long], impl_name: String, From: Int, To: Int)
  {
		val mask = List(1,1,2,4)

		val outFilenames = "flop" :: "tsc" :: "size" :: "bytes_transferred" :: "Counter3" :: "bytes_read" :: "bytes_write" :: Nil
		val outFiles = Map( (outFilenames map ( (x) => x -> (new PrintStream(x+"_"+impl_name+".txt")))):_* );


		def getFromTo ( f: (Int) => Long ) = for( i <- From to To ) yield f(i)

		for( i <- From to To ) {print( counters.getCounters(i) ); print("\n")}

		val outVals =
		getFromTo( (x) => (counters.getCounters(x),mask).zipped map(_*_) sum) ::
		getFromTo( counters.getTSC(_) ) ::
		sizes ::
		getFromTo( counters.getbytes_transferred(_) )::
    getFromTo( counters.getSCounter3.apply(_) )::
    getFromTo( counters.getbytes_read(_) )::
    getFromTo( counters.getbytes_write(_) )::
		Nil;
		
		print(outVals)
		val outStrings = outVals map ( _.mkString(" ") )

		(outFilenames zip outStrings) foreach ( (x) => (outFiles.get(x._1) getOrElse( System.out )).print(x._2) )

		outFiles foreach ( (x) => x._2.close )

  }
}
