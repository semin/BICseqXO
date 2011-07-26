package org.bioholic

import scala.math._
import scala.BigInt._
import scala.io.Source._
import scala.concurrent.stm._
import scala.collection.mutable._
import scala.util.matching.Regex
import scala.util.Sorting.quickSort
import scala.collection.JavaConversions._
import scala.collection.parallel.ForkJoinTasks.defaultForkJoinPool._

import scopt._
import java.util.{Date}
import net.sf.samtools._
import java.io.{File, FileWriter}
import com.google.common.collect.TreeMultimap
import org.apache.commons.io.{FileUtils, FilenameUtils}

object Logger {
  def currentTime: String = { 
    val date = new Date()
    val time = date.getTime()
    "%s %s".format(date, time)
  }
  def info(msg: String) = println("[%s]: %s".format(currentTime, msg))
  def apply(msg: String) = info(msg)
}

case class Command(command: String)

case class MapConfig(
  var minMappingQuality: Int  = 1,
  var maxMismatch: Int = 4,
  var refNames: String = "chr1-22,chrX,chrY",
  var outDir: String = ".",
  var ncpu: Int = Runtime.getRuntime().availableProcessors(),
  var bamTest: String = null,
  var bamControl: String = null
)

case class BinConfig(
  var refNames: String = "chr1-22,chrX,chrY",
  var bedFile: String = "",
  var outDir: String = ".",
  var ncpu: Int = Runtime.getRuntime().availableProcessors(),
  var mapDirTest: String = null,
  var mapDirControl: String = null
)

case class SegConfig(
  var refNames: String = "chr1-22,chrX,chrY",
  var lambda: Int = 2,
  var binDir: String = ".",
  var outDir: String = ".",
  var id: String = "N/A",
  var ncpu: Int = Runtime.getRuntime().availableProcessors()
)

case class Win(posStart: Int, posEnd: Int, var readCount: Int = 1) {
  def length = posEnd - posStart + 1
  override def toString = "Win[%d-%d:%d]".format(posStart, posEnd, readCount)
}

case class Bin( refName: String, posStart: Int, posEnd: Int,
                var readCountTest: Int = 0,
                var readCountControl: Int = 0,
                val gcContent: Int = 0) {
  val length = posEnd - posStart + 1
  val posRange = (posStart to posEnd)
  val readCountBoth = readCountTest + readCountControl
  val readCountRatio = readCountTest / readCountControl.toDouble
  val readCountRatioTest = readCountTest / readCountBoth.toDouble
  lazy val log2Ratio =
    if (readCountRatio != 0)
      log(readCountRatio) / log(2)
    else
      log(readCountRatio + 0.0000000001) / log(2)
  lazy val logL =
    if (readCountTest != 0 && readCountControl != 0)
      readCountTest * log(readCountRatioTest) + readCountControl * log(1.0 - readCountRatioTest)
    else
      0
  def contains(w: Win) = posRange.contains(w.posStart) && posRange.contains(w.posEnd)
  def + (operand: Bin): Bin = {
    val newGCContent = (gcContent * length + operand.gcContent * operand.length) / (length + operand.length)
    Bin(refName, posStart, operand.posEnd,
        readCountTest + operand.readCountTest,
        readCountControl + operand.readCountControl, 
        newGCContent)
  }
  override def toString = "Bin[%s:%d-%d:%d/%d:%d]".format(refName, posStart, posEnd, readCountTest, readCountControl, gcContent)
}

case class SegBin(tumor: Int, total: Int, freq: Double, start:Int, end: Int) {
  def length = end - start + 1
  def range = (start to end)
  def normal = total - tumor
  def log2Ratio =
    if (freq != 0)
      log(freq) / log(2)
    else
      log(freq + 0.000000000001) / log(2)
  def logL =
    if (tumor != 0 && tumor != total)
      tumor * log(freq) + (total - tumor) * log(1.0 - freq)
    else
      0.0
  def + (other: SegBin): SegBin = {
    val _tumor = tumor + other.tumor
    val _total = total + other.total
    val _freq = _tumor / _total.toDouble
    val _start = List(start, other.start).min
    val _end = List(end, other.end).max
    SegBin(_tumor, _total, _freq, _start, _end)
  }
  override def toString = "Bin[%d-%d:%d:%d]".format(start, end, tumor, total)
}

object BICseqXO {
  def printMainUsage(): Unit = {
    println("\nUsage: java -jar BICseqXO-<version>.jar <command> [options]\n")
    println("BICseqXO commands are:")
    println("   map     Create map files containing start positions of uniquely mapped reads")
    println("   bin     Create bin files using a fixed window size or coordinates in a BED file")
    println("   seg     Create seg files by merging bins based on Bayesian Information Criterion (BIC)\n")
    println("See 'java -jar BICseqXO-<version>.jar <commnad> -h/--help' for more information on a specific command.")
  }

  def parseRefNames(arg: String): Array[String] = {
    val refNames = arg.split(',')
    var newRefnames = new ArrayBuffer[String]()
    val RefNameRegex = """(\w*)(\d+)-(\d+)""".r
    for (refName <- refNames) {
      refName match {
        case RefNameRegex(prefix, start, end) => (start.toInt to end.toInt).foreach(r => newRefnames += (prefix + r.toString))
        case _ => newRefnames += refName
      }
    }
    newRefnames.toArray
  }

  def createMapFile(bamFile: String, refName: String, mapDir: String, minMappingQuality: Int, maxMismatch: Int): Unit = {
    val mapFile = FilenameUtils.concat(mapDir, refName + ".map")
    val mapFileWriter = new FileWriter(mapFile)
    val bamFileReader = new SAMFileReader(new File(bamFile))
    val reads = bamFileReader.query(refName, 0, 0, false)

    for (read <- reads) {
      if (read.getMappingQuality >= minMappingQuality) {
        var matCount = 0
        val matRegex = """(\d+)M""".r
        for (matRegex(mat) <- matRegex.findAllIn(read.getCigarString)) {
          matCount += mat.toInt
        }
        if ((read.getReadLength() - matCount) <= maxMismatch) {
          mapFileWriter.write(read.getAlignmentStart + "\n")
        }
      }
    }
    reads.close
    bamFileReader.close
    mapFileWriter.close
  }

  def createMapFiles(args: Array[String]): Unit = {
    SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT)

    var config = MapConfig()
      val parser = new OptionParser("java -jar BICseqXO-<version>.jar map") {
      intOpt("q", "minmapq", "<integer>", "mapping quality cutoff for filtering reads (default: 1)", {v: Int => config.minMappingQuality = v})
      intOpt("m", "maxmismatch", "<integer>", "alignment mismatch cutoff for filtering reads (default: 4)", {v: Int => config.maxMismatch = v})
      opt("r", "references", "list of reference names to parse out (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store results (default: '.')", {v: String => config.outDir = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: number of all available cpus)", {v: Int => config.ncpu = v})
      help("h", "help", "display this help and exit")
      arg("<test BAM file>", "BAM file for test (tumor) sample", {v: String => config.bamTest = v})
      arg("<control BAM file>", "BAM file for control (normal) sample", {v: String => config.bamControl = v})
    }

    if (parser.parse(args)) {
      setParallelism(config.ncpu)
      val refNames = parseRefNames(config.refNames)
      for (bam <- List(config.bamTest, config.bamControl)) {
        val mapDir = FilenameUtils.concat(config.outDir, FilenameUtils.getBaseName(bam))
        FileUtils.forceMkdir(new File(mapDir))
        refNames.par.foreach { r => createMapFile(bam, r, mapDir, config.minMappingQuality, config.maxMismatch) }
      }
    } else {
      // When arguments are bad, usage message will have been displayed
    }
  }

  def createBinsFromBEDFile(refName: String, bedFile: String, wins1bpTest: Array[Win], wins1bpControl: Array[Win]): Array[Bin] = {
    var bins = new ArrayBuffer[Bin]()
    val BinRegex = """^(\S+)\s+(\d+)\s+(\d+)\s+(\d*).*$""".r

    for (line <- fromFile(bedFile).getLines()) {
      line match {
        case BinRegex(chrName, posStart, posEnd, gcContent) =>
          if (chrName == refName) {
            if (gcContent.size > 0)
              bins += Bin(chrName, posStart.toInt + 1, posEnd.toInt, gcContent=gcContent.toInt)
            else 
              bins += Bin(chrName, posStart.toInt + 1, posEnd.toInt)
          }
        case _ => // do nothing
      }
    }
    // probably sorting process here
    Logger("%d bins in %s detected from %s".format(bins.size, refName, bedFile))

    var tIdx = 0
    var cIdx = 0

    for (b <- bins) {
      if (wins1bpTest.isDefinedAt(tIdx) && (wins1bpTest(tIdx).posStart <= b.posEnd)) {
        while (wins1bpTest.isDefinedAt(tIdx) && !b.contains(wins1bpTest(tIdx)) && (wins1bpTest(tIdx).posStart <= b.posEnd)) {
          tIdx += 1
        }
        while (wins1bpTest.isDefinedAt(tIdx) && b.contains(wins1bpTest(tIdx))) {
          b.readCountTest += wins1bpTest(tIdx).readCount
          tIdx += 1
        }
      }
      if (wins1bpControl.isDefinedAt(cIdx) && (wins1bpControl(cIdx).posStart <= b.posEnd)) {
        while (wins1bpControl.isDefinedAt(cIdx) && !b.contains(wins1bpControl(cIdx)) && (wins1bpControl(cIdx).posStart <= b.posEnd)) {
          cIdx += 1
        }
        while (wins1bpControl.isDefinedAt(cIdx) && b.contains(wins1bpControl(cIdx))) {
          b.readCountControl += wins1bpControl(cIdx).readCount
          cIdx += 1
        }
      }
    }
    bins.toArray
  }

  def create1bpWinsFromMapFile(mapFile: String): Array[Win] = {
    var windows = new ArrayBuffer[Win]()
    var positions = fromFile(mapFile).getLines.map(_.toInt).toArray
    quickSort(positions)
    windows += Win(positions.head, positions.head)

    for (pos <- positions.tail) {
      if (windows.last.posStart == pos) {
        windows.last.readCount += 1
      } else{
        windows += Win(pos, pos)
      }
    }
    windows.toArray
  }

  def createBinFiles(args: Array[String]) = {
    var config = BinConfig()
    val parser = new OptionParser("java -jar BICseqXO-<version>.jar bin") {
      opt("r", "references", "list of reference names to create bins (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store bin files (default: '.')", {v: String => config.outDir = v})
      opt("b", "bed", "<BED file>", "coordinates in BED format", {v: String => config.bedFile = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: number of all available cpus)", {v: Int => config.ncpu = v})
      help("h", "help", "display this help and exit")
      arg("<test map directory>", "directory containing map files for test sample", {v: String => config.mapDirTest = v})
      arg("<control map directory>", "directory containing map files for control sample", {v: String => config.mapDirControl = v})
    }

    if (parser.parse(args)) {
      setParallelism(config.ncpu)
      val refNames = parseRefNames(config.refNames)
      val refNameToBins = TMap[String, ArrayBuffer[Bin]]().single
      var totalReadCountTest = Ref(0)
      var totalReadCountControl = Ref(0)
      val binDir = config.outDir
      FileUtils.forceMkdir(new File(binDir))

      for (ref <- refNames.par) {
        val mapFileTest = FilenameUtils.concat(config.mapDirTest, ref + ".map")
        val mapFileControl = FilenameUtils.concat(config.mapDirControl,  ref + ".map")
        val wins1bpTest = create1bpWinsFromMapFile(mapFileTest)
        val wins1bpControl = create1bpWinsFromMapFile(mapFileControl)
        var totalRefReadCountTest = 0
        var totalRefReadCountControl = 0

        atomic { implicit txn => refNameToBins += ref -> new ArrayBuffer[Bin]() }

        if (config.bedFile != "") {
          val bedFile = new File(config.bedFile)
          if (bedFile.exists) {
            val binFile = FilenameUtils.concat(binDir, ref + ".raw.bin")
            val binFileWriter = new FileWriter(binFile)
            val bins = createBinsFromBEDFile(ref, config.bedFile, wins1bpTest, wins1bpControl)
            for (bin <- bins) {
              totalRefReadCountTest += bin.readCountTest
              totalRefReadCountControl += bin.readCountControl
              atomic { implicit txn => refNameToBins(ref) += bin }
              if (bin.readCountBoth != 0) {
                binFileWriter.write("%d\t%d\t%.6f\t%d\t%d\n".format(bin.readCountTest, bin.readCountBoth, bin.readCountRatioTest, bin.posStart, bin.posEnd))
              }
            }
            atomic { implicit txn =>
              totalReadCountTest() = totalReadCountTest() + totalRefReadCountTest
              totalReadCountControl() = totalReadCountControl() + totalRefReadCountControl
            }
            binFileWriter.close
            Logger("%d/%d (test/control) reads detected from %s".format(totalRefReadCountTest, totalRefReadCountControl, ref))
          } else {
            // exit program for the moment.
            System.err.println("Cannot find %s".format(config.bedFile))
            System.exit(1)
          }
        } else {
          // do regular BIC-seq without bins
          System.err.println("You must provide a BED file for the current version of BICseqXO.")
          System.exit(1)
        }
      }

      // normalization by total read count ratio
      var scale = totalReadCountTest.single.get / totalReadCountControl.single.get.toDouble
      Logger("Scale factor for read count normalization: %.6f (Test: %d / Control: %d)".format(scale, totalReadCountTest.single.get, totalReadCountControl.single.get))

      for (ref <- refNames.par) {
        val binFile = FilenameUtils.concat(binDir, ref + ".rc_norm.bin")
        val binFileWriter = new FileWriter(binFile)
        for (bin <- refNameToBins(ref)) {
          val readCountTestNorm = (scale * bin.readCountTest).ceil.toInt
          if (bin.readCountBoth != 0) {
            binFileWriter.write("%d\t%d\t%.6f\t%d\t%d\n".format(readCountTestNorm, bin.readCountBoth, bin.readCountRatioTest, bin.posStart, bin.posEnd))
          }
        }
        binFileWriter.close
        Logger("Read counts of %s normalized by scale factor, %.6f".format(ref, scale))
      }
    } else {
       //When arguments are bad, usage message will have been displayed
    }
  }
    
  def createSegFile(binFile: String, segFile: String, lambda: Int, id: String): Unit = {

    def calculateBIC(bins: ArrayBuffer[Bin], lambda: Int, totalReadCount: Int): Double =
      (bins.size + 1) * lambda * log(totalReadCount) - 2 * bins.par.map(_.logL).sum

    def calculateBICDiff(bins1: ArrayBuffer[Bin], bins2: ArrayBuffer[Bin], lambda: Int, totalReadCount: Int): Double =
      calculateBIC(bins1, lambda, totalReadCount) - calculateBIC(bins2, lambda, totalReadCount)

    var bins = new ArrayBuffer[Bin]()
    val refName = FilenameUtils.getBaseName(binFile).split('.').head
    val BinRegex = """^(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+).*$""".r

    for (line <- fromFile(binFile).getLines) {
      line match {
        case BinRegex(readCountTest, readCountBoth, readCountRatioTest, posStart, posEnd) =>
          val readCountControl = readCountBoth.toInt - readCountTest.toInt
          bins += Bin(refName, posStart.toInt, posEnd.toInt, readCountTest.toInt, readCountControl)
        case _ => // do nothing
      }
    }

    Logger("%d bins detected from %s".format(bins.size, binFile))

    // BIC-seq
    val totalReadCount = bins.par.map(_.readCountBoth).sum
    var originalBIC = calculateBIC(bins, lambda, totalReadCount)
    var noConsecutiveBins = 2
    var binPosRangeToBin = HashMap[(Int, Int), Bin]()
    var binPosRangeToDiffBIC = HashMap[(Int, Int), Double]()

    do {
      var anyBinMerged = false
      var noMergedBinInLevel = 0
      do {
        var minDiffBIC = 0.0
        var minDiffIdx = 0
        var minDiffBin = bins.head // could be any Bin object
        val untilIdx = bins.size - noConsecutiveBins

        (0 to untilIdx).foreach { fromIdx =>
          val binPosRange = (bins(fromIdx).posStart, bins(fromIdx + noConsecutiveBins - 1).posEnd)
          val (diffBIC, mergedBin) = {
            if (binPosRangeToBin.contains(binPosRange)) {
              (binPosRangeToDiffBIC(binPosRange), binPosRangeToBin(binPosRange))
            } else {
              val binsToMerge = bins.slice(fromIdx, fromIdx + noConsecutiveBins)
              val binsToMergeLogL = binsToMerge.map(_.logL).sum
              val _mergedBin = binsToMerge.reduceLeft(_ + _)
              val _diffBIC = 2 * (binsToMergeLogL - _mergedBin.logL) - (noConsecutiveBins - 1) * lambda * log(totalReadCount)
              binPosRangeToBin(binPosRange) = _mergedBin
              binPosRangeToDiffBIC(binPosRange) = _diffBIC
              (_diffBIC, _mergedBin)
            }
          }

          if (diffBIC < minDiffBIC) {
            minDiffBIC = diffBIC
            minDiffIdx = fromIdx
            minDiffBin = mergedBin
          }
        }

        if (minDiffBIC < 0) {
          anyBinMerged = true
          noMergedBinInLevel += 1
          originalBIC += minDiffBIC
          bins(minDiffIdx) = minDiffBin
          bins.remove(minDiffIdx + 1, noConsecutiveBins - 1)
        } else {
          anyBinMerged = false
        }
      } while (anyBinMerged)

      Logger("Level %d, Merged bins: %d, BIC: %.6f, Bin size: %d".format(noConsecutiveBins, noMergedBinInLevel, originalBIC, bins.size))

      if (noMergedBinInLevel > 0)
        noConsecutiveBins += 1
      else
        noConsecutiveBins -= 1
    } while (noConsecutiveBins > 1)

    // write segmentaion results in SEG (CBS) format
    val segFileWriter = new FileWriter(segFile)
    segFileWriter.write("'ID\tChromosome\tStart\tEnd\tRead_Count_Tumor\tRead_Count_Normal\tLog2Ratio\n")
    bins.foreach { b => segFileWriter.write("%s\t%s\t%d\t%d\t%d\t%d\t%.6f\n".format(id, b.refName, b.posStart, b.posEnd, b.readCountTest, b.readCountControl, b.log2Ratio)) }
    segFileWriter.close
  }

  def createSegFileOptimized(binFile: String, segFile: String, lambda: Int, id: String = "N/A") = {

    def calculateBIC(bins: ArrayBuffer[SegBin], lambda: Int, totalReadCount: Int): Double =
      (bins.size + 1) * lambda * log(totalReadCount) - 2 * bins.par.map(_.logL).sum

    def calculateBICDiff(bins1: ArrayBuffer[SegBin], bins2: ArrayBuffer[SegBin], lambda: Int, totalReadCount: Int): Double =
      calculateBIC(bins1, lambda, totalReadCount) - calculateBIC(bins2, lambda, totalReadCount)

    val bins = new ArrayBuffer[SegBin]()
    val refName = FilenameUtils.getBaseName(binFile).split('.').head
    val BinRegex = """^(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+).*$""".r

    for (line <- fromFile(binFile).getLines) {
      line match {
        case BinRegex(tumor, total, freq, start, end) =>
          bins += SegBin(tumor.toInt, total.toInt, freq.toDouble, start.toInt, end.toInt)
        case _ => // do nothing
      }
    }

    val totalReadCount = bins.map(_.total).sum
    Logger("Initial bins: %d, lamba: %d, BIC: %.6f".format(bins.size, lambda, calculateBIC(bins, lambda, totalReadCount)))

    case class SegWin(bins: ArrayBuffer[SegBin]) {
      val mergedBin = bins.reduceLeft(_+_)
      val bicDiff = calculateBICDiff(ArrayBuffer(mergedBin), bins, lambda, totalReadCount)
    }

    class DoubleOrdering extends Ordering[Double] {
      override def compare(x: Double, y: Double): Int = Ordering.Double.compare(x, y)
    }

    class SegWinOrdering extends Ordering[SegWin] {
      override def compare(x: SegWin, y: SegWin): Int = Ordering.Int.compare(x.bins(0).start, y.bins(0).start)
    }

    def runBicSeqLevel(noConsecutiveBins: Int): Int = {
      //Logger("\t%d-bins: total bins: %d".format(noConsecutiveBins, bins.size))
      //Logger("\t%d-bins: initial BIC: %.6f".format(noConsecutiveBins, calculateBIC(bins, lambda, totalReadCount)))

      val wins = new ArrayBuffer[SegWin]()
      val bicDiffToWin: TreeMultimap[Double, SegWin] = TreeMultimap.create(new DoubleOrdering, new SegWinOrdering)

      var i = 0
      while (i <= bins.size - noConsecutiveBins) {
        val w = SegWin(bins.slice(i, i + noConsecutiveBins))
        wins += w
        bicDiffToWin.put(w.bicDiff, w)
        i += 1
      }

      def mergeWin(w: SegWin) = {
        val binIdx = bins.indexOf(w.bins(0))
        bins(binIdx) = w.mergedBin
        bins.remove(binIdx + 1, noConsecutiveBins - 1)

        val winIdx = wins.indexOf(w)
        val fromIdx = List(0, winIdx - noConsecutiveBins + 1).max
        val toIdx = List(wins.size - 1, winIdx+noConsecutiveBins - 1).min

        var j = 0
        while (j <= toIdx - fromIdx) {
          bicDiffToWin.remove(wins(fromIdx).bicDiff, wins(fromIdx))
          wins.remove(fromIdx)
          j += 1
        }

        var k = winIdx
        while (k >= fromIdx) {
          val w = SegWin(bins.slice(k, k + noConsecutiveBins))
          wins.insert(fromIdx, w)
          bicDiffToWin.put(w.bicDiff, w)
          k -= 1
        }
      }

      var anyBinMerged = false
      var noMergedBins = 0

      do {
        try {
          val minKey = bicDiffToWin.asMap().firstKey()
          if (minKey < 0) {
            mergeWin(bicDiffToWin.get(minKey).first())
            noMergedBins += 1
            anyBinMerged = true
          } else
            anyBinMerged = false
        } catch {
          case e: NoSuchElementException =>
            anyBinMerged = false
        }
      } while (anyBinMerged)

      return noMergedBins
    }

    var noMergedBins = 0
    var noConsecutiveBins = 2

    do {
      noMergedBins = runBicSeqLevel(noConsecutiveBins)
      Logger("Level %d, Merged bins: %d, BIC: %.6f, Bin size: %d".format(noConsecutiveBins, noMergedBins, calculateBIC(bins, lambda, totalReadCount), bins.size))
      //Logger("Merged %d bins of size %d; %d remaining.".format(noMergedBins, noConsecutiveBins, bins.size))

      if (noMergedBins > 0)
        noConsecutiveBins += 1
      else
        noConsecutiveBins -= 1
    } while (noConsecutiveBins > 1)

    // write segmentaion results in SEG (CBS) format
    val segFileWriter = new FileWriter(segFile)
    segFileWriter.write("'ID\tChromosome\tStart\tEnd\tRead_Count_Tumor\tRead_Count_Normal\tLog2Ratio\n")
    bins.foreach { b => segFileWriter.write("%s\t%s\t%d\t%d\t%d\t%d\t%.6f\n".format(id, refName, b.start, b.end, b.tumor, b.normal, b.log2Ratio)) }
    segFileWriter.close
  }

  def createSegFiles(args: Array[String]): Unit = {
    var config = SegConfig()
    val parser = new OptionParser("java -jar BICseqXO-<version>.jar seg") {
      intOpt("l", "lambda", "<integer>", "penalty term of Bayesian Infromation Criterion (BIC) (default: 2)", {v: Int => config.lambda = v})
      opt("r", "references", "list of reference names to generate segments (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store results (default: '.')", {v: String => config.outDir = v})
      opt("i", "id", "ID for the first column of SEG file format (default: 'N/A')", {v: String => config.id = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: number of all available cpus)", {v: Int => config.ncpu = v})
      help("h", "help", "display this help and exit")
      arg("<bin directory>", "directory containing bin files", {v: String => config.binDir = v})
    }

    if (parser.parse(args)) {
      setParallelism(config.ncpu)
      val refNames = parseRefNames(config.refNames)
      val segDir = config.outDir
      FileUtils.forceMkdir(new File(segDir))

      for (refName <- refNames) {
        val binFile = FilenameUtils.concat(config.binDir, refName + ".bin")
        val segFile = FilenameUtils.concat(segDir, refName + ".seg")
        //createSegFile(binFile, segFile, config.lambda, config.id)
        //createSegFileOptimized2(binFile, segFile, config.lambda, config.id)
        createSegFileOptimized(binFile, segFile, config.lambda, config.id)
      }
    } else {
      // When arguments are bad, usage message will have been displayed
    }
  }

  def main(args: Array[String]) = {
    if (args.size == 0) {
      printMainUsage()
      System.exit(1)
    } else {
      val timeStart = System.nanoTime: Double
      val command = Command(args.head)
      command match {
        case Command("map") => createMapFiles(args.tail)
        case Command("bin") => createBinFiles(args.tail)
        case Command("seg") => createSegFiles(args.tail)
        case _ => {
            System.err.println("Error: '%s' is not a BICseqXO command.\n".format(args.head))
            printMainUsage()
            System.exit(1)
        }
      }
      val timeEnd = System.nanoTime: Double
      val timeElapsed = (timeEnd - timeStart) / 1e9
      Logger("Elapsed time: %.6f s".format(timeElapsed))
      System.exit(0)
    }
  }
}
