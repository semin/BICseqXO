package org.bioholic

import scala.math._
import scala.io.Source._
import scala.concurrent.stm._
import scala.collection.mutable._
import scala.util.matching.Regex
import scala.util.Sorting.quickSort
import scala.collection.JavaConversions._
import scala.collection.parallel.ForkJoinTasks.defaultForkJoinPool._

import scopt._
import net.sf.samtools._
import java.util.{Date}
import java.io.{File, FileWriter}
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

case class Window(posStart: Int, posEnd: Int, var readCount: Int = 1) {
  def length = posEnd - posStart + 1
  override def toString = "Window[%d-%d:%d]".format(posStart, posEnd, readCount)
}

case class Bin( refName: String, posStart: Int, posEnd: Int,
                var readCountTest: Int = 0,
                var readCountControl: Int = 0,
                val gcContent: Int = 0) {
  def length = posEnd - posStart + 1
  def posRange = (posStart to posEnd)
  def readCountBoth = readCountTest + readCountControl
  def readCountRatio = readCountTest / readCountControl.toDouble
  def readCountRatioTest = readCountTest / readCountBoth.toDouble
  def log2Ratio = 
    if (readCountRatio != 0) {
      log(readCountRatio) / log(2)
    } else {
      log(readCountRatio + 0.0000000001) / log(2)
    }
  def logL =
    if (readCountTest != 0 && readCountControl != 0) {
      readCountTest * log(readCountRatioTest) + readCountControl * log(1.0 - readCountRatioTest)
    } else {
      0
    }
  def contains(w: Window) = posRange.contains(w.posStart) && posRange.contains(w.posEnd)
  def + (operand: Bin): Bin = {
    val newGCContent = (gcContent * length + operand.gcContent * operand.length) / (length + operand.length)
    Bin(refName, posStart, operand.posEnd,
        readCountTest + operand.readCountTest,
        readCountControl + operand.readCountControl, 
        newGCContent)
  }
  override def toString = "Bin[%s:%d-%d:%d/%d:%d]".format(refName, posStart, posEnd, readCountTest, readCountControl, gcContent)
}

object BICseqXO {
  def printMainUsage(): Unit = {
    println("Usage: java -jar BICseqXO-<version>.jar <command> [options]\n")
    println("BICseqXO commands are:")
    println("   map     Create map files containing start positions of uniquely mapped reads")
    println("   bin     Create bin files using a fixed window size or bins in a BED file")
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
      intOpt("q", "minmapq", "<integer>", "mapping quality filtering cutoff (default: 1)", {v: Int => config.minMappingQuality = v})
      intOpt("m", "maxmismatch", "<integer>", "alignment match cutoff for filtering reads (default: 4)", {v: Int => config.maxMismatch = v})
      opt("r", "references", "list of reference names to parse out (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store results (default: '.')", {v: String => config.outDir = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: all)", {v: Int => config.ncpu = v})
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

  def createBinsFromBEDFile(refName: String, bedFile: String, wins1bpTest: Array[Window], wins1bpControl: Array[Window]): Array[Bin] = {
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

  def create1bpWindowsFromMapFile(mapFile: String): Array[Window] = {
    var windows = new ArrayBuffer[Window]()
    var positions = fromFile(mapFile).getLines.map(_.toInt).toArray
    quickSort(positions)
    windows += Window(positions.head, positions.head)

    for (pos <- positions.tail) {
      if (windows.last.posStart == pos) {
        windows.last.readCount += 1
      } else{
        windows += Window(pos, pos)
      }
    }
    windows.toArray
  }

  def createBinFiles(args: Array[String]) = {
    var config = BinConfig()
    val parser = new OptionParser("java -jar BICseqXO-<version>.jar bin") {
      opt("r", "references", "list of reference names to create bins (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store bin files (default: '.')", {v: String => config.outDir = v})
      opt("b", "bed", "<BED file>", "bins in BED format", {v: String => config.bedFile = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: all)", {v: Int => config.ncpu = v})
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
        val wins1bpTest = create1bpWindowsFromMapFile(mapFileTest)
        val wins1bpControl = create1bpWindowsFromMapFile(mapFileControl)
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
      // When arguments are bad, usage message will have been displayed
    }
  }
    
  def createSegFile(binFile: String, segFile: String, lambda: Int, id: String): Unit = {

    def calculateBIC(bins: Array[Bin], lambda: Int, totalReadCount: Int): Double = {
      (bins.size + 1) * lambda * log(totalReadCount) - 2 * bins.map(_.logL).sum
    }

    def diffBIC(bins1: Array[Bin], bins2: Array[Bin], lambda: Int, totalReadCount: Int): Double = {
      calculateBIC(bins1, lambda, totalReadCount) - calculateBIC(bins2, lambda, totalReadCount)
    }

    var tmpBins = new ArrayBuffer[Bin]()
    val refName = FilenameUtils.getBaseName(binFile).split('.').head
    val BinRegex = """^(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+).*$""".r

    for (line <- fromFile(binFile).getLines) {
      line match {
        case BinRegex(readCountTest, readCountBoth, readCountRatioTest, posStart, posEnd) =>
          val readCountControl = readCountBoth.toInt - readCountTest.toInt
          tmpBins += Bin(refName, posStart.toInt, posEnd.toInt, readCountTest.toInt, readCountControl)
        case _ => // do nothing
      }
    }

    Logger("%d bins detected from %s".format(tmpBins.size, binFile))

    // BIC-seq
    var bins = tmpBins.toArray
    val totalReadCount = bins.map(_.readCountBoth).sum
    var originalBIC = calculateBIC(bins, lambda, totalReadCount)
    var noConsecutiveBins = 2

    do {
      var anyBinMerged = false
      var noMergedBinInLevel = 0
      do {
        var minDiffBIC = Ref(0.0)
        var minDiffIdx = Ref(0)
        var minDiffBin = Ref(bins.head) // could be any Bin object
        val untilIdx = bins.size - noConsecutiveBins

        (0 to untilIdx).par.foreach { fromIdx =>
          val binRange = fromIdx to (fromIdx + noConsecutiveBins - 1)
          val mergedBin = binRange.map(bins(_)).reduceLeft(_ + _)
          val binRangeLogL = bins.slice(binRange.head, binRange.last + 1).map(_.logL).sum
          val bicDiff = 2 * (binRangeLogL - mergedBin.logL) - (noConsecutiveBins - 1) * lambda * log(totalReadCount)

          atomic { implicit txn =>
            if (bicDiff < minDiffBIC()) {
              minDiffBIC() = bicDiff
              minDiffIdx() = fromIdx
              minDiffBin() = mergedBin
            }
          }
        }

        atomic { implicit txn =>
          if (minDiffBIC() < 0) {
            originalBIC += minDiffBIC()
            anyBinMerged = true
            noMergedBinInLevel += 1
            val binsBeforeBreak = bins.splitAt(minDiffIdx())._1
            val binsAfterBreak = bins.splitAt(minDiffIdx() + noConsecutiveBins)._2
            bins = binsBeforeBreak ++ Array(minDiffBin()) ++ binsAfterBreak
          } else {
            anyBinMerged = false
          }
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

  def createSegFiles(args: Array[String]): Unit = {
    var config = SegConfig()
    val parser = new OptionParser("java -jar BICseqXO-<version>.jar seg") {
      intOpt("l", "lambda", "<integer>", "penalty term of Bayesian Infromation Criterion (BIC) (default: 2)", {v: Int => config.lambda = v})
      opt("r", "references", "list of reference names to generate segments (default: chr1-22,chrX,chrY)", {v: String => config.refNames = v})
      opt("o", "outdir", "<directory>", "directory to store results (default: '.')", {v: String => config.outDir = v})
      opt("i", "id", "ID for the first column of SEG file format (default: 'N/A')", {v: String => config.id = v})
      intOpt("n", "ncpu", "<integer>", "number of cpus to use (default: all)", {v: Int => config.ncpu = v})
      help("h", "help", "display this help and exit")
      arg("<bin directory>", "directory containing bin files", {v: String => config.binDir = v})
    }

    if (parser.parse(args)) {
      setParallelism(config.ncpu)

      val refNames = parseRefNames(config.refNames)
      val segDir = config.outDir
      FileUtils.forceMkdir(new File(segDir))

      for (refName <- refNames.par) {
        //val binFile = FilenameUtils.concat(config.binDir, refName + ".rc_norm.bin")
        val binFile = FilenameUtils.concat(config.binDir, refName + ".bin")
        val segFile = FilenameUtils.concat(segDir, refName + ".seg")
        createSegFile(binFile, segFile, config.lambda, config.id)
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
