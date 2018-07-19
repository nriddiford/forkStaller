#!/usr/bin/env python
from __future__ import division
import os, re, sys
import pysam
import random
import string

from collections import defaultdict
from collections import Counter
from optparse import OptionParser


def getNucFreqs(sequence):
    """Get nucleotide frequencies from a DNA sequence"""
    nuCounts = Counter(sequence)
    nucTotal = len(sequence)
    for i in nuCounts:
        print(" %s %s%%" % (i, round(nuCounts[i]/nucTotal * 100)))


def getSeq(fastaFile, bedFile, length):
    """For a bed file of replication origins, flip a coin to decide whether to fire F or R
        and extract genomic sequence length away"""
    genome = pysam.Fastafile(fastaFile)
    chroms = []
    origins = []
    seqs = []
    direction = []
    with open(bedFile, 'r') as repOris:
        for l in repOris:
            ## Decide whether this origin fires forwards or reverse
            ori = random.choice(['F', 'R'])
            parts = l.rstrip().split('\t')
            chrom = str(parts[0])
            start = int(parts[1])
            end = int(parts[2])

            if ori == 'F':
                pos = end
                start = end
                end = end + length
                window = genome.fetch(chrom, start, end)
            else:
                pos = start
                end = start
                start = start - length
                window = genome.fetch(chrom, start, end)
                window = ''.join(reversed(window))

            seqs.append(window)
            chroms.append(chrom)
            origins.append(pos)
            direction.append(ori)
        return(list(zip(chroms, origins, seqs, direction)))


def forkStaller(seqInfo, poolsize, length):
    """For each index position in 0..length shuffled sequences and remove a dNTP from the pool
        Once a pool is depleted, return the genomic coordinates of the fork stall and associated information"""

    dNTPS = dict.fromkeys(['A', 'T', 'C', 'G'], poolsize)

    nFreq = defaultdict(int)
    replicatedLength = 0

    for i in range(length):
        random.shuffle(seqInfo)
        for index, value in enumerate(list(seqInfo)):
            c, p, s, d = value
            if len(s) <= i+1:
                print("Overshot!. DNA synthesis reached the end of chromosome %s (%s:%s-%s)" % (c, c, p, p+i))
                seqInfo.pop(index)
                continue

            n = s[i:i + 1]

            if n == 'N':
                continue

            nFreq[n] += 1
            dNTPS[n] -= 1
            replicatedLength += 1

            if dNTPS[n] == 0:
                return (c, p, n, d, i, s, nFreq, replicatedLength, dNTPS)


def runSim(min, max, n, length, seqInfo, fasta, outfile):
    genome = pysam.Fastafile(fasta)
    print("Running %s simulations" % n)
    print("Randomly generating a dNTP poolsize within an upper and lower boundary of %s and %s" % (min, max))
    with open(outfile, 'w') as fastaOut:
        for x in range(n):
            print("-- Iteration %s --" % (x+1))
            pool = random.randint(min, max)
            print ("Using dNTP poolsize of %s" % pool)
            c, p, b, d, i, s, nFreq, replicatedLength, dNTPS = forkStaller(seqInfo, pool, length)
            if d == 'F':
                forkstallPos = p + i
                ori = '-->'
                location = 'downstream'
            else:
                forkstallPos = p - i
                ori = '<--'
                location = 'upstream'


            print("Starting from a replication origin of %s:%s and travelling %s, d%sTP exhausted and fork stalled %s bases %s (%s:%s)" % (c, p, ori, b, i, location, c, forkstallPos))
            print("o Final dNTP counts at fork stall:\n"
                  " dATPs = %s\n"
                  " dTTPs = %s\n"
                  " dCTPs = %s\n"
                  " sGTPs = %s" % (dNTPS['A'], dNTPS['T'], dNTPS['C'], dNTPS['G']))
            print("o Nucleotide frequency in replicated DNA:")
            for nucleotide in nFreq:
                f = round((nFreq[nucleotide] / replicatedLength), 2)
                print(" - %s : %s%% [%s]" % (nucleotide, f, nFreq[nucleotide]))

            print("o Sequence +/- 20bps from fork stall:")
            print(genome.fetch(c, forkstallPos - 20, forkstallPos), b, genome.fetch(c, forkstallPos + 1, (forkstallPos + 1) + 20))
            print("-----------")
            # fastaOut.write(">%s:%s_%s_%s\n%s\n" % (c, p, i, forkstallPos, genome.fetch(c, forkstallPos-100, forkstallPos+100)))


def fileLength(file):
    with open(file) as f:
        return len(f.read().split(b'\n')) - 1


def get_args():
    parser = OptionParser()

    parser.add_option("-b",
                      "--bed_file",
                      dest = "bed_file",
                      action = "store",
                      help = "A bed file containing the replication origins",
                      metavar = "FILE")

    parser.add_option("-g",
                      "--genome",
                      dest = "genome",
                      action = "store",
                      help = "Genome fasta file",
                      metavar = "FILE")

    parser.add_option("--pool_min",
                      dest = "min",
                      action = "store",
                      type = "int",
                      help = "Set the minimal dNTP pool size [Default: 1e4]")

    parser.add_option("--pool_max",
                      dest = "max",
                      action = "store",
                      type = "int",
                      help = "Set the maximum dNTP pool size [Default: 1e6]")

    parser.add_option("-n",
                      "--sim_num",
                      dest = "sim_num",
                      action = "store",
                      type = "int",
                      help = "Number of simulations to run")

    parser.add_option("-s",
                      "--seed",
                      dest = "seed",
                      action = "store",
                      type = "int",
                      help = "Seed to use")

    parser.add_option("-o",
                      "--outfile",
                      dest = "outfile",
                      action = "store",
                      help = "File to write sequences surrounding fork stall positions to")

    parser.set_defaults(bed_file = '/Users/Nick_curie/Desktop/misc_bed/RepOris/Bg3origins.mappable.bed',
                        genome = '/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Dmel_6.12.fasta',
                        outfile = 'out/simrepstalls.fa',
                        sim_num = 2,
                        min = int(1e4),
                        max = int(1e6))

    options, args = parser.parse_args()

    return (options, args)


def main():
    options, args = get_args()

    dir = os.path.dirname(options.outfile)
    if not os.path.exists(dir):
        os.makedirs(dir)

    originCount = fileLength(options.bed_file)
    length = (int(options.max/originCount)*4)+2
    info = getSeq(options.genome, options.bed_file, length)

    if options.seed:
        print("Using seed %s" % (options.seed))
        random.seed(options.seed)

    runSim(min = options.min, max = options.max, n = options.sim_num, length = length, seqInfo = info, fasta = options.genome, outfile = options.outfile)


if __name__ == "__main__":
    sys.exit(main())