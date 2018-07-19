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
    """Get genomewide nucleotide frequencies"""
    nuCounts = Counter(sequence)
    nucTotal = len(sequence)

    for i in nuCounts:
        print(" %s %s%%" % (i, round(nuCounts[i]/nucTotal * 100)))


def getSeq(fastaFile, repOrigins, length):
    genome = pysam.Fastafile(fastaFile)

    chroms = []
    pos = []
    seqs = []

    coordinates = defaultdict(dict)
    with open(repOrigins, 'r') as origins:
        for l in origins:
            parts = l.rstrip().split('\t')
            chrom = str(parts[0])
            start = int(parts[1])
            end = int(parts[2])

            # midpoint = int((end - start)/2)
            # start = int(start + midpoint)
            start = end
            end = end + length
            window = genome.fetch(chrom, start, end)

            coordinates[chrom][start] = window

            seqs.append(window)
            chroms.append(chrom)
            pos.append(start)

        return(list(zip(chroms, pos, seqs)))


def forkStaller(seqInfo, poolsize, length):
    dNTPS = dict.fromkeys(['A', 'T', 'C', 'G'], poolsize)

    nFreq = defaultdict(int)
    replicatedLength = 0


    for i in range(length):
        q = 0
        random.shuffle(seqInfo)
        for c, p, s in seqInfo:
            # c, p, s = value
            if len(s) <= i+1:
                # print("Overshot!. DNA synthesis reached the end of chromosome %s (%s:%s-%s)" % (c, c, p, p+i))
                # seqInfo.pop(index)
                del seqInfo[q]
                continue

            n = s[i:i + 1]

            if n == 'N':
                continue

            nFreq[n] += 1
            dNTPS[n] -= 1
            replicatedLength += 1
            q += 1

            if dNTPS[n] == 0:
                return (c, p, n, i, s, nFreq, replicatedLength, dNTPS)


def runSim(min, max, n, length, seqInfo, fasta, outfile):
    genome = pysam.Fastafile(fasta)
    print("Running %s simulations" % n)
    print("Randomly generating a dNTP poolsize within defined boundaries of %s and %s" % (min, max))
    with open(outfile, 'w') as fastaOut:
        for x in range(n):
            print("-- Iteration %s --" % (x+1))
            pool = random.randint(min, max)
            print ("Using dNTP poolsize of %s" % pool)
            c, p, b, i, s, nFreq, replicatedLength, dNTPS = forkStaller(seqInfo, pool, length)
            forkstallPos = p + i
            print("Starting from a replication origin of %s:%s, d%sTP exhausted and fork stalled %s bases away at %s:%s" % (c, p, b, i, c, forkstallPos))
            print("o Final dNTP counts at forkstall:\n"
                  " dATPs = %s\n"
                  " dTTPs = %s\n"
                  " dCTPs = %s\n"
                  " sGTPs = %s" % (dNTPS['A'], dNTPS['T'], dNTPS['C'], dNTPS['G']))
            print("o Nucleotide frequency in replicated DNA:")
            for nucleotide in nFreq:
                f = round((nFreq[nucleotide] / replicatedLength), 2)
                print(" - %s : %s%% [%s]" % (nucleotide, f, nFreq[nucleotide]))

            print("o Sequence +/- 20bps from forskstall:")
            print(genome.fetch(c, forkstallPos - 20, forkstallPos), b, genome.fetch(c, forkstallPos + 1, (forkstallPos + 1) + 20))
            print("-----------")
            fastaOut.write(">%s:%s_%s_%s\n%s\n" % (c, p, i, forkstallPos, genome.fetch(c, forkstallPos-100, forkstallPos+100)))


def fileLength(file):
    with open(file) as f:
        return len(f.read().split(b'\n')) - 1


def get_args():
    parser = OptionParser()

    parser.add_option("-b",
                      "--bed_file",
                      dest="bed_file",
                      action="store",
                      help="A bed file containing the replication origins",
                      metavar="FILE")

    parser.add_option("-g",
                      "--genome",
                      dest="genome",
                      action="store",
                      help="Genome fasta file",
                      metavar="FILE")

    parser.add_option("--pool_min",
                      dest="min",
                      action="store",
                      type="int",
                      help="Set the minimal dNTP pool size [Default: 1e4]")

    parser.add_option("--pool_max",
                      dest="max",
                      action="store",
                      type="int",
                      help="Set the maximum dNTP pool size [Default: 1e6]")

    parser.add_option("-n",
                      "--sim_num",
                      dest="sim_num",
                      action="store",
                      type="int",
                      help="Number of simulations to run")

    parser.add_option("-s",
                      "--seed",
                      dest="seed",
                      action="store",
                      type="int",
                      help="Seed to use")

    parser.add_option("-o", \
                      "--outfile", \
                      dest="outfile",
                      action="store",
                      help="File to write sequences surrounding fork stall positions to")

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
    length = int(options.max/originCount)*4

    info = getSeq(options.genome, options.bed_file, length)

    random.seed(42)

    runSim(min=options.min, max=options.max, n=options.sim_num, length=length, seqInfo=info, fasta=options.genome, outfile=options.outfile)


if __name__ == "__main__":
    sys.exit(main())