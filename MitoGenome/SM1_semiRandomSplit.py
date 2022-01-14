#!/usr/bin/env python
''' Random Split of sequence reads.
    Author: Pedro Tarroso
    Date: 28/06/2017
    Version 1.0

    This script accepts as input an gzip compressed fastaq file with interleaved
    paired-end sequence reads and outputs multiple fq.gz files with a random 
    split of reads. 
    The user define the number of files to be generated with random reads and
    the script will throw each read to a random file.

    Typical usage:
        semiRandomSplit.py INPUT_FASTA_FILE number_of_random_splits 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import gzip
import random
import sys

def randorder(size):
    x = range(size)
    random.shuffle(x)
    return(x)

def prepFiles(gzfile, splits):
    basename = ".".join(gzfile.split(".")[:-2])
    infile = gzip.open(gzfile, "r")
    outfiles = []
    for i in range(splits):
        partfmt = "{bn}_part{num:0{width}}.fastq"
        partfile = partfmt.format(bn=basename, num=i+1, width=len(str(splits)))
        outfiles.append(open(partfile, "w"))
    return([infile, outfiles])

def closeFiles(files):
    infile = files[0]
    outfiles = files[1]
    infile.close()
    for f in outfiles:
        f.close()

def getPairedSeqs(instream, nseqs=1):
    seqs = []
    for s in range(nseqs):
        paired = []
        for p in range(8):
            line = instream.readline()
            paired.append(line)
        seqs.append(paired)
    if seqs == [['' for x in range(8)] for y in range(nseqs)]:
        raise EOFError
    return(seqs)

def randomSeqsWrite(files, splits):
    end = False
    while not end:
        try:
            seqs = getPairedSeqs(files[0], splits)
            order = randorder(splits)
            for i in range(splits):
                files[1][i].writelines(seqs[order[i]])
        except EOFError:
            end = True
    closeFiles(files)


if __name__ == "__main__":

    inFile = sys.argv[1]
    nParts = int(sys.argv[2])

    files = prepFiles(inFile, nParts)
    randomSeqsWrite(files, nParts)

            



