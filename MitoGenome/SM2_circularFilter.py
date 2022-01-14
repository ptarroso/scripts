#!/usr/bin/env python
''' Circular extraction script.
    Author: Pedro Tarroso
    Date: 28/06/2017
    Version 1.0

    A simple script to find and extract the circular part of a dna sequence by
    removing the superfluous sequences at both sides. The algorithm searchs a 
    pattern that is repeated and extracts the sequence within the repeated
    pattern, preserving the pattern at the left side, but removing from the 
    right side. There are two modes of operation: (1) each sequence in the input 
    fasta file is analysed independently, resulting in likely different patterns
    in each sequence or (2) all sequences are searched simultaneously for the 
    best pattern. The second mode might take more time but all sequences will
    begin in the same sequence. Due to the lack of aligment and rotation of 
    sequences, some can be eliminated when using the second method. 

    Typical usage:
        circularFilter.py [-l] [-m] [-h] INPUT_FASTA > OUTPUT_FASTA

    Type "circularFilter.py -h" for help.

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

import re
import argparse

parser = argparse.ArgumentParser(description='Circular extraction')

parser.add_argument('input', help="Input fasta file")
parser.add_argument('--length', '-l', action="store", type=int, default=50,
                    help="Length of the pattern to search (default=50)")
parser.add_argument('--multi', '-m', action="store_true",
                    help="Look for same pattern in all files?")
args = parser.parse_args()

class sequence:
    def __init__(self, name, sequence):
        self.__name = name
        self.__seq = sequence
    def seq(self, seq=None):
        if seq:
            self.__seq = seq
        else:
            return(self.__seq)
    def name(self, name=None):
        if name:
            self.__name = name
        else:
            return(self.__name)
    def fasta(self):
        return(">{0}\n{1}".format(self.__name, self.__seq))

def readSeqs(fastafile):
    seqs = []    
    stream = open(fastafile, "r")
    while True:
        seq = readSeq(stream)
        if not seq:
            break
        seqs.append(seq)
    stream.close()
    return(seqs)

def readSeq(stream):
    header = stream.readline().strip()
    if header != '':
        seq = stream.readline().strip()
        seq = sequence(header[1:], seq)
        return(seq)
    else :
        return(False)

def extractCircle(seq, pattern):
    ''' Given a pattern and a sequence, this function extacts the circle if the
        pattern is found at least twice. If the pattern is found more then two
        times, it seachs for the highest possible length. It writes the number 
        of times the pattern was found in the sequence header
    '''
    dna = seq.seq()
    p = re.compile(pattern)
    f = p.findall(dna)
    result = None
    if len(f) > 1:        
        p2 = re.compile("{0}.*{0}".format(pattern))
        name = "{0}   circular with {1} matches".format(seq.name(), len(f))
        result = sequence(name, p2.findall(dna)[0][:-len(pattern)])
    return(result)

def singleCircle(seq, length=50):
    ''' This algorithm searchs independently each sequence in the fasta file
        for a sequence of user defined length (default 50 nucleotides). The 
        search pattern is defined from the left and moving a single nucleotide
        each iteration to look for the highest length of circular dna.
    '''
    dna = seq.seq()
    best = sequence('','')
    i = 0
    while i < (len(dna)/2):
        pattern = dna[i:(i+length)]
        circ = extractCircle(seq, pattern)
        if circ:
            if len(best.seq()) < len(circ.seq()):
                best = circ
        i += 1
    return(best)

def multiCircle(seqs, length=50):
    ''' Detects circularity by finding the same sequence with user defined
        length. The same pattern is tested in all sequences, so they should
        begin in the same position. The best result might have less sequences
        than the original file. 
        The algorithm can be further improved to avoid testing all combinations.
    '''
    lengths = [len(x.seq()) for x in seqs]
    avglen = sum(lengths)/len(lengths)
    
    results = []
    maxlen = 0
    i = 0
    while i < (avglen/2):
        for seq in seqs:
            dna = seq.seq()
            pattern = dna[i:(i+length)]
            circs = [extractCircle(s, pattern) for s in seqs]
            circlen = [len(x.seq()) for x in circs if x is not None]
            if sum(circlen) > maxlen:
                results = circs
                maxlen = sum(circlen)
        i += 1
        
    results = [x for x in results if x is not None]
    return(results)


if __name__ == '__main__':
    seqs = readSeqs(args.input)
    if args.multi:
        circularSeqs = multiCircle(seqs, args.length)
    else:
        circularSeqs = [singleCircle(seq, args.length) for seq in seqs]

    for seq in circularSeqs:
        print(seq.fasta())
