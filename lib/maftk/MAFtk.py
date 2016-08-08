# Copyright 2016 by Tristan Bitard-Feildel. All rights reserved
"""MAFtk module is a set of function to be used with the MafIO parser
provided in biopython (supported in the branch MafIO of 
https://github.com/T-B-F/biopython)

This module contains read and write function to create a interval tree
based on the genomic position of the aligned blocks in MAF files.
A function is also provided to query the tree for a given species,
returning the alignments found from the query's start and stop.
The returned aligned are sliced according to the query's position.

"""

from __future__ import print_function

import sys
from itertools import islice

try:
    from intervaltree import IntervalTree
except:
    print("intervaltree not installed")
    print("please see lib/intervaltree or used pip install intervaltree")
    sys.exit(1)
    
try:
    from Bio.AlignIO import MafIO
except:
    print("MafIO not available")
    print("please installed MafIO branch from "
          "https://github.com/T-B-F/biopython")
    sys.exit(1)

msa_characters = ["-", "?", "!", "*", "."]

def compute_offset_pos(seq, pos):
    """ compute the offset from a position in a MSA to the normal sequence
    
    Parameters
    ==========
    seq : string
        the sequence from the MSA 
    pos : int
        the position to convert
        
    Return
    ======
    k : int
        the position in the MSA
    """

    k = 0 
    cnt = 0 if seq[k] not in msa_characters else -1
    while cnt != pos and k+1 < len(seq):
        k += 1 
        if seq[k] not in msa_characters:
            cnt += 1
    return k


def seq2msa_startstop(seq, start, stop):
    """ Compute new starting positions, taking into account gaps and other 
    non amino acids characters
    
    Parameters
    ==========
    seq : string
        the sequence from the MSA
    start : int
        the start index
    stop : int
        the stop index
        
    Return
    ======
    offstart : int
        the new start
    offstop : int
        the new stop
    size : int
        the size of the segment
    """
    size = stop - start
    # need to add 1 to start, converting list indexes to positions
    offstart = compute_offset_pos(seq, start)
    offstop = compute_offset_pos(seq, stop)
    return offstart, offstop+1, size

class MafTK(object):
    """ toolkit class to handle file in MAF (Multiple Alignment Format
    https://genome.ucsc.edu/FAQ/FAQformat.html#format5
    """
    
    def __init__(self):
        self.tree = dict()

    def make_maf_index(self, pathfiles, output):
        """ create a tabular file to store maf index from a list of files
        """
        if self.tree:
            print("Warning, tree already filled, used MAFtk.clean() "
                " before creating a new index or use a new instance")
            return
        
        with open(output, "w") as outf:
            for c, pathfile in enumerate(pathfiles):
                with open(pathfile) as inf:
                    for i, align in enumerate(MafIO.MafIterator(inf)):
                        block = align._block_lines
                        for record in align:
                            start, size = record.annotations["start"], record.annotations["size"]
                            outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id, start, start+size, block[0], block[1], pathfile))
                            if record.id not in self.tree:
                                self.tree[record.id] = IntervalTree()
                            self.tree[record.id][start: start+size] = (block[0], block[1], pathfile)
        return self.tree

    def read_index(self, pathindex):
        """ read an interval tree from a index file
        """
        if self.tree:
            print("Warning, tree already filled, used MAFtk.clean() "
                " before reading a new index or use a new instance")
            return
        
        with open(pathindex) as inf:
            for line in inf:
                tmp = line.split()
                name, start, stop, block_start, block_end, pathfile = line.strip().split("\t")
                start, stop = int(start), int(stop)
                block_start, block_end = int(block_start), int(block_end)
                if name not in self.tree:
                    self.tree[name] = IntervalTree()
                self.tree[name][start: stop] = (block_start, block_end, pathfile)
        return self.tree
    
    def write_tree(self, output):
        """ store tree into file
        """
        with open(output, "w") as outf:
            for sp in self.tree:
                for iv in self.tree[sp]:
                    # name, start, stop, block start, block end, path
                    start, stop, (block_start, block_end, pathfile) = iv
                    outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, start, stop, block_start, block_end, pathfile))
            
    def get_alignments(self, species, start, stop):
        """ get spanning alignments for a specific range of a given species
        """
        files = dict()
        for iv in self.tree[species][start: stop]:
            block_start, block_stop, pathfile = iv.data
            # store all blocks related to a file
            files.setdefault(pathfile, []).append((block_start, block_stop))
            #print(">", block_start, block_stop)
        alignments = list()
        for pathfile in files:
            positions = sorted(files[pathfile])
            #print(positions)
            prev_pos = 0
            prev_stop = -1
            with open(pathfile) as handle:
                for pos_start, pos_stop in positions:
                    assert(pos_start > prev_stop) # block cannot overlap, can they?
                    pos_size = pos_stop - pos_start
                    # we remove prev because we only move of the number of lines between the two blocks
                    line = next(islice(handle, pos_start - 1 - prev_pos , pos_start - prev_pos)) 
                    align = next(MafIO.MafIterator(handle))
                    #print(align)
                    # look for the sequence of the species
                    found = False
                    for record in align:
                        if record.id == species:
                            found = True
                            break
                    # only keep part of the alignment that we are interest of
                    if found:
                        start_seq = record.annotations["start"]
                        seq_size = record.annotations["size"]                      
                        
                        max_start = max(start_seq, start)
                        max_stop = min(start_seq+seq_size, stop)
                        
                        abs_start = max_start - start_seq
                        abs_stop = max_stop - start_seq
                        # get new start of alignment
                        
                        msa_start, msa_stop, msa_size = seq2msa_startstop(str(record.seq), abs_start, abs_stop)
                        #print(pos_start, pos_stop, start_seq, seq_size, max_start, max_stop, abs_start, abs_stop)
                        #print(msa_start, msa_stop)
                        #offset_start = align_start - start_seq + cnt_starting_gap
                        # get new stop of alignment
                        #ali_size = len(record.seq)
                        alignments.append(align[:, msa_start: msa_stop])
                    else:
                        print("Unable to find SeqRecord for species {} in alignment:".format(species))
                        print(align)
                    prev_pos = pos_start + pos_size # + ali_size # ?
                    prev_stop = pos_stop
        return alignments
    
    def clean(self):
        """ erase tree 
        """
        self.tree = dict()

        