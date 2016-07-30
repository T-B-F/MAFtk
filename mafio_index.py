#!/usr/bin/env python

from intervaltree import IntervalTree
from itertools import islice
from Bio.AlignIO import MafIO
import sys

def make_maf_index(pathfiles, output):
    """ create a tabular file to store maf index
    """
    tree = dict()
    with open(output, "w") as outf:
        for c, pathfile in enumerate(pathfiles):
            with open(pathfile) as inf:
                for i, align in enumerate(MafIO.MafIterator(inf)):
                    block = align._block_lines
                    for record in align:
                        start, size = record.annotations["start"], record.annotations["size"]
                        outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record.id, start, start+size, block[0], block[1], pathfile, i))
                        if record.id not in tree:
                            tree[record.id] = IntervalTree()
                        tree[record.id][start: start+size] = (block[0], block[1], pathfile)
    return tree

def create_index(pathindex):
    """ create an interval tree from a index file
    """
    tree = dict()
    with open(pathindex) as inf:
        for line in inf:
            tmp = line.split()
            name, start, stop, block_start, block_end, pathfile, i = line.strip().split("\t")
            start, stop = int(start), int(stop)
            block_start, block_end = int(block_start), int(block_end)
            if name not in tree:
                tree[name] = IntervalTree()
            tree[name][start: stop] = (block_start, block_end, pathfile)
    return tree

def get_alignment(tree, species, start, stop):
    files = dict()
    for iv in tree[species][start: stop]:
        block_start, block_stop, pathfile = iv.data
        # store all blocks related to a file
        files.setdefault(pathfile, []).append((block_start, block_stop))
    alignments = []
    for pathfile in files:
        positions = sorted(files[pathfile])
        #print(positions)
        prev = 0
        prev_stop = -1
        with open(pathfile) as handle:
            for pos_start, pos_stop in positions:
                assert(pos_start > prev_stop) # block cannot overlap, can they?
                size = pos_stop - pos_start
                #print(pos_start, size)
                line = next(islice(handle, pos_start - 1 - prev , pos_start - prev))
                align = next(MafIO.MafIterator(handle))
                #print(align)
                # only keep part of the alignment that we are interest of
                found = False
                for record in align:
                    if record.id == species:
                        found = True
                        break
                if found:
                    start_seq = record.annotations["start"]
                    # count trailing gap
                    cnt_starting_gap = 0
                    if record.seq[0] == "-":
                        i = 0
                        while True:
                            if record.seq[i] == "-":
                                cnt_starting_gap += 1
                            else:
                                break
                            i += 1
                    # get new start of alignment
                    align_start = max(start_seq, start)
                    offset_start = align_start - start_seq + cnt_starting_gap
                    # get new stop of alignment
                    seq_size = record.annotations["size"]
                    ali_size = len(record.seq)
                    offset_stop = offset_start
                    cur_seq = align_start
                    # gaps are not taken into account, only increasing if no gap caracters
                    # to find right stop position
                    while offset_stop < ali_size:
                        if cur_seq == stop:
                            break
                        elif record.seq[offset_stop] != "-":
                            cur_seq += 1
                        offset_stop += 1
                    alignments.append(align[:, offset_start: offset_stop])
                else:
                    print("Unable to find SeqRecord for species {} in alignment:".format(species))
                    print(align)
                prev = pos_start + size
                prev_stop = pos_stop
    return alignments
            
    
if __name__ == "__main__":
    maf_text = """track name=euArc visibility=pack
##maf version=1 scoring=tba.v8 
# tba.v8 (((human chimp) baboon) (mouse rat)) 
                   
a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

"""
    with open("test.maf", "w") as outf:
        outf.write(maf_text)
        
    # create index
    tree1 = make_maf_index(["test.maf"], "test.maf.idx")
    
    # read index
    tree2 = create_index("test.maf.idx")
    
    # check that the two trees are identical
    print(tree1 == tree2)
    print()
    
    # access alignments
    #aligns = get_alignment(tree1, "hg18.chr7",  27578838, 27707221+3)
    #for align in aligns:
        #print(align)
    #print()    
    
    # access alignments
    aligns = get_alignment(tree1, "mm4.chr6",  53215344+5, 53310102+10)
    for align in aligns:
        print(align)
        
    sys.exit(0)
