#!/usr/bin/env python
# Copyright 2016 - Tristan Bitard-Feildel.  All rights reserved.

"""Tests for maftk"""

import unittest
from itertools import islice
import sys, os

from intervaltree import IntervalTree, Interval

from maftk.MAFtk import MafTK

# maf text3 from UCSC Mouse chr10 first two alignment
maf_text1 = """##maf version=1 scoring=autoMZ.v1                                                                                                        
a score=11761.000000                                                                                                                     
s hg19.chr10      60000 45 + 135534747 GAATTCCTTGAGGCCTAAATGCATCGGGGTGCTCTGGTTTTGTTG                                                     
s ponAbe2.chrUn 3569800 45 +  72422247 GAATTCCTTGAGGCCTAATTGCATCAGGGTGCTCTGGTTTTGTTG                                                     
q ponAbe2.chrUn                        999999999999999999999999999999999999999999999                                                     
i ponAbe2.chrUn N 0 C 0                                                                                                                  
s panTro2.chr18    5870 45 +  77261746 GAATTCCTTGAGGCCTAAATGCATCGGGGAGCTCTGGTTTTGTTG                                                     
q panTro2.chr18                        999999888999978999999999889979967847889999999                                                     
i panTro2.chr18 N 0 C 0

a score=45884.000000                                                                                                                     
s hg19.chr10           60045 108 + 135534747 TTGTTATTTCTGAATGACATTTACTTTGGTGCTCTTTATTTTGCGTATTTAAAACTATTAGATCGTGTGATTATATTTGACAGGTCTTAATTGACGCGCTGTTCAGCC                                                                                                                         
s panTro2.chr18         5915 108 +  77261746 TTGTTATTTCTGAATGACATTGACTTTGGTGCTTTTTATTTTGCATATTTAAAACTATTAGATCGTGTGATTATATTTGACAAGTCTTAATTGACGCGCTGTTCAGCG                                                                                                                         
q panTro2.chr18                              999999999999999999999999999995999999999999988999999999999999999963999999999999999999999999999999999999999999                                                                                                                         
i panTro2.chr18      C 0 C 0                                                                                                             
s ponAbe2.chrUn      3569845 108 +  72422247 TTGCTATTTCTGAATGACATTGACTTTGGTGCTCTTTATTTTGCATATTTAAAACTATTAGATCATGTGATTATATTTGACAGGTCTTAATTGATGCACTCTTCAGCG                                                                                                                         
q ponAbe2.chrUn                              999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999                                                                                                                         
i ponAbe2.chrUn      C 0 C 0                                                                                                             
s calJac1.Contig5564   47630 102 -    188978 TTGCTATTTCTGAGTGACATTGACTTTGGTGCACTTTATTTGGAATATTTAAAA-TATTAGATTGTG---TTGTATTTGAAAGGCCTTAATTGATGTGCTGTTCAG--                                                                                                                         
q calJac1.Contig5564                         999999999999999999999999999999999999999999999999999999-999999999999---999999999999999999999999999999999999--                                                                                                                         
i calJac1.Contig5564 N 0 I 14  

##eof maf
"""

# from https://genome.ucsc.edu/FAQ/FAQformat.html#format5
maf_text2 = """track name=euArc visibility=pack
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

test_tree1 = {'calJac1.Contig5564': IntervalTree([Interval(47630, 47732, (10, 21, 'text.maf'))]), 
              'ponAbe2.chrUn': IntervalTree([Interval(3569800, 3569845, (1, 9, 'text.maf')), 
                                             Interval(3569845, 3569953, (10, 21, 'text.maf'))]), 
              'hg19.chr10': IntervalTree([Interval(60000, 60045, (1, 9, 'text.maf')), 
                                          Interval(60045, 60153, (10, 21, 'text.maf'))]), 
              'panTro2.chr18': IntervalTree([Interval(5870, 5915, (1, 9, 'text.maf')), 
                                             Interval(5915, 6023, (10, 21, 'text.maf'))])
              }


test_line = ["hg19.chr10	60000	60045	1	9	text.maf\n",
             "ponAbe2.chrUn	3569800	3569845	1	9	text.maf\n",
             "panTro2.chr18	5870	5915	1	9	text.maf\n",
             "hg19.chr10	60045	60153	10	21	text.maf\n",
             "panTro2.chr18	5915	6023	10	21	text.maf\n",
             "ponAbe2.chrUn	3569845	3569953	10	21	text.maf\n",
             "calJac1.Contig5564	47630	47732	10	21	text.maf\n"]

class TestMafIO(unittest.TestCase):

    def test_one_two(self):    
        
        with open("text.maf", "w") as outf:
            outf.write(maf_text1)
            
        # create index
        maf = MafTK()
        tree1 = maf.make_maf_index(["text.maf"], "text.maf.idx")            
        self.assertEqual(tree1, test_tree1)
        
        with open("text.maf.idx", "r") as inf:
            for i, line in enumerate(inf):
                self.assertEqual(line, test_line[i])
        
        # check that the two trees are identical
        maftk2 = MafTK()
        tree2 = maftk2.read_index("text.maf.idx")
        for k in tree1:
            self.assertEqual(tree1[k], tree2[k])
        
        # cleaning up
        if os.path.isfile("text.maf"):
            os.remove("text.maf")
        if os.path.isfile("text.maf.idx"):
            os.remove("text.maf.idx")
        
    def test_three(self):
        with open("text.maf", "w") as outf:
            outf.write(maf_text2)
            
        # create index
        maf = MafTK()
        tree = maf.make_maf_index(["text.maf"], "text.maf.idx")
                
        # access alignments
        aligns = maf.get_alignments("hg18.chr7",  27578838, 27707221+3)
        print(str(aligns[0][0].seq))
        print(str(aligns[1][0].seq))
        print(str(aligns[2][0].seq))
        self.assertEqual(len(aligns), 3)
        self.assertEqual(len(aligns[2]), 4)
        self.assertEqual(str(aligns[1][0].seq), "TAAAGA")    
        
        ## access alignments
        aligns = maf.get_alignments("mm4.chr6",  53215344+5, 53310102+10)
        self.assertEqual(len(aligns), 3)
        self.assertEqual(len(aligns[1]), 5)
        self.assertEqual(str(aligns[0][4].seq), "GGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG") 
        
        #cleaning up
        if os.path.isfile("text.maf"):
            os.remove("text.maf")            
        if os.path.isfile("text.maf.idx"):
            os.remove("text.maf.idx")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
    