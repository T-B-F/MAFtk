# MAFtk
Small toolkit to handle and query Multiple Sequence Alignment in MAF format

## python example

Create a test file in MAF.

```python
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
```

Use the MAFtk toolkit functions to parse and create an store positions in an interval tree

```python
from mafio_index inport make_maf_index, create_index
        
# create index
tree1 = make_maf_index(["test.maf"], "test.maf.idx")
    
# read index
tree2 = create_index("test.maf.idx")
    
# check that the two trees are identical
print(tree1 == tree2)
print()
```

Use the MAFtk toolkit functions to access specific location of an alignment

```python    
from mafio_index inport get_alignment

# access alignments
aligns = get_alignment(tree1, "hg18.chr7",  27578838, 27707221+3)
for align in aligns:
    print(align)
print()    
    
# access alignments
aligns = get_alignment(tree1, "mm4.chr6",  53215344+5, 53310102+10)
for align in aligns:
    print(align)
```
