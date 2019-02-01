# PyBlast

PyBlast is a small library using mostly BioPython, designed to process easily local blast requests and global alignment on multiple cores. 

## Local alignment

A blast command line can be launched using `BCLine`. By default only one core is used, unless `ncore` is provided (see below for more info about how multiprocessing is done).

```python3
from pyblast import BCLine

query = "path/to/your/query"
db = "path/to/your/blastdb"

bcl = BCLine("blastp", query=query, subject=db, word_size=3)
print (bcl.run(ncore=4, chunksize=20))
```

Most of the process is similar to usual biopython [blastcommandline functions](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc98). Additionally, PyBlast includes `BCLine6`, which automatically process blast table output (`outfmt 6`) and returns a pandas dataframe.

```python3
from pyblast import BCLine6

query = "path/to/your/query"
db = "path/to/your/blastdb"

bcl = BCLine6("blastp", query=query, subject=db)
print (bcl.run(ncore=4, chunksize=20))
```

The dataframe includes several columns, such as identity, coverage and percentage of each of these columns compared to the query and subject size.

## Global alignment

Global alignment is provided by the BioPython [`pairwise2`](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc86) module. By default, the `pairwise2.align.globalxx` function is used. Before using this function, I highly encourage to read the biopython documentation about this module since most of the results while depend of the input function and parameters. Within PyBlast, local alignments can be launched using the `GlobalAlignment` object which took as parameter an iterable of `SeqRecord` pairs. In this example, we will just make local alignments of every pairs of a single fasta file.

```python3
from Bio import SeqIO
from itertools import combinations
from pyblast import GlobalAlignment

records = SeqIO.parse(query, "fasta")
pairs = combinations(records, 2)

gal = GlobalAlignment(pairs)
print (gal.run(ncore=4))
```

By default, results are only output alignments. Tou can also use the callback function `run_pair_identity` which will automatically determine query and subject similarity of each alignments.

```python3
records = SeqIO.parse(query, "fasta")
pairs = combinations(records, 2)

gal = GlobalAlignment(pairs)
callback = GlobalAlignment.run_pair_identity
print (gal.run(ncore=4, callback=callback))
```

Note that the biopython alignment function and paramaters can be changed using the `gfun`, `args` and `kwargs` attributes.

```python3
args = (5, -4, -2, -0.5)
gal = GlobalAlignment(pairs)
callback = GlobalAlignment.run_pair_identity
print (gal.run(* args, ncore=4, callback=callback, gfun="localms", one_alignment_only=True))
```

## Combining local and global alignments

PyBlast provides a way to combine (blast) local and (biopython) global alignment into a single object : `Global2Local`. The idea is to select only the best local alignments for each query and produce a global alignments for the selected pairs. By default, the `BCLine6.get_best` function is used, which took the 10 best local alignments based on the query positive value (`qpos`). Both local and global functions can be modified.

```python3
from pyblast import Local2Global

l2g = Local2Global(query=query, subject=db)
res = l2g.run("blastp", ncore=4, chunksize=20, gkwargs={"one_alignment_only" : True, "gfun" : None})
print (res)
```

## Some notes about multiprocessing

For multicore (`ncore` > 1) Blast analysis, the query fasta file is splitted into several sub fasta files with a predifined sequences number called `chunksize` and the results are afterward merged together. This option can be combined with the blast option `num_threads`. However, while an increasing number of threads allows a lower memory consomption and increase execution speed, splitting the fasta file and running multiple blast tasks is way faster. 

Below is shown the time used to blast 3000 random *S. cerevisiae* protein sequences against the *S. cerevisiae* proteom ([SGD](https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/)), full data (50x50) can be found in the `blast_time` directory. Note that while `num_thread` does not speed up that much execution time, increasing only the number of ncore can lead to memory overflow.

![blast_performance](/Users/jsgounot/Documents/Script/SFTP_ST/PyLibs/PyPan/blast_time/time.50.png)
