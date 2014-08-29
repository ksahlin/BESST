BESST
======

OBS:
----
Common pitfall:
BESST requires reads to be mapped in FR mode, i.e. --->  <---, matepairs thus need to be reverse complemented.

 
INPUT:
------
Required arguments:

* -c < path to a contig file >  

* -f < path to bamfiles >  (increasing order of insert size)

* -o < path to location for the output >

EXAMPLE RUN:
-----------
```sh
runBESST -c /path/to/contigfile.fa -f /path/to/file1.bam /path/to/file2.bam -o /path/to/output 
```
Optional arguments:
-------------------

The following arguments are computed internally by BESST. It is however good to specify mean and standard deviation if your assembly is very fragmented compared to the library insert size (not enough large contains to compute library statistics on).

* -m < the means of the insert sizes of the library, one integer number for each library > (integer numbers)

* -s < standard deviation of the libraries, one integer number for each library> (integer numbers)
 
* -T < Thresholds that are some upper levels that you think not too many PE/MP will have a longer insert size than (in the end of the mode) > (integer numbers) 

* -r < Mean read length for each of the libraries > (integer number) 

* -e < The least amount of witness links that is needed to create a link edge in graph (one for each library) > (integer number) 

* -k < Minimum contig size to include in the scaffolding in each scaffolding step >  (One number for each library (default 0) (integer numbers))

* -z < Coverage cutoff for repeat classification > ( e.g. -z 100 says that contigs with coverages over 100 will be discarded from scaffolding. ) (integer numbers, one for each library) 

* -g < Haplotype detection function on or off > ( default = 0 (off) <0 or 1> )

* -a < Maximum length difference ratio for merging of haplotypic regions> (float nr)
 
* -b < Nr of standard deviations over mean/2 of coverage to allow for clasification of haplotype> (integer value) 
 
* -d < check for sequencing duplicates and count only one of them (when computing nr of links) if they are occurring> (default on = 1). <0 or 1>  

* -y < 0 or 1 > Extend scaffolds with smaller contains (default on).

* -q < optinal flag > Parallellize work load of path finder module in case of multiple processors available using multiprocessing library for pyhton.


NOTE:
-------

1. Definition of insert size: BESST assumes the following definitions: 
  * insert size = fragment length. 
  * For PE: 
```
   s                    t
   ------>      <-------
```
from s to t, that is, insertsize = readlen1 + gap + readlen2. So when mean and/or threshold is supplied, it should be mean and threshold of this distance.

2. Mapping reads: If you want to map a mate pair library, you will need to map them as paired end library i.e. forward-reverse mode. All read pair libraries should be in this order.

3. Order of scaffolding: It is crucial for the algorithm that you give the libraries in increasing order of the insert size.


