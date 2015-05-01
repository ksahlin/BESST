BESST v1.3
======

OBS:
----
#### Common pitfall: ####

If --orientation is not specified, BESST assumes that all libraries was aligned in fr orientation.
In versions less than 1.2 BESST cannot parse rf orientations. Thus, BESST requires reads to be mapped in FR mode, i.e. --->  <---, matepairs thus need to be reverse complemented.


#### Time and memory requirements: ####
Version 1.3 and later have implemented several major improvements in runtime and memory requirements. These fixes has the most effect on large fragmented assemblies with hundereds of thousands to millions of contigs.

#### Bam files and mapping  ####
BESST requires sorted and indexed bamfiles as input. Any read aligner + samtools can be used to obtain such files. Read pairs needs to be aligned in paired read mode. BESST provides a script (https://github.com/ksahlin/BESST/blob/master/scripts/reads_to_ctg_map.py) for obtaining sorted and indexed bamfiles with BWA-sampe or BWA-mem in one go. An example call for mapping with this script is

```sh
python reads_to_ctg_map.py /path/to/lib1_A.fq /path/to/lib1_A.fq /path/to/contigs.fasta --threads N
```
where N is an integer specifying how many threads BWA-mem should use. --nomem can be specified to the above call to use BWA-sampe as the paired read alignment pipeline.
 
INPUT:
------
Required arguments:

* -c < path to a contig file >  

* -f < path to bamfiles >  (increasing order of insert size)

* -o < path to location for the output >

Highly reccomended argument:

* --orientation < fr/rf one for each library >

EXAMPLE RUN:
-----------
For scaffolding with one PE and one MP library:
```sh
runBESST -c /path/to/contigfile.fa -f /path/to/file1.bam /path/to/file2.bam -o /path/to/output --orientation fr rf
```
If the mate pair library was reversed complemented before it was aligned, '--orientation fr fr' should be specified.

Optional arguments:
-------------------

The following arguments are computed internally / set by BESST. It is however good to specify mean and standard deviation if your assembly is very fragmented compared to the library insert size (not enough large contains to compute library statistics on).

#### Library parameters: ####

* -m < the means of the insert sizes of the library, one integer number for each library > (integer numbers)

* -s < standard deviation of the libraries, one integer number for each library> (integer numbers)
 
* -T < Thresholds that are some upper levels that you think not too many PE/MP will have a longer insert size than (in the end of the mode) > (integer numbers) 

* -r < Mean read length for each of the libraries > (integer number) 

* -z <ints> Coverage cutoff for repeat classification ( e.g. -z 100 says that contigs with coverages over 100 will be discarded from scaffolding). Integer numbers, one for each library) 

* -d Check for sequencing duplicates and count only one of them (when computing nr of links) if they are occurring.

#### Contig parameters: ####

* -k <ints> Minimum contig size to be seen as a "large contigs" (for statistical scoring only). One number for each library.

* --filter_contigs <int> Remove contigs smaller than this value from all scaffolding. These contigs are not even incuded in the outpus of BESST.


#### Algorithm parameters: ####

* --iter <int> Maximum number of iterations in BFS search for paths between scaffolds.

* --max_extensions <int> Maximum number of extensions performed between contig/scaffold ends.

* -e <int> The least amount of witness links that is needed to create a link edge in graph (one for each library)

* -y Extend scaffolds with smaller contigs (default on).

* --no_score Do not perform statistical scoring, only run path search between contigs.

* --score_cutoff <float> Only consider paths with score > score_cutoff (default 1.5)



#### Under construction / proven unstable ####

* -g < Haplotype detection function on or off > ( default = 0 (off) <0 or 1> )

* -a < Maximum length difference ratio for merging of haplotypic regions> (float nr)
 
* -b < Nr of standard deviations over mean/2 of coverage to allow for clasification of haplotype> (integer value) 

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


