
# Material and methods

"Enhancer distance matters: forced linear recruitment of a 
super-enhancer strongly reactivates the developmentally silenced 
fetal HBG globin genes"

 
## ATAC-seq 

ATAC-seq was performed using the Omni-ATAC protocol. Briefly, 200,000 cells were lysed with 0.1% NP-40, 0.1% Tween-20 and 0.01% digitonin, and incubated with homemade Tagment DNA Enzyme for 30 min at 37 °C. DNA was purified with QIAGEN MinElute Reaction Cleanup Kit. Library fragments were amplified using Phusion High-Fidelity PCR Master Mix with HF Buffer (Thermo Fisher Scientific, catalog no. F531S) and customized primers with unique single or dual indexes. Libraries were purified using AMPure XP beads (Beckman Coulter, catalog no. A63881) at a 1:1 ratio and according to the manufacturer’s instructions. Constructed libraries were evaluated using the Agilent Bioanalyzer 2100 with the DNA 7500 kit (catalog no. 5067-1504). 

Libraries were then pooled and sequenced in single (SE) or paired-end (PE) mode using a P2 flow cell on the NextSeq 2000 using Illumina-supplied kits as appropriate. 

Reads were trimmed with cutadapt 4.4 with --minimum-length 10 and removal of adapters, 
aligned to the GRCh38.p14 reference genome with bowtie2 with parameters 
--local --very-sensitive-local --no-unal --no-mixed --no-discordant --dovetail -X 1000.

Non-primary alignments, not-properly paired reads (in case of PE reads) and alignments with mapping quality 
lower than 15 were filtered out with
samtools 1.17 (view - F 260) and blacklisted against ENCODE Unified GRCh38 Blacklist ENCFF356LFX
(https://www.encodeproject.org/files/ENCFF356LFX/) with bedtools 2.31.0.

For visualization bigwigs were created with a method similar to the one in the ENCODE atac-seq-pipeline 
(https://github.com/ENCODE-DCC/atac-seq-pipeline). 
First we converted the bams to tagAlign BED files with bamtobed (bedtools 2.31.1) and 
created treatment bedgraphs during macs2 (2.2.9.1) peak calling, 
with artificial reads centralized 75bp up- and downstream of Tn5 transposase cutting sites,
with parameters -f BED -g hs -p 0.01 --shift -75 --extsize 150 --nomodel --bdg --SPMR --keep-dup all.
Treatment bedgraphs were converted to bigwigs with bedGraphToBigWig (ucsc tools 3.3.5).

Bigwigs were displayed using rtracklayer 1.62.2 and GenomicRanges 1.54.1 packages.
Bigwigs for clones being biological replicates were averaged with bigwigAverage (deeptools=3.5.2) 
and the final signal track was normalized against the average signal around transcription start sites (500bp up- and downstream).

Snakemake pipeline will be publicly available on GITHUB later. 

  - cutadapt=4.4
  - bedtools=2.31.0
  - bowtie2=2.5.1
  - macs2=2.2.9.1
  - samtools=1.17


## CUT&RUN 

CUT&RUN was performed using the CUT&RUN Kit from Cell Signaling (catalog no. 86652S) according to the manufacturer’s instructions. Briefly, 0.5 million cells were harvested and washed twice in 1× wash buffer at room temperature and bound to 15 μl of Concanavalin A beads. Primary antibody incubation was performed in PCR tubes for 2 h at 4 °C with rotation. H3K4me3 antibody supplied by the kit was used (1:xyz). After two washes, pAG-MNase was added and incubated for 1 h at 4 °C with rotation. After three washes, cells were resuspended in 1× wash buffer and pAG-MNase digestion was activated by adding 2 mM CaCl2. Samples were incubated on ice for 45 min, then 2× stop buffer was added to end the reaction. The samples were incubated at 37 °C for 10 min to release the captured chromatin fragments. The sequencing libraries were prepared using the NEBNext Ultra II DNA Library (New England Biolabs, catalog no. E7645) with 15–30 ng of input DNA. 

Sequencing libraries were pooled and pair-end (2 × 50 bp) sequenced on an Illumina Nextseq 2000 platform.

Reads were trimmed with cutadapt 4.4 with --minimum-length 10 and removal of adapters, 
aligned to the GRCh38.p14 reference genome with bowtie2 with parameters 
--local --very-sensitive-local --no-unal --no-mixed --no-discordant --dovetail -X 1000.

Non-primary alignments, not-properly paired reads (in case of PE reads) and 
alignments with mapping quality lower than 1 were filtered out with
samtools 1.17 (view - F 260) and blacklisted against ENCODE Unified GRCh38 Blacklist ENCFF356LFX
(https://www.encodeproject.org/files/ENCFF356LFX/) with bedtools 2.31.0. 

Duplicates where removed with picard MarkDuplicates --REMOVE_DUPLICATES true (picard 2.27.5)  

Bigwigws were created using bamCoverage (deeptools 3.5.5) with parameters
bamCoverage --binSize 10 --normalizeUsing RPGC --ignoreForNormalization chrX chrM --extendReads

Bigwigs were displayed using rtracklayer 1.62.2 and GenomicRanges 1.54.1 packages.
Bigwigs for clones being biological replicates were averaged with bigwigAverage (deeptools=3.5.2) 
and the final signal track was normalized against the signal at the HBA1 and HBA2 loci.

  - cutadapt=4.4
  - bedtools=2.31.0
  - bowtie2=2.5.1
  - deeptools=3.5.2
  - picard=2.27.5
  - samtools=1.17


Snakemake pipeline may be publicly available on GITHUB later.  
