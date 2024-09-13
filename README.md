
# ATAC and CR coverage 

ATAC and C&R coverage plots for 
"Enhancer distance matters: forced linear recruitment of a 
super-enhancer strongly reactivates the developmentally silenced 
fetal HBG globin genes"


## Figure ATAC and Cut&Run coverage

See also [README methods](doc/methods.md)


### ATAC coverage coverage plots

- Normalisation of bigwigs: --SPMR by MACS2 ("MACS will SAVE signal per million reads for fragment pileup profiles") 
- Scaling ATAC on the signal on all TSS sites (500bp up- and downstream)
- Figure specifications [ATAC_coverage_figure.tsv](data/ATAC_coverage_figure.tsv) and
[ATAC_coverage_figure_sup.tsv](data/ATAC_coverage_figure_sup.tsv)
- bigwigs: ```~/Documents/projects/HbF_react/bigwigs atac/output_HbF_react_no_dedup```
  - laat4 ~/source/ATAC-seq/output_HbF_react_no_dedup/bigwig
  - laat2 /storage/shared/HbF_react/ATAC-seq/output_HbF_react_no_dedup/bigwig


### Cut & Run coverage plots

- Normalisation of bigwigs: deeptools bamCoverage --normalizeUsing RPGC reads per
                        genomic content (1x normalization)  
- Alignment with mapping quality 1.
- Scaling on the signal at HBA1 and HBA2
- Figure specifications [CR_coverage_figure.tsv](data/CR_coverage_figure.tsv) and
[CR_coverage_figure_sup.tsv](data/CR_coverage_figure_sup.tsv)

- bigwigs: ```~/Documents/projects/VER9574/CR/low_q_mapping/```
    - laat2 /storage/shared/HbF_react/Cut-and-Run/output_low_q_mapping_new/bigwig

