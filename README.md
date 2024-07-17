# mouse_drip-seq_liver_hepa
DRIP-seq analysis in mouse HEPA1-6 and liver tissue, as described in https://www.biorxiv.org/content/10.1101/2024.05.07.592855v1 which is currently in revision.

Puzzo_Crossley_DRIP-seq_2024.sh ;
This is the script used to align raw sequencing reads to the mm10 genome and generate DRIP-seq tracks.
This script also includes includes some further data analysis: peak calling, peak annotation, aggregating signal over peaks and genes, analyzing nucleotide content of peaks and the density of non-B-DNA elements around DRIP peaks.

DRIP peaks from mouse HEPA1-6 cells and liver tissue are included (annotated with strand from gene strand they intersect) ;

hepa_annotated_peaks.bed
liver_annotated_peaks.bed


mm10_gencodev25_basic.6cols.bed ;
Also are included the mm10 genes used for analysis (e.g. RNA-seq, TSS and TES plots).

drip_seq_peaks_analysis_puzzo_crossley_et_al.ipynb ;
Python code used for Figure 1

G4_around_peaks_puzzo_crossley_et_al.ipynb ;
Python code used for Extended Fig 2


