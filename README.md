# mouse_drip-seq_liver_hepa
DRIP-seq analysis in mouse HEPA1-6 and liver tissue

Puzzo_Crossley_DRIP-seq_2024.sh 
This is the script used to align raw sequencing reads to the mm10 genome and generate DRIP-seq tracks.
This script also includes includes some further data analysis: peak calling, peak annotation, aggregating signal over peaks and genes, analyzing nucleotide content of peaks and the density of non-B-DNA elements around DRIP peaks.

hepa_annotated_peaks.bed
liver_annotated_peaks.bed

DRIP peaks from mouse HEPA1-6 cells and liver tissue are included (annotated with strand from gene strand they intersect)

mm10_gencodevm23.6cols.bed
Also are included the mm10 genes used for analysis (e.g. RNA-seq, TSS and TES plots).
