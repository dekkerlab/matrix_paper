module load bedtools/2.28.0

# Extract the loop anchor information from the loops

## 2D to 1D 
cat loop_list_A.bedpe | awk '{if (NR!=1) {print}}' | \
awk -F "\t" '{print $1"\t"$2"\t"$3 } ' > loop_list_A_1.bed
cat loop_list_A.bedpe | awk '{if (NR!=1) {print}}' | \
awk -F "\t" '{ print $4"\t"$5"\t"$6 } ' > loop_list_A_2.bed
cat loop_list_A_1.bed loop_list_A_2.bed | sortBed > loop_list_A_both_anchors.bed
# bedtools to merge same anchors into one list
bedtools merge -i loop_list_A_both_anchors.bed -c 1 -o count > loop_list_A_both_anchors.merged.bed 


# Get the number of loops for each condition 
 cat loop_list_A.bedpe | wc -l
# Get the number of anchors for each conditions
 cat loop_list_A_both_anchors.merged.bed | wc -l
 
# To plot the loops vs loop anchors use dots_vs_anchors_scatter.py

# Used anchor_diff.ipynb to plot the # of anchors vs # of loops per anchor


## 2D to 1D for all loop lists

dot_file_list=(cloops_U54-ESC4DN-FA-DpnII-R1-R2_hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	cloops_U54-ESC4DN-DSG-DpnII-R1-R2_hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	cloops_U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	cloops_U54-HFFc6-FA-DpnII-R1-R2_hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	cloops_U54-HFFc6-DSG-DpnII-R1-R2_hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	cloops_U54-HFFc6-FA-DSG-MNase-R1-R3.hg38.mapq_30.1000.mcool.combined.bedpe.postproc
	)

for file in ${dot_file_list[@]}; do
	filename_temp=$(echo $file | cut -f11 -d"/");
	filename=$(echo $filename_temp | cut -f1 -d".");
	#echo $filename;done
	cat $file | awk '{if (NR!=1) {print}}' | awk -F "\t" '{print $1"\t"$2"\t"$3 } ' > $filename'_1.bed';
	cat $file | awk '{if (NR!=1) {print}}' | awk -F "\t" '{ print $4"\t"$5"\t"$6 } ' > $filename'_2.bed';
	cat  $filename'_1.bed' $filename'_2.bed' | sortBed > $filename'.bed';
	cat $filename'.bed' | wc -l;
	bedtools merge -i $filename'.bed' -c 1 -o count > $filename'.merged.bed';
	cat $filename'.merged.bed' | wc -l;done



# Find the anchor overlaps between 3 conditions

# Steps
# 1. find pairwise overlaps between 3 conditions, A, B, C
# output files= A/B, B/A, AnB,AuB,  A/C, C/A, AnC, AuC,  B/C, C/B, BnC, BuC
# 2. AnB overlap C 
# output files= (AnB)nC, (AnB)/C, C/(AnB)
# 3. AnC overlap B 
# output files= (AnC)nB, (AnC)/B, B/(AnC)
# 4. BnC overlap C 
# output files= (BnC)nA, (BnC)/A, A/(BnC)
# 5. A overlap (BuC)
# output files= A/(BuC), An(BuC), (BuC)/A
# 6. B overlap (AuC)
# output files=B/(AuC), Bn(AuC), (AuC)/B
# 7. C overlap (AuB)
# output files=C/(AuB), Cn(AuB), (AuB)/C

# final output for overlapping 3 conditions
# (AnB)/C, (AnC)/B, (BnC)/A, (AnB)nC, A/(BuC), B/(AuC), C/(AuB)


# to find intersect 
bedtools intersect -wa -wb -a $file1 -b $file2 > $file1 n $file2

# to find the difference
bedtools intersect -v -a $file1 -b $file2 > $file1 / $file2

# to find the union
cat $file1 $file2 | sortBed > $file1 u $file2

# Get 2 lists from the analysis above 
# list1= intersect of 3 anchors lists
# list2=FA+DSG-MNase only anchor list 

# overlap list1 and list 2 with ATAC Seq peaks of the specificed cell type 
# Then pileup the list1_overlap_ATAC_Seq_peaks in CTCF, SMC1, H3K4me3, H3K27ac, RNA Pol II and YY1 using deeptools (https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html) 


HFFc6_ATACSeq_peaks=HFFc6_Atacseq.mRp.clN_peaks.narrowPeak
H1-hESC_ATACSeq_peaks=hESC_Atacseq.mRp.clN_peaks.narrowPeak

# compute the matrix for deeptools 

computeMatrix reference-point -S input.bw \
-R list1_overlap_ATAC_Seq_peaks list2_overlap_ATAC_Seq_peaks \
--beforeRegionStartLength 1000 \
--referencePoint "center" \
--afterRegionStartLength 1000 \
--skipZeros -o output.matrix.gz

# plot the enrichment profiles using deeptools

plotProfile -m $output.matrix.gz \
--averageType mean \
--legendLocation upper-left \
--regionsLabel Intersection FA+DSG-MNase-Specific \
-out $filename.pdf;done

## Candidate Cis regulatory Elements

# downloaded the Candidate Cis Regulatory Elements from ENCODE (https://screen.encodeproject.org), Select H1-hESC (H1-hESC_cCRE.bed) and HFF (HFF_cCRE.bed). 

cat H1-hESC_cCRE.bed | awk -F "\t" '{if ($10!="Low-DNase") {print $0}}' > H1-hESC_cCRE_subset.bed
cat HFF_cCRE.bed | awk -F "\t" '{if ($10!="Low-DNase") {print $0}}' > HFF_cCRE_subset.bed

# get the loop anchors

cat Anchor_list_A.bedpe | awk '{if (NR!=1) {print}}' | \
awk -F "\t" '{print $1"\t"$2"\t"$3 } ' > Anchor_list_A_1.bed
cat Anchor_list_A.bedpe | awk '{if (NR!=1) {print}}' | \
awk -F "\t" '{ print $4"\t"$5"\t"$6 } ' > Anchor_list_A_2.bed

# example of HFF anchors

HFF_Anchor2=Anchor_list_A_2.bed
HFF_Anchor1=Anchor_list_A_1.bed


hff=HFF_cCRE_subset.bed


bedtools intersect -wa -wb -a $hff -b $HFF_Anchor1 | sort -u -k4,4 > HFF_Anchor1_overlap_cCRE.bed
bedtools intersect -wa -wb -a $hff -b $HFF_Anchor2 | sort -u -k4,4 > HFF_Anchor2_overlap_cCRE.bed


# To extract the number of only Promoter or Enahncer loops

cat HFF_Anchor1_overlap_cCRE.bed | sort -u -k1,1 -k2,2 | awk -F "\t" '{if (($13=="PLS") || ($13=="pELS") || $13=="dELS" ) {print $0}}' | wc -l
cat HFF_Anchor2_overlap_cCRE.bed | sort -u -k1,1 -k2,2 | awk -F "\t" '{if (($13=="PLS") || ($13=="pELS") || $13=="dELS" ) {print $0}}' | wc -l



#### Plot loop anchors seperately 

# Scanned for CTCF motifs using http://meme-suite.org/tools/fimo and downloaded ctcf_fimo_hg38_sorted.bed
# (1) Overlap CTCF motifs (ctcf_fimo_hg38_sorted.bed) with CTCF peaks (HFFc6_CTCF_CT.mRp.clN_peaks.narrowPeak)
# (2) Overlap (1) with Anchor 1 and Anchor 2 seperately
# Then plot the enrichment profile of this (2) in CTCF, SMC1, H3K4me3, H3K27ac, PolII and YY1. 


ctcf_hff=HFFc6_CTCF_CT.mRp.clN_peaks.narrowPeak

motif=ctcf_fimo_hg38_sorted.bed

# overlap wit motif_peak

bedtools intersect -wa -wb -a $motif -b $ctcf_hff > HFF_CTCF_with_motifs.bed



input_path=path2anchorlist

for i in $(ls $input_path | grep ".bed$" | grep "HFF" ) ; do \
filename=$(echo $i | cut -f1 -d"b");
echo $filename;
bedtools intersect -wa -wb -a HFF_CTCF_with_motifs.bed -b $input_path/$i > $filename'motif_peak.bed';done

# Match Anchor 1 with Anchor 2 and select the anchors that have convergent CTCF motifs enriched for CTCF peak

Rscript merge_func.R filename'_Anchor1_motif_peak.bed' $filename'Anchor2_motif_peak.bed' filename'_12.motif_peak.bed'

cat filename'_12.motif_peak.bed'  | awk '{if (NR!=1) {print}}'  | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$27"\t"$28"\t"$29"\t"$30}' | uniq |  awk '{if (($4=="+")&&($8=="-")) {print $0}}' | uniq > $filename'.motif_peak.deeptools_input.bed'
cat $filename'.motif_peak.deeptools_input.bed' | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > $filename'.motif_peak.deeptools_input1.bed'
cat $filename'.motif_peak.deeptools_input.bed' | awk '{print $5"\t"$6"\t"$7"\t"$8}' | uniq > $filename'.motif_peak.deeptools_input2.bed'

# 

computeMatrix reference-point -S HFFc6_H3K4me3.bigWig \
-R $filename'.motif_peak.deeptools_input1.bed' $filename'.motif_peak.deeptools_input2.bed' \
--beforeRegionStartLength 1000 \
--referencePoint "center" \
--afterRegionStartLength 1000 \
--skipZeros -o $filename'.2000.motif_peak.matrix.gz'


for i in $(ls $input_path | grep ".2000.motif_peak.matrix.gz$" ) ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $filename;done
plotProfile -m $i \
--averageType mean \
--legendLocation upper-left \
--regionsLabel Micro-C-Anchor1 Micro-C-Anchor2 \
-out $filename.2000.prof.motif_peak.pdf;done






















