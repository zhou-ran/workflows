#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# @Date    : 2020-01-08 21:15:28
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

function showHelp { 
cat << EOF
Usage: PSI.sh Command
Description:
    Derived from http://www.currentprotocols.com/protocol/hg1116, adding <path_to_bedtools2.23> for specify bedtools2.23

Usage:
    bash ExonicPartPSI.sh <path_to_bedtools2.23> <ExonicPart.gff> <alignment_file.bam> <readLength> <SJ.out.tab> <baseName>

Mandatory:
    <path_to_bedtools2.23>      Path to bedtools, version 2.23
    <ExonicPart.gff>            ExonicPart gff, deriverd from dexseq_prepare_annotation.py
    <alignment_file.bam>        RNAseq aligned bam files
    <readLength>                RNAseq reads length, only one side.
    <SJ.out.tab>                SJ.out.tab file.
    <baseName>                  Prefix of outputfiles.

For details see: "Schafer, S., et al. 2015. Alternative splicing signatures in RNA-seq data: percent spliced in (PSI). Curr. Protoc. Hum. Genet."
EOF
exit 1 
} 

set -e
set -u
set -o pipefail

[ $# -lt 6 ] && showHelp
BEDTOOLS=$1
GFF=$2
INBAM=$3
readLength=$4
JUNCTIONS=$5
PREFIX=$6

echo "Make junctions bed file...."
sort -k1,1 -k2,2n -k3,3n ${JUNCTIONS} | awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300" }' > ${PREFIX}_junctions.bed


echo "Counting exon coverage...."
${BEDTOOLS} coverage -split -abam ${INBAM} -b ${GFF} | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$5-$4+1,$9,$10}' | sort -k5,5 > ${PREFIX}_exonic_parts.inclusion

sed 's/,/\t/g' ${PREFIX}_junctions.bed | grep -v description | awk '{OFS="\t"}{print $1,$2+$13, $3-$14,$4,$5,$6}' > ${PREFIX}_intron.bed

echo "Counting exclusion...."
${BEDTOOLS} intersect -wao -f 1.0 -s -a ${GFF} -b ${PREFIX}_intron.bed | awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | sort -k1,1 > ${PREFIX}_exonic_parts.exclusion
rm ${PREFIX}_intron.bed

echo "Calculating PSI value..."
cut -f5 ${PREFIX}_exonic_parts.inclusion > ${PREFIX}_exonID1.txt 
cut -f1 ${PREFIX}_exonic_parts.exclusion > ${PREFIX}_exonID2.txt 
diff ${PREFIX}_exonID1.txt ${PREFIX}_exonID2.txt > /dev/null || (echo "Unsorted exonID exit" &&  return 9)
rm ${PREFIX}_exonID1.txt ${PREFIX}_exonID2.txt

paste ${PREFIX}_exonic_parts.inclusion ${PREFIX}_exonic_parts.exclusion | awk -v "len=$readLength" 'BEGIN{OFS="\t"; print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' > ${PREFIX}_exonic_parts.psi
rm ${PREFIX}_exonic_parts.inclusion ${PREFIX}_exonic_parts.exclusion

echo "ALL done!"