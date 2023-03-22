#!/bin/bash
#1: file prefix

#if file has header then error can happen
cat $1.bed | sed 's/^chr//' > $1.chr.removed.bed #if bed is GRCH37 type then error can happen
cat $1.chr.removed.bed | awk '{if ($1 ~ /^[0-9]/) print $0}' > $1.chr.removed.bed.num
cat $1.chr.removed.bed | awk '{if ($1 !~ /^[0-9]/) print $0}' > $1.chr.removed.bed.sex
sort -k1n,1 -k2n,2 -k3n,3 $1.chr.removed.bed.num > $1.chr.removed.bed.num.sorted
sort -k1,1 -k2n,2 -k3n,3 $1.chr.removed.bed.sex > $1.chr.removed.bed.sex.sorted #-k6,6
#wc -l "$1"*
cat $1.chr.removed.bed.num.sorted $1.chr.removed.bed.sex.sorted > $1.bed
rm $1.chr.removed.bed $1.chr.removed.bed.num $1.chr.removed.bed.sex $1.chr.removed.bed.num.sorted $1.chr.removed.bed.sex.sorted
