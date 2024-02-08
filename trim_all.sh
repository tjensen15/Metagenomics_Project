#!/bin/bash

ADAPTER_FWD="ACTCCTACGGGAGGCAGCAG"
ADAPTER_REV="GGACTACHVGGGTWTCTAAT"

for i in $(seq -w 28 51); do
    cutadapt -a $ADAPTER_FWD -A $ADAPTER_REV -o trimmed_SRR229935${i}_1.fastq -p trimmed_SRR229935${i}_2.fastq SRR229935${i}_1.fastq SRR229935${i}_2.fastq
done
