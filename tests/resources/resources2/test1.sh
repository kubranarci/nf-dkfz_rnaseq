#!/bin/bash
#
# Copyright (c) 2020 DKFZ.
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/nf-arriba/blob/master/LICENSE).
#

set -ue
set -o pipefail

outDir="${1:?No outDir set}"
environmentDir="${2:-"$outDir/test-environment"}"
workflowDir=$(readlink -f $(dirname "$BASH_SOURCE")"/..")

TEST_TOTAL=0
TEST_ERRORS=0

runTest() {
	# delete output files from previous run
	rm -f "$outDir/"{fusions.tsv,fusions_discarded.tsv.gz,fusions.pdf,fusions_alignments.bam,fusions_alignments.bam.bai,virus_expression.tsv,virus_alignments.bam,virus_alignments.bam.bai}

	# run test
	(set -x; nextflow run "$workflowDir/main.nf" \
		-profile dkfzLsf,singularity \
		-ansi-log \
		--sophiaSVs="$workflowDir/test/sophia_svs.tsv" \
		--outputDir="$outDir" \
		"$@")

	# check if expected fusion is found and if it is annotated properly
	local ERROR=$(awk -F '\t' '
		NR==1 { for (i=1; i<=NF; i++) col[$i] = i } # get column names
		NR==2 && ($col["split_reads1"] != 2 || $col["split_reads2"] != 2) { print "ERROR: expected fusion not found"; exit(0) }
		NR==2 && $col["closest_genomic_breakpoint1"] == "." { print "ERROR: Sophia SVs are not used properly" }
		NR==2 && $col["tags"] != "Mitelman" { print "ERROR: known fusion not annotated" }
		NR==2 && $col["protein_domains"] !~ /kinase/ { print "ERROR: protein domains not annotated" }
	' "$outDir/fusions.tsv")

	# check if there is expected text in the PDF file
	pdftotext "$outDir/fusions.pdf" - | grep -q "SUPPORTING READ COUNT" || ERROR="$ERROR\nERROR: invalid PDF file"

	# check if the expected viral genome is listed in the virus expression file
	grep -q NC_001357 "$outDir/virus_expression.tsv" || ERROR="$ERROR\nERROR: missing virus expression"

	# check if interesting reads were extracted
	[ $(samtools view "$outDir/fusions_alignments.bam" | wc -l) -ne 0 ] || ERROR="$ERROR\nERROR: missing fusion alignments"
	[ $(samtools view "$outDir/virus_alignments.bam" | wc -l) -ne 0 ] || ERROR="$ERROR\nERROR: missing viral alignments"

	# summarize errors
	if [ -n "$ERROR" ]; then
		echo -e "$ERROR" >>/dev/stderr
		let TEST_ERRORS=($TEST_ERRORS + 1)
	fi
	let TEST_TOTAL=($TEST_TOTAL + 1)
}

testFinished() {
	echo "" >>/dev/stderr
	echo "$TEST_ERRORS of $TEST_TOTAL tests failed." >>/dev/stderr
	if [[ $TEST_ERRORS > 0 ]]; then
		exit 1
	else
		exit 0
	fi
}

# Setup the test environment
mkdir -p "$outDir"
if [[ ! -d "$environmentDir" ]]; then
 	conda env create -f "$workflowDir/test-environment.yml" -p "$environmentDir"
fi
set +ue
source activate "$environmentDir"
set -ue

# run all tests
runTest --fastq1="$workflowDir/test/all_reads_paired-end_1.fastq.gz" --fastq2="$workflowDir/test/all_reads_paired-end_2.fastq.gz"
runTest --fastq1="$workflowDir/test/half1_reads_paired-end_1.fastq.gz,$workflowDir/test/half2_reads_paired-end_1.fastq.gz" --fastq2="$workflowDir/test/half1_reads_paired-end_2.fastq.gz,$workflowDir/test/half2_reads_paired-end_2.fastq.gz"
runTest --fastq1="$workflowDir/test/all_reads_single-end.fastq.gz"
runTest --fastq1="$workflowDir/test/half1_reads_single-end.fastq.gz,$workflowDir/test/half2_reads_single-end.fastq.gz"
runTest --bam="$workflowDir/test/all_reads_paired-end.bam"
runTest --bam="$workflowDir/test/all_reads_single-end.bam"
testFinished

