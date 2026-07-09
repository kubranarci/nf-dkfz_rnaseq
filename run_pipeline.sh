#!/bin/bash
nextflow run main.nf -profile singularity -c conf/test_cluster_arriba.config --outdir results2 -resume