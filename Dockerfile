FROM dhspence/docker-genomic-analysis
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL \
  description="Haloplex QC scripts"

COPY CalculateCoverageQC.071417.sh /usr/local/bin/
COPY CoveragePlots.R /usr/local/bin/
