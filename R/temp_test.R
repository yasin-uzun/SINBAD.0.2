
base_counts = foreach(i=1:5, .export= ls(globalenv()) ) %dopar%
{
  library(Rsamtools)

  bam_file = bam_files[i]
  print(paste(i, 'Computing coverage rate for', bam_file))
  res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
  print(paste(i, 'Computed coverage rate for', bam_file))
  base_count = nrow(res)
  print(base_count)
  base_count
}
