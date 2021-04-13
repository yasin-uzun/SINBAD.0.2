library(doSNOW)
library(vioplot)
library(readr)
library(Rsamtools)
library(ShortRead)


align_sample <- function(read_dir,
                         genomic_sequence_path,
                         alignment_dir,
                         aligner,
                         num_cores,
                         mapq_threshold,
                         main_log_dir,
                         pattern = '')
{
  log_dir= paste0(main_log_dir, '/alignment/')
  dir.create(log_dir, recursive = T)
  setwd(alignment_dir)

  fastq_files_1 = list.files(read_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(read_dir, pattern = "*.fastq")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[grepl(pattern = pattern, fastq_files)]

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')

    clusterExport(cl, ls(.GlobalEnv))
    registerDoSNOW(cl)
    clusterExport(cl, ls(.GlobalEnv))

    foreach(i=1:length(fastq_files)) %dopar%
    {
      fastq_file = fastq_files[i]
      cat('** fastq_file: ', fastq_file, '\n')
      align_cell(read_dir = read_dir,
                 fastq_file = fastq_file,
                 aligner = aligner,
                 genomic_sequence_path = genomic_sequence_path,
                 alignment_dir = alignment_dir,
                 log_dir = log_dir)

    }#foreach
  }else#if(num_cores > 1)
  {
    for(i in 1:length(fastq_files))
    {
      fastq_file = fastq_files[i]
      align_cell(read_dir = read_dir,
                 fastq_file = fastq_file,
                 aligner = aligner,
                 genomic_sequence_path = genomic_sequence_path,
                 alignment_dir = alignment_dir,
                log_dir)

    }#foreach

  }#else


}#run_alignment




align_cell <- function(read_dir, fastq_file, aligner, genomic_sequence_path,
                       alignment_dir, log_dir){

  print('align_cell')
  cat(' -- fastq_file: ', fastq_file, '\n')

  cell_id = gsub('.fastq.gz', '', fastq_file)

  #Align
  if(aligner == 'bismark')
  {

    run_bismark_aligner(read_dir = read_dir,
                        fastq_file_left = fastq_file,
                        genomic_sequence_path = genomic_sequence_path,
                        alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)

  }else if(aligner == 'bsmap')
  {

    run_bsmap_aligner(read_dir = read_dir,
                        fastq_file = fastq_file,
                        genomic_sequence_path = genomic_sequence_path,
                        alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)
  }else if(aligner == 'bs_seeker')
  {

    run_bs_seeker_aligner(read_dir = read_dir,
                      fastq_file = fastq_file,
                      genomic_sequence_path = genomic_sequence_path,
                      alignment_dir = alignment_dir,
                      cell_id = cell_id, log_dir = log_dir)
  }else
  {
    msg = paste0('Unrecognized aligner name: ', aligner, ' . Should be one of: bismark, bs_map, bs_seeker')
    stop(msg)
  }

  gc()
  #
  #
  #Mapq filtering
  if(mapq_threshold > 0)
  {
    filter_mapq(alignment_dir, cell_id, mapq_threshold = mapq_threshold, log_dir = log_dir)
  }#if(mapq_threshold > 0)
  #
  use_mapq_filtered_reads = F
  if(mapq_threshold > 0)
  {
    use_mapq_filtered_reads = T
  }

  gc()


  remove_duplicate_reads(alignment_dir = alignment_dir,
                         cell_id = cell_id,
                         use_mapq_filtered_reads = use_mapq_filtered_reads,
                         duplicate_remover = duplicate_remover, log_dir = log_dir)
  gc()


  filter_non_conversion(alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)

  gc()


  split_lambda(alignment_dir = alignment_dir, cell_id = cell_id, log_dir = log_dir)

  gc()




}



filter_cell <- function(read_dir, fastq_file, aligner, genomic_sequence_path,
                       alignment_dir, log_dir){

  use_mapq_filtered_reads = T
  print('Filtering cell')
  cat(' -- fastq_file: ', fastq_file, '\n')

  cell_id = gsub('.fastq.gz', '', fastq_file)


  print(paste('*Duplicate remover: ', duplicate_remover))


  remove_duplicate_reads(alignment_dir = alignment_dir,
                         cell_id = cell_id,
                         use_mapq_filtered_reads = use_mapq_filtered_reads,
                         duplicate_remover = duplicate_remover, log_dir = log_dir)
  gc()


  filter_non_conversion(alignment_dir = alignment_dir,
                        cell_id = cell_id, log_dir = log_dir)

  gc()


  split_lambda(alignment_dir = alignment_dir, cell_id = cell_id, log_dir = log_dir)

  gc()




}


run_bismark_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                                genomic_sequence_path, alignment_dir,
                                cell_id, log_dir)
{

  setwd(alignment_dir)

  log_sub_dir = paste0(log_dir, '/bismark_aligner/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  sys_command = paste0('bismark --bowtie2 ', bismark_aligner_param_settings
                             ,' --fastq  '
                             ,' --basename ', cell_id
                             ,' --bam ', genomic_sequence_path
                             ,'  ', read_dir, fastq_file_left
                             ,' > ', log_file
                       )

  cat('Aligning ', fastq_file_left, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:Bismark aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment



run_bs_seeker_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                                genomic_sequence_path, alignment_dir,
                                cell_id, log_dir, is_paired = F)
{

  #setwd(alignment_dir)
  setwd(bs_seeker_path)

  log_sub_dir = paste0(log_dir, '/align/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  if(is.null(fastq_file_right)) #single ended
  {
    sys_command = paste0('./bs3-align  '
                         ,' -i ', read_dir, fastq_file_left
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bs_seeker_aligner_param_settings
                         ,' -g  ', bs_seeker_genome_dir, bs_seeker_genome_fasta
                         ,' --db=', bs_seeker_genome_dir
                         #,' > ', log_file
    )



  }else #pair ended
  {
    sys_command = paste0('bs3-align  '
                         ,' -1', read_dir, fastq_file_left
                         ,' -2', read_dir, fastq_file_right
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bs_seeker_aligner_param_settings
                         ,' -g  ', bs_seeker_genome_dir, bs_seeker_genome_fasta
                         ,' --db=', bs_seeker_genome_dir
                         #,' > ', log_file
    )


  }




  cat('Aligning ', fastq_file, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:BS_Seeker aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment

run_bsmap_aligner <- function(read_dir, fastq_file_left, fastq_file_right = NULL,
                                  genomic_sequence_path, alignment_dir,
                                  cell_id, log_dir, is_paired = F)
{

  #setwd(alignment_dir)

  log_sub_dir = paste0(log_dir, '/align/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, fastq_file_left, '.log')

  if(is.null(fastq_file_right)) #single ended
  {
    sys_command = paste0(bsmap_path, '/bsmap  '
                         ,' -a ', read_dir, fastq_file_left
                         ,' -d  ', reference_genome_fasta
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bsmap_aligner_param_settings

                         #,' > ', log_file
                        )





  }else #pair ended
  {
    sys_command = paste0(bsmap_path, '/bsmap  '
                         ,' -a ', read_dir, fastq_file_left
                         ,' -b ', read_dir, fastq_file_right
                         ,' -d  ', reference_genome_fasta
                         ,' -o ', alignment_dir, cell_id, '.bam'
                         ,' ', bsmap_aligner_param_settings

                         #,' > ', log_file
              )



  }




  cat('Aligning ', fastq_file, '\n')
  print(sys_command)
  command_result = system(sys_command)

  sink(log_file, append = T)
  if(command_result == 0)
  {
    cat('Alignment is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR:Bismark aligner failed to run for cell ', cell_id ,'. Exiting the pipeline. Please see the output log.')
  }
  sink()

  return(command_result)
}#run_alignment



filter_mapq <- function(alignment_dir, cell_id, mapq_threshold = 10, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  log_sub_dir = paste0(log_dir, '/filter_mapq/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)

  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  if(mapq_threshold > 0)
  {
    cat('Filtering alignments mapq > ', mapq_threshold, ' :' , cell_id, '\n')
    sys_command = paste0('samtools view -bhq ',mapq_threshold,' ',cell_id,'.bam > ',cell_id,'.mapq_filtered.bam')
    command_result = system(sys_command)
  }else
  {
    print(paste0('Mapq threshold is 0. No mapq filtering done for ', cell_id))
  }
  sink()

  if(command_result == 0)
  {
    cat('Read filtering is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR: Read filtering failed to run ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }


  return(command_result)
}





remove_duplicate_reads <- function(alignment_dir, cell_id,
                                   use_mapq_filtered_reads = T,
                                   duplicate_remover = 'samtools'
                                   , log_dir)
{
   bam_file = paste0(cell_id,'.bam')
   if(use_mapq_filtered_reads)
   {
      bam_file = paste0(cell_id,'.mapq_filtered.bam')
   }
   rmdup_file = paste0(cell_id,'.rmdup.bam')

   log_sub_dir = paste0(log_dir, '/rmdup/')
   dir.create(log_sub_dir, recursive = T, showWarnings = F)

   log_file = paste0(log_sub_dir, cell_id, '.log')

   cat('Using ', duplicate_remover, ' to remove duplicates...')

   if(duplicate_remover == 'picard')
   {
     sorted_file = paste0(cell_id,'.sorted.bam')
     metric_file = paste0(cell_id,'.picard_metrics.txt')

     sys_command = paste0('java -jar ',picard_path,' SortSam ',
                          ' INPUT=', bam_file,
                          ' OUTPUT=', sorted_file,
                          ' SORT_ORDER=coordinate ',
                          ' >> ', log_file)

     command_result = system(sys_command)

     if(command_result == 0)
     {
       cat('Sorted alignments for cell ', cell_id, '\n')
     }else
     {
       stop(paste('ERROR: Picard sort failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
     }

     sys_command = paste0('java -jar ',picard_path,' MarkDuplicates ',
                          ' INPUT=', sorted_file,
                          ' OUTPUT=', rmdup_file,
                          ' REMOVE_DUPLICATES=true ',
                          ' METRICS_FILE=',metric_file)

     command_result = system(sys_command)

     if(command_result == 0)
     {
       cat('Removed duplicate reads for cell ', cell_id, '\n')
     }else
     {
       stop(paste('ERROR: Picard remove duplicates failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
     }

   }#if(duplicate_remover == 'picard')

   if(duplicate_remover == 'samtools')
   {
     sys_command = paste0('samtools rmdup -S ', bam_file, '  ', rmdup_file)
     command_result = system(sys_command)

     if(command_result == 0)
     {
       cat('Removed duplicate reads for cell ', cell_id, '\n')
     }else
     {
       stop(paste('ERROR: Samtools remove duplicates failed. Could not sort alignments (bam file) for cell ', cell_id, '\n'))
     }

   }

}#remove_duplicate_reads



filter_non_conversion <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  rmdup_file = paste0(cell_id,'.rmdup.bam')

  log_sub_dir = paste0(log_dir, '/filter_non_conversion/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  cat('Filtering nonconversion reads  :' , cell_id, '\n')
  sys_command = paste0(bismark_path, '/filter_non_conversion  '
                       , ' --samtools_path ', samtools_path
                       , ' ', filter_non_conversion_param_settings
                       , ' ', rmdup_file
                       , ' > ', log_file)
  command_result = system(sys_command)



  if(command_result == 0)
  {
    cat('Non-conversion filtering is successful for cell ', cell_id, '\n')
  }else
  {
    stop('ERROR: Non-conversion filtering failed to run for ', cell_id, '.\n Exiting the pipeline. Please see the output log.')
  }

  return(command_result)
}


split_lambda <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  noncon_file = paste0(cell_id,'.rmdup.nonCG_filtered.bam')
  sorted_file = paste0(cell_id,'.rmdup.nonCG_filtered.sorted.bam')

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/split_lambda/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  cat('Splitting lambda control  :' , cell_id, '\n')
  sys_command = paste0('samtools sort -o ', sorted_file, ' -T tmp_sort_', cell_id, ' ' , noncon_file)
  print(sys_command)
  command_result_0 = system(sys_command)
  sys_command = paste0('samtools index -b ', sorted_file)
  command_result_1 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', sorted_file, ' ' , lambda_chrom_numbers , ' > ', lambda_file)
  command_result_2 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', sorted_file, ' ' , organism_chrom_numbers , ' > ', organism_file)
  command_result_3 = system(sys_command)
  sink()

  sorted_organism_file = gsub('organism.bam', 'organism.sorted.bam', organism_file)
  sys_command = paste0('samtools sort -o ', sorted_organism_file, ' -T tmp_sort_', cell_id, ' ' , organism_file)
  print(sys_command)
  command_result_4 = system(sys_command)

  sys_command = paste0('samtools index -b ', sorted_organism_file)
  print(sys_command)
  command_result_5 = system(sys_command)


  command_result = command_result_1 + command_result_2 + command_result_3

  if(command_result == 0)
  {
    cat('Lambda control split is successful for cell ', cell_id, '\n')
    cat('Output (organism): ', organism_file, '\n')
  }else
  {
    stop('Lambda control splitting failed to run for ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }

  return(command_result)

}



split_lambda_old <- function(alignment_dir, cell_id, log_dir)
{
  setwd(alignment_dir)
  command_result = 0
  noncon_file = paste0(cell_id,'.rmdup.nonCG_filtered.bam')
  sorted_file = paste0(cell_id,'.rmdup.nonCG_filtered.sorted.bam')

  lambda_file = paste0(cell_id,'.lambda.bam')
  organism_file = paste0(cell_id,'.organism.bam')

  log_sub_dir = paste0(log_dir, '/split_lambda/')
  dir.create(log_sub_dir, recursive = T, showWarnings = F)
  log_file = paste0(log_sub_dir, cell_id, '.log')

  sink(log_file)
  cat('Splitting lambda control  :' , cell_id, '\n')
  sys_command = paste0('samtools sort ', sorted_file, ' -T tmp_sort_', cell_id, ' ' , noncon_file)
  command_result_0 = system(sys_command)
  sys_command = paste0('samtools index -b ', noncon_file)
  command_result_1 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', noncon_file, ' ' , lambda_chrom_numbers , ' > ', lambda_file)
  command_result_2 = system(sys_command)
  sys_command = paste0('samtools view -hb  ', noncon_file, ' ' , organism_chrom_numbers , ' > ', organism_file)
  command_result_3 = system(sys_command)
  sink()

  command_result = command_result_1 + command_result_2 + command_result_3

  if(command_result == 0)
  {
    cat('Lambda control split is successful for cell ', cell_id, '\n')
  }else
  {
    stop('Lambda control splitting failed to run for ', cell_id, '. Exiting the pipeline. Please see the output log.')
  }

  return(command_result)

}


process_bismark_alignment_reports <- function(alignment_dir)
{
  setwd(alignment_dir)
  report_files = list.files(alignment_dir, pattern = "*SE_report.txt")
  row_names = c()
  result_list = list()
  for(report_file in report_files)
  {
    #report_file = "GSM3444126_R2_SE_report.txt"
    print(paste('Reading ', report_file))

    cell_id = gsub('_SE_report.txt', '', report_file)

    lines = readLines(report_file)
    informative_lines = lines[grep('\t', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]
    informative_lines = informative_lines[!grepl('C methylated in Unknown context', informative_lines)]

    if( length(informative_lines) > 0 )
    {
      df_report = read.delim(text = informative_lines, sep = '\t',  header = F)
      class(df_report)
      row_names = gsub(':', '', df_report$V1)
      result_list[[cell_id]] = as.character(df_report$V2)
    }#if
  }#for

  row_names[1] = "Total_reads"
  row_names[2] = "Alignments"
  row_names[3] = "Alignment_rate"

  row_names[4] = "Sequences_with_no_alignments"
  row_names[6] = "Discarded_sequences"

  #unlist(lapply(result_list, length))
  table(unlist(lapply(result_list, length)))


  df_result = do.call(cbind, result_list)
  rownames(df_result) = row_names
  t_df_result =  t(df_result)
  final_table = data.frame(Cell_ID = rownames(t_df_result),  t_df_result, stringsAsFactors = F)

  return(final_table)

}#process_bismark_alignment_reports

#df_alignment_stats = final_table

plot_alignment_stats <- function(sample_name, df_alignment_stats)
{

    par(mfrow = c(3,2))


    #Plot alignment rates
    if(is.character( df_alignment_stats$Alignment_rate) )
    {
      df_alignment_stats$Alignment_rate =   parse_number(df_alignment_stats$Alignment_rate)
    }

    alignment_rates = df_alignment_stats$Alignment_rate
    roof = max(alignment_rates) * 1.2
    hist(alignment_rates, col = '#ffffb3', main = 'Alignment rates', xlim =  c(0, 100), breaks = 20,
         xlab ='%Alignment rate', ylab = 'Number of cells', cex.lab = 1.25,  cex.axis = 1.25)

    median_alignment_rate = median(alignment_rates)
    abline(v = median_alignment_rate, col = 'red', lty = 2)
    legend('topright', legend = c('Median'), col = c('red'), lty = 2, bty = "n")

    #Plot methylation rates

    if(is.character(df_alignment_stats$Total.methylated.C.s.in.CHG.context))
    {
      attach(df_alignment_stats)

      Total.methylated.C.s.in.CH.context = parse_number(Total.methylated.C.s.in.CHG.context) + parse_number(Total.methylated.C.s.in.CHH.context)
      Total.unmethylated.C.s.in.CH.context = parse_number(Total.unmethylated.C.s.in.CHG.context) + parse_number(Total.unmethylated.C.s.in.CHH.context)
      met_CH = 100 * Total.methylated.C.s.in.CH.context / (Total.methylated.C.s.in.CH.context + Total.unmethylated.C.s.in.CH.context)
      met_CpG = parse_number(C.methylated.in.CpG.context)
      met_CHG = parse_number(C.methylated.in.CHG.context)
      met_CHH = parse_number(C.methylated.in.CHH.context)
      detach(df_alignment_stats)
    }else
    {
      attach(df_alignment_stats)

      Total.methylated.C.s.in.CH.context = Total.methylated.C.s.in.CHG.context + Total.methylated.C.s.in.CHH.context
      Total.unmethylated.C.s.in.CH.context = Total.unmethylated.C.s.in.CHG.context + Total.unmethylated.C.s.in.CHH.context
      met_CH = 100 * Total.methylated.C.s.in.CH.context / (Total.methylated.C.s.in.CH.context + Total.unmethylated.C.s.in.CH.context)
      met_CpG = C.methylated.in.CpG.context
      met_CHG = C.methylated.in.CHG.context
      met_CHH = C.methylated.in.CHH.context
      detach(df_alignment_stats)
    }
    #x_labels = c('mCpG', 'mCHG', 'mCHH')
    #roof = max(c(met_CpG, met_CHG, met_CHH))
    #pcts = list(met_CpG, met_CHG, met_CHH)

    x_labels = c('mCpG', 'mCH')
    roof = max(c(met_CpG, met_CH), na.rm = T) * 1.1
    pcts = list(met_CpG, met_CH)

    #pcts[[1]] = pcts[[1]][!is.na(pcts[[1]])]
    #pcts[[2]] = pcts[[2]][!is.na(pcts[[2]])]

    vioplot(pcts, ylim = c(0,  roof),
            col = "#fb8072", names = x_labels, cex.main = 1.25, cex.names = 1.25,
            main = 'Methylation rate distributions', cex.lab = 1.25, cex.axis = 1.25,
            ylab = 'Methylation rate (%)'
    )



    #Plot number of reads
    total_reads = as.numeric(df_alignment_stats$Total_reads)
    aligned_reads = as.numeric(df_alignment_stats$Alignments)
    mapq_read_counts = df_alignment_stats$mapq_read_counts
    rmdup_read_counts = df_alignment_stats$rmdup_read_counts
    nc_filtered_read_counts = df_alignment_stats$nc_filtered_read_counts
    organism_read_counts = df_alignment_stats$organism_read_counts

    roof = max(total_reads/1000000)

    x_labels = c('Total', 'Aligned', 'Mapq. filt.', 'Rm-dup.', 'NC filt.', 'Org.')


    vioplot(
    #boxplot(
            total_reads/1000000,
            aligned_reads/1000000,
            mapq_read_counts/1000000,
            rmdup_read_counts/1000000,
            nc_filtered_read_counts/1000000,
            organism_read_counts/1000000,
            ylim = c(0, roof),
            col = "#8dd3c7", names = x_labels,
            ylab = 'million reads',
            main = 'Read statistics', las = 2,
            cex.axis = 1.25, cex.names = 1.25,
            cex.main = 1.25, cex.lab = 1.25
            #outline =F
            )    #Plot number of reads









    #Plot total Cs
    total_cpgs = as.numeric(df_alignment_stats$Total.methylated.C.s.in.CpG.context) +
                 as.numeric(df_alignment_stats$Total.unmethylated.C.s.in.CpG.context)

    total_chgs = as.numeric(df_alignment_stats$Total.methylated.C.s.in.CHG.context) +
      as.numeric(df_alignment_stats$Total.unmethylated.C.s.in.CHG.context)

    total_chhs = as.numeric(df_alignment_stats$Total.methylated.C.s.in.CHH.context) +
      as.numeric(df_alignment_stats$Total.unmethylated.C.s.in.CHH.context)

    total_chs = total_chgs + total_chhs

    #x_labels = c('CpG', 'CHG', 'CHH')
    #vioplot(total_cpgs, total_chgs, total_chhs, #ylim = c(0,  roof),
    #        col = "lightblue", names = x_labels, main = 'Read statistics')

    #hist(total_cpgs/1000000, col = '#bebada',
    #     main = 'CpG sites called',
    #     xlab = 'million CpG sites',
    #     ylab = 'Number of cells', cex.lab = 1.25, , cex.axis = 1.25)




    #Coverage rates
    coverage_rates =  df_alignment_stats$coverage_rate
    coverage_rates = coverage_rates[!is.na(coverage_rates)]
    median_coverage_rate = median(coverage_rates)
    hist(coverage_rates,
         main = 'Genomic coverage by cell',
         xlab = "% Coverage rate",
         breaks = 20, cex.lab = 1.25, cex.axis = 1.25,
         col = 'aquamarine', xlim = c(0, max(coverage_rates) * 1.1)
         )
    abline(v = median_coverage_rate, col = 'red', lty = 2)
    legend('topright', legend = c('Median'), col = c('red'), lty = 2, bty = "n")

    #df_alignment_stats[, c('Cell_ID', 'coverage_rate')]

    ###Scatter plot filtering
    filter_aln_rate =  df_alignment_stats$Alignment_rate  > alignment_rate_threshold
    #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
    filter_read_count = df_alignment_stats$nc_filtered_read_counts > minimum_filtered_read_count

    passing = filter_aln_rate &  filter_read_count

    df_alignment_stats$pass = 0
    df_alignment_stats$pass[passing] = 1

    Alignment.rate = df_alignment_stats$Alignment_rate[passing]
    Filtered.reads = df_alignment_stats$nc_filtered_read_counts[passing] / 1000000
    plot(Alignment.rate,
         Filtered.reads,
         xlim = c(0, max(Alignment.rate)),
         ylim = c(0, max(Filtered.reads)),
         pch = 20, col = 'darkgreen', xlab = '% Alignment', ylab = 'Filtered reads (million)'
         , cex.lab = 1.25, cex.axis = 1.25
         , main = "Cell QC"
         )

    Alignment.rate = df_alignment_stats$Alignment_rate[!passing]
    Filtered.reads = df_alignment_stats$nc_filtered_read_counts[!passing] / 1000000
    points(Alignment.rate,
           Filtered.reads,
           pch = 20, col = 'red')

    abline(v = alignment_rate_threshold, lty = 2)
    abline(h = minimum_filtered_read_count/ 1000000 + 0.01, lty = 2)

    #Filtering pie chart
    cnt_failed_low_aln = sum(!filter_aln_rate & filter_read_count)
    cnt_failed_low_reads = sum(filter_aln_rate & !filter_read_count)
    cnt_failed_low_aln_and_low_reads = sum(!filter_aln_rate & !filter_read_count)

    cnt_passed = sum(passing)
    cnt_failed = sum(!passing)

    #slices = c(cnt_passed, cnt_failed_low_aln, cnt_failed_low_reads, cnt_failed_low_aln_and_low_reads)
    slices = c(cnt_passed, cnt_failed)
    #lbls = c('Passed', 'Low Aln', 'Low Reads', 'Low Aln and Reads')
    lbls = c('Passed', 'Discarded')

    lbls <- paste0(lbls, " (", slices, ")") # add percents to labels
    lbls <- paste(lbls,"",sep="") # ad % to labels

    nice.pie(slices, labels = lbls, col=c('lightgreen', 'pink'), cex = 1.25,
             text_col ='black', main = 'Cell filtering'  )

    #To be added: Genomic Coverage -> Done


    title(paste(sample_name, ': Alignment results'), line = -1, outer = TRUE)

}


count_bam_files <- function(alignment_dir)
{
   setwd(alignment_dir)
   mapq_files = list.files(pattern = '*.mapq_filtered.bam')
   rmdup_files = list.files(pattern = '*.rmdup.bam')
   nc_filtered_files = list.files(pattern = '.*.nonCG_filtered.bam$')
   organism_bam_files = list.files(pattern = '.*.organism.bam$')

   #cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
   #registerDoSNOW(cl)

   #temp = mapq_files[1:10]

   #t1 = Sys.time()
   #temp1 = mcsapply(temp, FUN = countBam, mc.cores = 10)
   #t2 = Sys.time()
   #print(t2-t1)

   library(ShortRead)

   print('Counting mapq filtered bam files')
   df_mapq_read_counts = mcsapply(mapq_files, FUN = countBam, mc.cores = num_cores)
   mapq_read_counts = unlist(df_mapq_read_counts['records', ] )
   names(mapq_read_counts) = gsub('.mapq_filtered.bam', '', names(mapq_read_counts))


   print('Counting rmdup filtered bam files')
   df_rmdup_read_counts = mcsapply(rmdup_files, FUN = countBam, mc.cores = num_cores)
   rmdup_read_counts = unlist(df_rmdup_read_counts['records', ] )
   names(rmdup_read_counts) = gsub('.rmdup.bam', '', names(rmdup_read_counts))


   print('Counting nonconversion filtered bam files')
   df_nc_filtered_read_counts = mcsapply(nc_filtered_files, FUN = countBam, mc.cores = num_cores)
   nc_filtered_read_counts = unlist(df_nc_filtered_read_counts['records', ] )
   names(nc_filtered_read_counts) = gsub('.rmdup.nonCG_filtered.bam', '', names(nc_filtered_read_counts))


   print('Counting organism bam files')
   df_organism_read_counts = mcsapply(organism_bam_files, FUN = countBam, mc.cores = num_cores)
   organism_read_counts = unlist(df_organism_read_counts['records', ] )
   names(organism_read_counts) = gsub('.organism.bam', '', names(organism_read_counts))
   print('Finished counting bam files')

   all_names = Reduce(union,
                     list(
                       names(mapq_read_counts),
                       names(rmdup_read_counts),
                       names(nc_filtered_read_counts),
                       names(organism_read_counts)
                        )
                    )

   df_bam_read_counts <- data.frame(mapq_read_counts = mapq_read_counts[all_names],
                                    rmdup_read_counts = rmdup_read_counts[all_names],
                                    nc_filtered_read_counts = nc_filtered_read_counts[all_names],
                                    organism_read_counts = organism_read_counts[all_names]
                                    )
   head(df_bam_read_counts)

   rownames(df_bam_read_counts) = gsub('.mapq_filtered.bam', '', rownames(df_bam_read_counts))

   #df_bam_read_counts = df_bam_read_counts[!grepl('sorted', rownames(df_bam_read_counts)), ]

   return(df_bam_read_counts)
}#count_bam_files <- function(alignment_dir)


#bam_file = 'Lane1_ACTTGA.rmdup.nonCG_filtered.bam'
compute_coverage_rates <- function(alignment_dir, parallel = T)
{
   setwd(alignment_dir)
   print('***********************')
   print('Computing coverage rates')

   df_chrom_sizes = read.table(chrom_sizes_file, header = F)

   df_chrom_ranges = data.frame(chrom = df_chrom_sizes$V1, start = 1, end = df_chrom_sizes$V2)
   total_genome_size = sum(df_chrom_ranges$end)

   chrom_ranges = makeGRangesFromDataFrame(df_chrom_ranges,
                                      keep.extra.columns=FALSE,
                                      ignore.strand=T,
                                      seqinfo=NULL,
                                      starts.in.df.are.0based=FALSE)

   sbp <- ScanBamParam(which=chrom_ranges )

   p_param <- PileupParam(distinguish_nucleotides=FALSE,distinguish_strands=FALSE,
                        min_base_quality=10, min_nucleotide_depth=1)

   bam_files = list.files(alignment_dir,  pattern = '*organism.sorted.bam$')
   base_counts = c()

   if(parallel)
   {
     cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
     #clusterExport(cl, varlist = ls(), envir = environment())
     #clusterExport(cl, list = ls(), envir = environment())
     #clusterExport(cl, ls(.GlobalEnv))

     registerDoSNOW(cl)
     #clusterExport(cl, varlist = ls(), envir = environment())
     clusterExport(cl, ls(.GlobalEnv))

     #base_counts = foreach(i=1:10, .export= ls(globalenv()) ) %dopar%
     base_counts = foreach(i=1:length(bam_files), .export= ls(globalenv()) ) %dopar%
     {
       library(Rsamtools)

       bam_file = bam_files[i]
       print(paste(i, 'Computing coverage rate for', bam_file))
       #res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)

       #res = NA

       res = tryCatch({
         pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
       }, warning = function(w) {

       }, error = function(e) {
         NA
       }, finally = {

       }
       )

       print(paste(i, 'Computed coverage rate for', bam_file))
       if(is.null(res))
       {
         base_count = NA

       }else
       {
         base_count = nrow(res)

       }

       print(base_count)
       base_count
     }#foreach

  }else{

       base_counts = c()

       for(i in 1:length(bam_files) )
       {
         library(Rsamtools)

         bam_file = bam_files[i]
         print(paste(i, 'Computing coverage rate for', bam_file))
         #res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)

         #res = NA

         res = tryCatch({
           pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
         }, warning = function(w) {

         }, error = function(e) {
           NA
         }, finally = {

         }
         )

         print(paste(i, 'Computed coverage rate for', bam_file))
         if(is.null(res))
         {
           base_count = NA

         }else
         {
           base_count = nrow(res)

         }

         print(base_count)
         base_counts[i] = base_count

       }#for

     }#else





   cell_ids = bam_files
   cell_ids = sub('.organism.sorted.bam', '', cell_ids)
   cell_ids = sub('.organism.bam', '', cell_ids)
   cell_ids = sub('.bam', '', cell_ids)

   names(base_counts) = cell_ids
   base_counts = unlist(base_counts)

   # for(bam_file in bam_files)
   # {
   #   print(bam_file)
   #   res <- pileup(bam_file, scanBamParam=sbp, pileupParam=p_param)
   #   class(res)
   #   dim(res)
   #   head(res)
   #   tail(res)
   #   cell_id = sub('.rmdup.nonCG_filtered.bam', '', bam_file)
   #   base_counts[cell_id] = nrow(res)
   # }

   coverage_rates = round(base_counts / total_genome_size * 100, 2)
   #names(coverage_rates) = gsub('.rmdup.nonCG_filtered.bam', '', names(coverage_rates))

   #hist(coverage_rates, breaks = 20)
  df_coverage_rates = data.frame(base_count = base_counts, coverage_rate = coverage_rates)
  print('Computed coverage rates.')
  return(df_coverage_rates)
}


merge_r1_and_r2_bam_for_cell <- function(r1_bam, r2_bam, merged_bam)
{
  sys_command = paste0('samtools merge -h '
                       , ' --samtools_path ', samtools_path
                       , ' ', out_bam
                       , ' ', r1_bam
                       , ' ', r2_bam
                       )
  system(sys_command)
}



merge_r1_and_r2_bam_for_sample <- function(alignment_dir_in,  alignment_dir_out)
{
  r1_bam_files = list.files(alignment_dir_in, pattern = 'R1')

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
    clusterExport(cl, ls(.GlobalEnv))
    registerDoSNOW(cl)
    clusterExport(cl, ls(.GlobalEnv))


    foreach(i=1:length(r1_bam_files)) %dopar%
    {
      r1_bam_file = r1_bam_files[i]
      r2_bam = paste0(alignment_dir_in,  gsub('R1', 'R2', r1_bam) )
      merged_bam = paste0(alignment_dir_out,  gsub('R1', '', r1_bam) )
      r1_bam = paste0(alignment_dir_in, r1_bam)

      merge_r1_and_r2_bam_for_cell(r1_bam, r2_bam, merged_bam)
    }#foreach

    stopCluster(cl)

  }else
  {

    for(r1_bam in r1_bam_files)
    {
      r2_bam = paste0(alignment_dir_in,  gsub('R1', 'R2', r1_bam) )
      merged_bam = paste0(alignment_dir_out,  gsub('R1', '', r1_bam) )
      r1_bam = paste0(alignment_dir_in, r1_bam)

      merge_r1_and_r2_bam_for_cell(r1_bam, r2_bam, merged_bam)
    }#for

  }####else

}#merge_r1_and_r2_bam_for_sample






