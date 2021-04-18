library(doSNOW)


test <- function()
{
  cat('SINBAD installation is ok.\n')
}
#

read_configs <- function(config_dir)
{
  print('Reading config.general.R for program paths...')
  source(paste0(config_dir, '/config.general.R') )
  print('Reading config.genome.R for genome paths...')
  source(paste0(config_dir, '/config.genome.R') )
  print('Reading config.project.R for project settings...')
  source(paste0(config_dir, '/config.project.R') )
  #system('echo $PATH')

}



#sample_name = 'Test'
#readRenviron(path = 'variables.env')
#system('source ~/.bashrc')

#image_file = paste0(working_dir, '/Image_2020_06_03.img')
#image_file = paste0(working_dir, '/Image_2020_06_11.img')
#save.image(image_file)

construct_sinbad_object <- function(raw_fastq_dir  = NA,
                                    demux_index_file = NA,
                                    working_dir ,
                                    sample_name)
{

  sample_working_dir = paste0(working_dir, '/', sample_name, '/')
  dir.create(sample_working_dir, recursive = T)

  main_log_dir <<- paste0(sample_working_dir, '/logs/')
  demux_fastq_dir <<- paste0(sample_working_dir, '/demux_fastq/')
  trimmed_fastq_dir <<- paste0(sample_working_dir, '/trimmed_fastq/')
  alignment_dir <<- paste0(sample_working_dir, '/alignments/')

  if(sequencing_type == 'paired')
  {
    alignments_r1_and_r2_dir <<- paste0(sample_working_dir, '/alignments_r1_and_r2/')
    r2_meta_fastq_dir <<- paste0(sample_working_dir, '/r2_meta_fastq/')

    dir.create(alignments_r1_and_r2_dir, showWarnings = F, recursive = T)
    dir.create(r2_meta_fastq_dir, showWarnings = F, recursive = T)

  }

  methylation_calls_dir <<- paste0(sample_working_dir, '/methylation_calls/')
  summary_dir <<- paste0(sample_working_dir, '/summary/')
  plot_dir <<- paste0(sample_working_dir, '/plots/')
  matrix_dir <<- paste0(sample_working_dir, '/matrices/')
  object_dir <<- paste0(sample_working_dir, '/objects/')
  alignment_summary_file <<- paste0(sample_working_dir, '/r2_meta_fastq/')

  if(trimmer == "NoTrimming")
  {
    trimmed_fastq_dir = demux_fastq_dir
  }


  dir.create(main_log_dir, showWarnings = F, recursive = T)
  dir.create(trimmed_fastq_dir, showWarnings = F, recursive = T)
  dir.create(alignment_dir, showWarnings = F, recursive = T)
  dir.create(methylation_calls_dir, showWarnings = F, recursive = T)
  dir.create(summary_dir, showWarnings = F, recursive = T)
  dir.create(plot_dir, showWarnings = F, recursive = T)
  dir.create(demux_fastq_dir, showWarnings = F, recursive = T)
  dir.create(matrix_dir, showWarnings = F, recursive = T)
  dir.create(object_dir, showWarnings = F, recursive = T)

  sinbad_object = list('sample_working_dir' = sample_working_dir,
                       'main_log_dir' = main_log_dir,
                       'raw_fastq_dir' = raw_fastq_dir,
                       'demux_index_file' = demux_index_file,
                       'demux_fastq_dir' = demux_fastq_dir,
                       'trimmed_fastq_dir' =  trimmed_fastq_dir,
                       'alignment_dir' = alignment_dir,
                       'methylation_calls_dir'  = methylation_calls_dir,
                       'summary_dir' = summary_dir,
                       'plot_dir' = plot_dir,
                       'sample_name' = sample_name,
                       'matrix_dir' = matrix_dir,
                       'object_dir' = object_dir,
                       'r2_meta_fastq_dir' = r2_meta_fastq_dir

                       )

  class(sinbad_object) = 'Sinbad'

  return(sinbad_object)


}


wrap_demux_fastq_files <- function(sinbad_object, flag_r2_index_embedded_in_r1_reads = FALSE)
{
  #Demultiplex fastq files
  #TODO: Check the input fastq dir and demux file exists, give error and stop otherwise


  if(sequencing_type == 'paired' )
  { #Pair ended
    demux_fastq_files(sinbad_object$raw_fastq_dir,
                      sinbad_object$demux_index_file,
                      demux_index_length,
                      sinbad_object$demux_fastq_dir,
                      sinbad_object$main_log_dir,
                      sinbad_object$sample_name,
                      read_type = 'R1')

    if(flag_r2_index_embedded_in_r1_reads)
    {

      demux_fastq_files(sinbad_object$r2_meta_fastq_dir,
                      sinbad_object$demux_index_file,
                      demux_index_length,
                      sinbad_object$demux_fastq_dir,
                      sinbad_object$main_log_dir,
                      sinbad_object$sample_name,
                      read_type = 'R2')
    }else{
      demux_fastq_files(sinbad_object$raw_fastq_dir,
                        sinbad_object$demux_index_file,
                        demux_index_length,
                        sinbad_object$demux_fastq_dir,
                        sinbad_object$main_log_dir,
                        sinbad_object$sample_name,
                        read_type = 'R2')
    }# if(flag_r2_index_embedded_in_r1_reads)
  }else
  { #single ended
    demux_fastq_files(sinbad_object$raw_fastq_dir,
                      sinbad_object$demux_index_file,
                      demux_index_length,
                      sinbad_object$demux_fastq_dir,
                      sinbad_object$main_log_dir,
                      sinbad_object$sample_name)
  }

  return(sinbad_object)

}

wrap_demux_stats <- function(sinbad_object)
{

  sinbad_object$df_demux_reports = read_demux_logs(sinbad_object$main_log_dir)

  demux_summary_file = paste0(summary_dir, '/Demux_statistics.tsv')
  write.table(sinbad_object$df_demux_reports,
              file = demux_summary_file,
              sep = '\t', quote = F,
              row.names = F, col.names = T)

  #Count demuxd reads
  sinbad_object$demux_read_counts =  count_fastq_reads(sinbad_object$demux_fastq_dir)

  return(sinbad_object)

}

wrap_trim_fastq_files <- function(sinbad_object)#Trim adapters
{
  trim_fastq_files(sinbad_object$demux_fastq_dir,
                   sinbad_object$trimmed_fastq_dir,
                   sinbad_object$main_log_dir)

  return(sinbad_object)
}

wrap_trim_stats <- function(sinbad_object)#Trim adapters
{
  sinbad_object$trimmed_read_counts = count_fastq_reads(sinbad_object$trimmed_fastq_dir)



  return(sinbad_object)
}

wrap_plot_preprocessing_stats <- function(sinbad_object)
{
  plot_file = paste0(sinbad_object$plot_dir, '/Preprocessing_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = T, title = sample_name)
  plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/Preprocessing_statistics.png')
  png(plot_file, widt = 800, height = 600)
  plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()
}

wrap_align_sample <- function(sinbad_object, pattern = '')
{
  pbat_flag = 0

  #Run aligner
  if(sequencing_type == 'paired')
  {
    if(grepl('--pbat', bismark_aligner_param_settings) ) {pbat_flag = 1}
  }


  if(pbat_flag == 0)
  {
    align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
               genomic_sequence_path = genomic_sequence_path,
               alignment_dir = sinbad_object$alignment_dir,
               aligner = aligner,
               num_cores= num_cores,
               mapq_threshold =mapq_threshold,
               main_log_dir = sinbad_object$main_log_dir,
               pattern = pattern)
  }else
  {
    print('Bismark-paired-pbat')
    #align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
     #            genomic_sequence_path = genomic_sequence_path,
    #             alignment_dir = sinbad_object$alignment_dir,
     #            aligner = aligner,
      #           num_cores= num_cores,
       #          mapq_threshold =mapq_threshold,
        #         main_log_dir = sinbad_object$main_log_dir,
         #        pattern = 'R1')

    bismark_aligner_param_settings = gsub('--pbat', '', bismark_aligner_param_settings)

    align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
                 genomic_sequence_path = genomic_sequence_path,
                 alignment_dir = sinbad_object$alignment_dir,
                 aligner = aligner,
                 num_cores= num_cores,
                 mapq_threshold =mapq_threshold,
                 main_log_dir = sinbad_object$main_log_dir,
                 pattern = 'R2')

  }


}



wrap_generate_alignment_stats <- function(sinbad_object)
{
  sinbad_object$df_alignment_reports = process_bismark_alignment_reports(alignment_dir = sinbad_object$alignment_dir)
  sinbad_object$df_bam_read_counts = count_bam_files(alignment_dir = sinbad_object$alignment_dir)
  dim(sinbad_object$df_alignment_reports)
  dim(sinbad_object$df_bam_read_counts)


  sinbad_object$df_alignment_stats = base::merge(sinbad_object$df_alignment_reports,
                                                 sinbad_object$df_bam_read_counts, by = 0)
  dim(sinbad_object$df_alignment_stats)
  head(sinbad_object$df_alignment_stats)


  return(sinbad_object)

}

wrap_compute_coverage_rates <- function(sinbad_object)
{
  df_coverage_rates = compute_coverage_rates(alignment_dir = sinbad_object$alignment_dir)
  df_coverage_rates_ordered = df_coverage_rates[as.character(sinbad_object$df_alignment_stats$Cell_ID), ]

  sinbad_object$df_coverage_rates =df_coverage_rates_ordered
  sinbad_object$df_alignment_stats$base_count = df_coverage_rates_ordered[,1]
  sinbad_object$df_alignment_stats$coverage_rate = df_coverage_rates_ordered[,2]



  return(sinbad_object)
}


merge_r1_and_r2_alignment_stats <- function(df_alignment_stats)
{
  df_alignment_stats$Row.names = NULL
  stat_columns = colnames(df_alignment_stats)
  rate_columns = stat_columns[grepl('rate', stat_columns) | grepl('^C\\.', stat_columns)]
  count_columns = setdiff(stat_columns, rate_columns)
  df_count_stats = df_alignment_stats[, count_columns]

  cell_id_wo_r   = df_count_stats$Cell_ID
  cell_id_wo_r = gsub('_R1', '', cell_id_wo_r)
  cell_id_wo_r = gsub('_R2', '', cell_id_wo_r)

  df_count_stats$Cell_ID = cell_id_wo_r

  for(count_column in setdiff(count_columns, 'Cell_ID' ) )
  {
    df_count_stats[, count_column] = as.integer(df_count_stats[, count_column])
  }

  attach(df_count_stats)
  df_count_stats_sum = aggregate(as.matrix(df_count_stats[, 2:ncol(df_count_stats)]),
                                 by=list(Cell_ID), FUN = sum, na.rm = TRUE)
  detach(df_count_stats)
  colnames(df_count_stats_sum)[1] = 'Cell_ID'

  df_count_stats_sum$Alignment_rate  = round(df_count_stats_sum$Alignments / df_count_stats_sum$Total_reads * 100)
  df_chrom_sizes = read.table(chrom_sizes_file, header = F)
  genome_size = sum(df_chrom_sizes$V2)

  df_count_stats_sum$coverage_rate  = -1

  attach(df_count_stats_sum)
  df_count_stats_sum$C.methylated.in.CpG.context = round(100 * Total.methylated.C.s.in.CpG.context / (Total.methylated.C.s.in.CpG.context + Total.unmethylated.C.s.in.CpG.context), 1)
  df_count_stats_sum$C.methylated.in.CHG.context = round(100 * Total.methylated.C.s.in.CHG.context / (Total.methylated.C.s.in.CHG.context + Total.unmethylated.C.s.in.CHG.context), 1)
  df_count_stats_sum$C.methylated.in.CHH.context = round(100 * Total.methylated.C.s.in.CHH.context / (Total.methylated.C.s.in.CHH.context + Total.unmethylated.C.s.in.CHH.context), 1)

  detach(df_count_stats_sum)

  rownames(df_count_stats_sum) = df_count_stats_sum$Cell_ID
  head(df_count_stats_sum)



  ###Scatter plot filtering
  filter_aln_rate =  df_count_stats_sum$Alignment_rate  > alignment_rate_threshold
  #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
  filter_read_count = df_count_stats_sum$nc_filtered_read_counts > minimum_filtered_read_count

  passing = filter_aln_rate &  filter_read_count

  df_count_stats_sum$pass = 0
  df_count_stats_sum$pass[passing] = 1

  return(df_count_stats_sum)
}

wrap_plot_alignment_stats <- function(sinbad_object)
{
  ###Scatter plot filtering
  filter_aln_rate =  sinbad_object$df_alignment_stats$Alignment_rate  > alignment_rate_threshold
  #filter_read_count = df_alignment_stats$organism_read_counts > organism_minimum_filtered_read_count
  filter_read_count = sinbad_object$df_alignment_stats$nc_filtered_read_counts > minimum_filtered_read_count

  passing = filter_aln_rate &  filter_read_count

  sinbad_object$df_alignment_stats$pass = 0
  sinbad_object$df_alignment_stats$pass[passing] = 1

  sinbad_object$df_alignment_stats$Row.names = NULL
  rownames(sinbad_object$df_alignment_stats) = sinbad_object$df_alignment_stats$Cell_ID
  sinbad_object$alignment_summary_file = paste0(sinbad_object$summary_dir, '/Alignment_statistics.tsv')
  write.table(sinbad_object$df_alignment_stats,
              file = sinbad_object$alignment_summary_file,
              sep = '\t', quote = F, row.names = F, col.names = T)

  dir.create(paste0(sinbad_object$plot_dir, '/QC/') )

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Alignment_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = F, title = sinbad_object$sample_name)
  plot_alignment_stats(sample_name = sinbad_object$sample_name,
                       df_alignment_stats = sinbad_object$df_alignment_stats)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Alignment_statistics.png')
  png(plot_file, width = 600, height = 800)
  plot_alignment_stats(sinbad_object$sample_name, sinbad_object$df_alignment_stats)
  dev.off()

  return(sinbad_object)
}

wrap_merge_r1_and_r2_bam <- function(sinbad_object)
{
  alignment_dir_01 = sinbad_object$alignment_dir
  alignment_dir_02 = gsub('alignments', 'alignments_r1_and_r2', alignment_dir_01)
  command = paste('mv', alignment_dir_01, alignment_dir_02)
  system(command)
  command = paste('mkdir', alignment_dir_01)
  system(command)

  merge_r1_and_r2_bam_for_sample(alignment_dir_02, alignment_dir_01)
}

wrap_call_methylation_sites <- function(sinbad_object)
{
  call_methylation_sites_for_sample(sinbad_object$alignment_dir,
                                    sinbad_object$methylation_calls_dir,
                                    sinbad_object$main_log_dir,
                                    bme_param_settings )

  return(sinbad_object)


}

wrap_generate_methylation_stats <- function(sinbad_object)
{

  sinbad_object$df_org_split_reports = process_bismark_split_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                     genome_type = 'organism')


  sinbad_object$list_org_bias_reports = process_bismark_bias_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                     genome_type = 'organism')

  sinbad_object$df_lambda_split_reports = NULL
  sinbad_object$list_lambda_bias_reports = NULL

  if(exists('lambda_control') )
  {
    if(lambda_control)
    {
      sinbad_object$df_lambda_split_reports = process_bismark_split_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                            genome_type = 'lambda')

      sinbad_object$list_lambda_bias_reports = process_bismark_bias_reports(methylation_calls_dir =  sinbad_object$methylation_calls_dir,
                                                                            genome_type = 'lambda')
    }
  }

  list_met_call_counts = list()
  met_types = c('CpG', 'CHG', 'CHH')

  for(met_type in met_types)
  {
    list_met_call_counts[[met_type]] = get_met_call_counts(methylation_calls_dir, met_type)
  }

  sinbad_object$list_met_call_counts = list_met_call_counts

  return(sinbad_object)
}

wrap_plot_met_stats <- function(sinbad_object)
{
  dir.create(paste0(sinbad_object$plot_dir, '/QC/') )

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Met_call_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = F, title = sinbad_object$sample_name)
  plot_split_reports(df_org_split_reports = sinbad_object$df_org_split_reports,
                     df_lambda_split_reports = sinbad_object$df_lambda_split_reports,
                     list_org_bias_reports = sinbad_object$list_org_bias_reports,
                     list_met_call_counts = sinbad_object$list_met_call_counts,
                     lambda_flag = lambda_control)
  dev.off()

  plot_file = paste0(sinbad_object$plot_dir, '/QC/Met_call_statistics.png')
  png(plot_file, width = 600, height = 800)
  plot_split_reports(sinbad_object$df_org_split_reports,
                     sinbad_object$df_lambda_split_reports,
                     sinbad_object$list_org_bias_reports,
                     list_met_call_counts = sinbad_object$list_met_call_counts,
                     lambda_flag = lambda_control)
  dev.off()

  return(sinbad_object)

}

wrap_read_annots <- function(sinbad_object)
{
  #Read regions
  print('Reading gene coordinates')
  annot_genes = read_region_annot(gene_annot_file, format_file)
  head(annot_genes)
  annot_genes_2kb = annot_genes
  annot_genes_2kb$start = annot_genes$start - 2000
  annot_genes_2kb$start[annot_genes_2kb$start < 0 ] = 0
  annot_genes_2kb$end = annot_genes$end + 2000

  print('Computing promoters')
  annot_promoters = get_promoters(df_gene_annot = annot_genes)
  head(annot_promoters)
  print('Reading 100k bins')
  annot_100k_bins = read_region_annot(bins_100k_file, format_file)
  head(annot_100k_bins)
  print('Reading 10k bins')
  annot_10k_bins = read_region_annot(bins_10k_file, format_file)
  head(annot_10k_bins)

  sinbad_object$annot_list = list()

  sinbad_object$annot_list[['100k_bins']] = annot_100k_bins
  sinbad_object$annot_list[['10k_bins']] = annot_10k_bins
  sinbad_object$annot_list[['genes']] = annot_genes
  sinbad_object$annot_list[['genes_2kb']] = annot_genes_2kb

  sinbad_object$annot_list[['promoters']] = annot_promoters

  return(sinbad_object)

}


wrap_quantify_regions <- function(sinbad_object, methylation_type = 'CpG')
{


  if(  !('annot_list' %in% names(sinbad_object) )  | length(sinbad_object$annot_list) == 0 )
  {
    print('Annotation regions are not computed. Call wrap_read_annots() function first. Exiting.')
    return(sinbad_object)
  }

  sinbad_object$met_matrices = list()
  attach(sinbad_object$df_alignment_stats)
  exclude_cells = rownames(sinbad_object$df_alignment_stats[pass == 0, ])
  detach(sinbad_object$df_alignment_stats)

  sinbad_object$list_of_list_call_count_matrices = list()
  sinbad_object$list_aggr_rates = list()

  for(annot_name in names(sinbad_object$annot_list))
  {
    print(annot_name)

    list_call_count_matrices = compute_call_count_matrices(df_region = sinbad_object$annot_list[[annot_name]],
                                                 methylation_calls_dir = sinbad_object$methylation_calls_dir,
                                                 methylation_type = methylation_type,
                                                 exclude_cells = exclude_cells)

    aggr_rate = compute_aggr_met_rate(list_call_count_matrices)

    sinbad_object$list_of_list_call_count_matrices[[methylation_type]][[annot_name]] = list_call_count_matrices
    sinbad_object$list_aggr_rates[[methylation_type]][[annot_name]] = aggr_rate

    met_mat = compute_region_met_matrix(list_call_count_matrices = list_call_count_matrices,
                                       min_call_count_threshold = min_call_count_threshold)
    sinbad_object$met_matrices[[annot_name]] = met_mat


    met_rate_file = paste0(sinbad_object$matrix_dir, 'met_rate_mat.' , annot_name, '.',methylation_type ,'.rds')
    saveRDS(met_mat, file = met_rate_file)

    aggr_rate_file = paste0(sinbad_object$matrix_dir, 'aggr_rate.' , annot_name, '.',methylation_type ,'.rds')
    saveRDS(aggr_rate, file = aggr_rate_file)

    met_count_mat_file = paste0(sinbad_object$matrix_dir, 'met_count_mat.' , annot_name, '.',methylation_type ,'.rds')
    saveRDS(list_call_count_matrices$full_met_call_count_matrix, file = met_count_mat_file)

    total_count_mat_file = paste0(sinbad_object$matrix_dir, 'total_count_mat.' , annot_name, '.',methylation_type ,'.rds')
    saveRDS(list_call_count_matrices$full_total_call_count_matrix, file = total_count_mat_file)


  }#for

  return(sinbad_object)

}



wrap_impute_nas <- function(sinbad_object,  max_ratio_of_na_cells = max_ratio_of_na_cells)
{
  sinbad_object$imputed_matrices = list()

  for(annot_name in names(sinbad_object$met_matrices))
  {
    met_mat = sinbad_object$met_matrices[[annot_name]]
    imputed_matrix = impute_nas(met_mat = met_mat
                                , max_ratio_of_na_cells = max_ratio_of_na_cells)

    imputed_file = paste0(sinbad_object$matrix_dir, 'imputed.' , annot_name, '.rds')
    saveRDS(imputed_matrix, file = imputed_file)

    sinbad_object$imputed_matrices[[annot_name]] = imputed_matrix
  }#for

  return(sinbad_object)

}



wrap_dim_red <- function(sinbad_object,
                         annot_type = '100k_bins',
                         legend_title = 'Methylation',
                         methylation_type = 'CpG')
{
  #Reduce dimensionality
  #marker_genes = get_marker_genes('leuk')

  if( ! ('dim_red_objects' %in% names(sinbad_object)  )   )
  {
    sinbad_object$dim_red_objects = list()
  }
  name_for_dim_red = annot_type

  dim_red_object = reduce_dims_for_sample(
                         met_mat_for_dim_red = sinbad_object$met_matrices[[annot_type]]  ,
                         name_for_dim_red = name_for_dim_red  ,
                         plot_dir = sinbad_object$plot_dir  ,
                         methylation_type = methylation_type
                         )

  sinbad_object$dim_red_objects[[annot_type]] = dim_red_object

  title = paste0(sample_name, ' - ' , methylation_type ,
                 '\nDR region: ', annot_type
                 #, '\nFeature region: ', name_for_features
                 )

  dir.create(paste0(plot_dir, '/DimRed/') )

  gg1 = plot_dim_red(dim_red_object$umap, title = title )
  print(gg1)

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')


  clusters = compute_clusters(dim_red_object$umap)


  gg1 = plot_dim_red(dim_red_object$umap, title = title, groups = clusters )
  print(gg1)

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.clusters.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  plot_file = paste0(plot_dir, '/DimRed/UMAP.',name_for_dim_red,'.clusters.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')


  return(sinbad_object)

}


wrap_plot_features <- function(sinbad_object, features, dim_red_annot_type = '100k_bins', features_annot_type = 'promoters')
{

  dim_red_object = sinbad_object$dim_red_objects[[dim_red_annot_type]]

  feature_matrix = sinbad_object$met_matrices[[features_annot_type]]
  feature_matrix[1:5, 1:5]
  feature_matrix['PAX5', 1:5]

  legend_title = ''
  if(features_annot_type == 'promoters')
  {
    met_mat_for_features = 1 - sinbad_object$met_mat_for_features
    legend_title = 'De-methylation'
  }
  if(features_annot_type == 'MAPLE')
  {
    met_mat_for_features = sinbad_object$met_mat_for_features
    legend_title = 'MAPLE'
  }

  plot_features(umap = dim_red_object$umap,
                feature_matrix = feature_matrix,
                features = features,
                name_for_dim_red = dim_red_annot_type,
                name_for_features = features_annot_type,
                legend_title = legend_title)


}


wrap_dmr_analysis <- function(sinbad_object)
{
  #DM Analysis

  sinbad_object$dm_stat_list_for_clusters = dm_stat_test_for_clusters(sinbad_object$dim_red_list)
  sinbad_object$dm_result_file = paste0(sinbad_object$summary_dir, '/DM_Analysis.xlsx')
  write.xlsx(sinbad_object$dm_stat_list_for_clusters$dm_result_list_with_summary,
             file = sinbad_object$dm_result_file)

  sinbad_object$all_dms = do.call("rbind", sinbad_object$dm_stat_list_for_clusters$dm_result_list_with_summary[2:15] )
  sinbad_object$all_dms_sorted = all_dms[order(sinbad_object$all_dms$p.value), ]
  print(head(sinbad_object$all_dms_sorted))

  sinbad_object$all_dms_sig = sinbad_object$all_dms_sorted[sinbad_object$all_dms_sorted$adjusted.p.value < dmr_adj_p_value_cutoff, ]
  sinbad_object$heatmap_dm_regions = as.character(sinbad_object$all_dms_sorted$region)

  if(nrow(sinbad_object$all_dms_sig) > dm_num_heatmap_regions)
  {
    sinbad_object$heatmap_dm_regions = sinbad_object$heatmap_dm_regions[1:dm_num_heatmap_regions]
  }
  #sig_regions = dm_stat_list_for_clusters$sig_ids
  sinbad_object$sorted_clusters = sort(sinbad_object$dim_red_list$clusters)

  sinbad_object$met_mat = replace_nas_by_column_mean(sinbad_object$dim_red_list$met_mat_for_features)
  sinbad_object$sig_mat = as.matrix(met_mat[sinbad_object$heatmap_dm_regions, names(sinbad_object$sorted_clusters) ])
  sinbad_object$df_annot = data.frame(cluster = sinbad_object$sorted_clusters)

  head(sinbad_object$df_annot)

  ph <- pheatmap::pheatmap(sinbad_object$sig_mat,
                           cluster_cols = F,
                           cluster_rows = T,
                           show_colnames = F,
                           fontsize = 13,
                           annotation_col = sinbad_object$df_annot,
                           color = viridis(100),
                           #annotation_colors = ann_colors,
                           fontsize_row = 9)

  plot_file = paste0(sinbad_object$plot_dir, '/Heatmap.DMR_',sinbad_object$dim_red_list$name_for_features,'.clusters.eps')
  ggsave(ph, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  return(sinbad_object)


}

process_with_Seurat <- function(sinbad_object, vmr_count = 2000, max_ratio_of_na_cells = 0.75)
{


  library(Seurat)

  if(sinbad_object$name_for_features == 'Promoters')
  {
    met_mat_for_features = 1 - sinbad_object$met_mat_for_features
    legend_title = 'De-methylation'
  }
  if(sinbad_object$name_for_features == 'MAPLE')
  {
    met_mat_for_features = sinbad_object$met_mat_for_features
    legend_title = 'MAPLE'
  }

  met_mat_for_dim_red = impute_nas(sinbad_object$met_mat_for_dim_red)
  met_mat_for_dim_red = impute_nas(met_mat_for_features)

  seurat_object_for_umap = CreateSeuratObject(met_mat_for_dim_red)

  seurat_object_for_umap <- NormalizeData(seurat_object_for_umap, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object_for_umap <- FindVariableFeatures(seurat_object_for_umap, selection.method = "vst", nfeatures = vmr_count)
  all.genes <- rownames(seurat_object_for_umap)
  seurat_object_for_umap <- ScaleData(seurat_object_for_umap, features = all.genes)
  seurat_object_for_umap <- RunPCA(seurat_object_for_umap, features = VariableFeatures(object = seurat_object_for_umap))

  DimPlot(seurat_object_for_umap, reduction = "pca")

  seurat_object_for_umap <- FindNeighbors(seurat_object_for_umap, dims = 1:10)
  seurat_object_for_umap <- FindClusters(seurat_object_for_umap, resolution = 0.1)

  seurat_object_for_umap <- RunUMAP(seurat_object_for_umap, dims = 1:10, min.dist = 0.01)


  gg1 = DimPlot(seurat_object_for_umap, reduction = "umap")

  print(gg1)
  plot_file = paste0(sinbad_object$plot_dir, '/Seurat.UMAP.DR_',sinbad_object$name_for_dim_red,'.clusters.eps')
  ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

  DimPlot(seurat_object_for_umap, reduction = "umap")
  plot_file = paste0(sinbad_object$plot_dir, '/Seurat.UMAP.DR_',sinbad_object$name_for_dim_red,'.clusters.png')
  ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')


  seurat_object_for_features = CreateSeuratObject(met_mat_for_features)

  seurat_object_for_features <- NormalizeData(seurat_object_for_features,
                                              normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object_for_features <- ScaleData(seurat_object_for_features, features = all.genes)

  seurat_object_for_features@reductions = seurat_object_for_umap@reductions
  seurat_object_for_features@meta.data =  seurat_object_for_umap@meta.data

  #FeaturePlot(seurat_object_for_umap, features = 'CD3G')

  features = get_marker_genes('leuk')
  feature_plot_dir = paste0(sinbad_object$plot_dir, '/seurat_regions_',sinbad_object$name_for_dim_red,'/')
  dir.create(feature_plot_dir)
  dir.create(paste0(feature_plot_dir, '/eps/'))
  dir.create(paste0(feature_plot_dir, '/png/'))


  counts_genes = rownames(seurat_object_for_features@assays$RNA@data)
  head(counts_genes)
  features       = intersect(features, counts_genes)

  if( length(features) > 0 )
  {
    for(feature in features)
    {
      #print(feature)
      #cairo_ps(plot_file, fallback_resolution = 2400)
      #postscript(plot_file, onefile = F, width = 7, height = 6)
      gg1 = FeaturePlot(seurat_object_for_umap, features = feature)
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', feature, '.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', feature, '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 20, height = 20, units = 'cm')
      #dev.off()
    }#for(gene in marker_gene s)

  }

  n = length(features)
  if(n <= 10)
  {
    gg1 = StackedVlnPlot(obj = seurat_object_for_features, features = features, group.by = 'seurat_clusters') #, group.by = 'cell_type', assay = 'SCT')
    print(gg1)
    plot_file = paste0(feature_plot_dir, '/eps/', 'Violin', '.eps')
    ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
    plot_file = paste0(feature_plot_dir, '/png/', 'Violin', '.png')
    ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
  }else
  {
    for(i in 1:ceiling(n/10))
    {
      start = (i-1)*10 + 1
      end = min(i*10, n)
      gg1 = StackedVlnPlot(obj = seurat_object_for_features,
                           features = features[start:end],
                           group.by = 'seurat_clusters', y.max = 0.5) #, group.by = 'cell_type', assay = 'SCT')
      print(gg1)
      plot_file = paste0(feature_plot_dir, '/eps/', 'Violin.',i ,'.eps')
      ggsave(gg1, filename = plot_file, device = 'eps', width = 10, height = 28, units = 'cm')
      plot_file = paste0(feature_plot_dir, '/png/', 'Violin.',i, '.png')
      ggsave(gg1, filename = plot_file, device = 'png', width = 10, height = 28, units = 'cm')
    }



  }#else

  seurat_object_file = paste0(sinbad_object$plot_dir, 'seurat_object.',sinbad_object$name_for_features, '.rds')
  saveRDS(seurat_object_for_features, seurat_object_file )

  return(seurat_object_for_features)

}

process_sample_wrapper <- function(sinbad_object)
{


  par(mfrow = c(2,2))

  sinbad_object = wrap_demux_fastq_files(sinbad_object)

  sinbad_object = wrap_demux_stats(sinbad_object)

  sinbad_object = wrap_trim_fastq_files(sinbad_object)

  sinbad_object = wrap_trim_stats(sinbad_object)

  sinbad_object = wrap_align_sample(sinbad_object)

  sinbad_object = wrap_generate_alignment_stats(sinbad_object)

  sinbad_object = wrap_call_methylation_sites(sinbad_object)

  sinbad_object = wrap_generate_methylation_stats(sinbad_object)

  sinbad_object = wrap_read_regions(sinbad_object)

  sinbad_object = wrap_quantify_regions(sinbad_object)

  sinbad_object = wrap_dim_red(sinbad_object)

  sinbad_object = wrap_dmr_analysis(sinbad_object)


  return(sinbad_object)


}


