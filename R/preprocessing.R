library(ShortRead)

get_fast_files <- function(fastq_dir, pattern = '')
{
  fastq_files_1 = list.files(fastq_dir, pattern = "\\.fastq.gz$")
  fastq_files_2 = list.files(fastq_dir, pattern = "\\.fastq$")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[grepl(pattern = pattern, fastq_files)]

  return(fastq_files)
}

get_r2_indeces_from_r1 <- function(r1_fastq_dir, r2_input_fastq_dir, r2_output_fastq_dir, sample_name = '')
{
  setwd(r1_fastq_dir)
  dir.create(r2_output_fastq_dir, recursive = T)

  r1_fastq_files = get_fast_files(fastq_dir = r1_fastq_dir, pattern = '_R1')
  r1_fastq_files = r1_fastq_files[!grepl(pattern = 'Undetermined', r1_fastq_files)]

  r2_fastq_files = get_fast_files(fastq_dir = r2_input_fastq_dir, pattern = '_R2')
  r2_fastq_files = r2_fastq_files[!grepl(pattern = 'Undetermined', r2_fastq_files)]

  r1_fastq_files = r1_fastq_files[grepl(sample_name, r1_fastq_files)]
  r2_fastq_files = r2_fastq_files[grepl(sample_name, r2_fastq_files)]


  r1_fastq_file_count = length(r1_fastq_files)




  if(num_cores > 1)
  {
    thread_count = min(r1_fastq_file_count, num_cores)
    cl <- makeCluster(thread_count, outfile="", type = 'SOCK')
    registerDoSNOW(cl)

    foreach(i=1:r1_fastq_file_count, .export = ls(globalenv()) ) %dopar%
    {
      r1_fastq_file = r1_fastq_files[i]
      print(r1_fastq_file)

      r2_fastq_file = gsub('_R1', '_R2', r1_fastq_file)

      command = paste0('perl ', perl_index_transfer_path, ' ',
                       ' --input_r1_fastq_file ',  r1_fastq_dir, '/', r1_fastq_file ,
                       ' --input_r2_fastq_file ', r2_input_fastq_dir, '/', r2_fastq_file ,
                       ' --demux_index_length ', demux_index_length,
                       ' --output_r2_fastq_file ', r2_output_fastq_dir, '/', r2_fastq_file )

      print(command)
      system(command)

    }#foreach(i=1:length(raw_fastq_files))

    stopCluster(cl)



  }else
  {

    for(i in 1:r1_fastq_file_count)
    {
      r1_fastq_file = r1_fastq_files[i]
      print(r1_fastq_file)

      r2_fastq_file = gsub('_R1', '_R2', r1_fastq_file)

      command = paste0('perl ', perl_index_transfer_path, ' ',
                       ' --input_r1_fastq_file ',  r1_fastq_dir, '/', r1_fastq_file ,
                       ' --input_r2_fastq_file ', r2_input_fastq_dir, '/', r2_fastq_file ,
                       ' --demux_index_length ', demux_index_length,
                       ' --output_r2_fastq_file ', r2_output_fastq_dir, '/', r2_fastq_file )

      print(command)
      system(command)

    }#foreach(i=1:length(raw_fastq_files))


  }



}#demux_fastq_files

demux_fastq_files <- function(raw_fastq_dir, demux_index_file, demux_index_length, demux_fastq_dir,
                              main_log_dir, read_type = '', sample_name = '')
{
  setwd(raw_fastq_dir)

  fastq_files_1 = list.files(raw_fastq_dir, pattern = "\\.fastq.gz$")
  fastq_files_2 = list.files(raw_fastq_dir, pattern = "\\.fastq$")

  raw_fastq_files = union(fastq_files_1, fastq_files_2)
  raw_fastq_files = raw_fastq_files[!grepl('Undetermined', raw_fastq_files)]

  raw_fastq_files = raw_fastq_files[grepl(sample_name, raw_fastq_files)]
  raw_fastq_files = raw_fastq_files[grepl(read_type, raw_fastq_files)]

  demux_log_dir = paste0(main_log_dir, '/demux/')
  dir.create(demux_log_dir, showWarnings = F, recursive = T)

  if(num_cores > 1)
  {
    raw_fastq_file_count = length(raw_fastq_files)
    thread_count = min(raw_fastq_file_count, num_cores)
    cl <- makeCluster(thread_count, outfile="", type = 'SOCK')
    registerDoSNOW(cl)
    foreach(i=1:raw_fastq_file_count, .export = ls(globalenv()) ) %dopar%
    {
      raw_fastq_file = raw_fastq_files[i]
      print(raw_fastq_file)
      output_prefix = gsub('.fastq.gz', '', raw_fastq_file)
      #ssdemux_log_file = paste0(demux_log_dir, '/', raw_fastq_file, '.log')

      demux_command = paste0('perl ', perl_demux_path, ' ',
                             ' --demux_index_file ',  demux_index_file ,
                             ' --raw_fastq_file ', raw_fastq_dir, '/', raw_fastq_file ,
                             ' --demux_index_length ', demux_index_length,
                             ' --output_dir ', demux_fastq_dir ,
                             ' --output_prefix ', output_prefix,
                             ' --log_dir ', demux_log_dir)

      print(demux_command)
      system(demux_command)

    }#foreach(i=1:length(raw_fastq_files))

    stopCluster(cl)

  }#if(num_cores > 1)



}#demux_fastq_files

read_demux_logs <- function(main_log_dir)
{

  demux_log_dir = paste0(main_log_dir, '/demux/')
  demux_log_files = list.files(demux_log_dir, pattern = "\\.log$")
  demux_log_files = demux_log_files[!grepl('Undetermined', demux_log_files)]
  setwd(demux_log_dir)
  list_df_demux_combined = list()

  for(demux_log_file in demux_log_files)
  {
    print(demux_log_file)
    lane_id = gsub('.log', '', demux_log_file)
    list_df_demux_combined[[lane_id]] = read.table(demux_log_file, header = T)
  }#for

  df_demux_combined = do.call('rbind', list_df_demux_combined)
  df_demux_combined = data.frame(Group = rownames(df_demux_combined), df_demux_combined)
  head(df_demux_combined)

  return(df_demux_combined)
}#read_demux_logs

#source("/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/preprocessing.R")

trim_fastq_files <- function(demux_fastq_dir, trimmed_fastq_dir, main_log_dir)
{

  print('Trimming fastq files')
  setwd(demux_fastq_dir)

  fastq_files_1 = list.files(demux_fastq_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(demux_fastq_dir, pattern = "*.fastq")

  demux_fastq_files = union(fastq_files_1, fastq_files_2)
  demux_fastq_files = demux_fastq_files[!grepl('No_matching_index', demux_fastq_files)]
  print('demux_fastq_files :')
  #print( demux_fastq_files)
  trimming_log_dir = paste0(main_log_dir, '/trimming/')
  dir.create(trimming_log_dir, showWarnings = F, recursive = T)

  if(num_cores > 1)
  {
    cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
    registerDoSNOW(cl)
    clusterExport(cl, ls(), envir = environment())



    print('*******************Trimming....')
    foreach(i=1:length(demux_fastq_files), .export = ls(globalenv())) %dopar%
    {
      fastq_file = demux_fastq_files[i]
      print(paste('Input fastq: ', fastq_file) )
      log_file = paste0(trimming_log_dir, fastq_file, '.log')

      if(trimmer == 'cutadapt')
      {
           command = paste0(cutadapt_path, '/cutadapt ',cutadapt_param_settings,
                             ' -o ', trimmed_fastq_dir, '/', fastq_file,
                             '  ',  demux_fastq_dir, '/',fastq_file ,
                             ' >  ',  log_file
                             )
          print(command)
          system(command)


      }else if(trimmer == 'trim_galore')
      {

        command = paste0(trim_galore_path, '/trim_galore ',trim_galore_param_settings,
                           ' -o ', trimmed_fastq_dir,
                           '  ',  demux_fastq_dir, '/',fastq_file ,
                           ' >  ',  log_file
                          )
        print(command)
        system(command)
        trimmed_fq_name = sub('.fastq.gz', '_trimmed.fq.gz',  fastq_file)
        mv_command = paste0('mv ', trimmed_fastq_dir, trimmed_fq_name, ' ',  trimmed_fastq_dir, fastq_file)
        system(mv_command)

      }else if(trimmer == 'Trimmomatic')
      {

        command = paste0('java -jar ', Trimmomatic_jar_path,
                         ' ', Trimmomatic_param_settings,
                         ' -trimlog ', trimming_log_dir, '/',log_file,
                         ' ', demux_fastq_dir, '/',fastq_file ,
                         ' ',  fastq_file
                         )

        system(command)


      }else
      {
        stop('Invalid  trimmer name, exiting. Trimmer must be one of these: cutadapt trim_galore Trimmomatic')
      }


      #ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25



    }#foreach(i=1:length(raw_fastq_files))

  }#if(num_cores > 1)



}


process_cutadapt_logs <- function(main_log_dir)
{
  cutadapt_log_dir = paste0(log_dir, '/cutadapt/')
  cutadapt_log_files = list.files(cutadapt_log_dir, pattern = "*.log")
  setwd(cutadapt_log_dir)

  row_names = c()
  result_list = list()
  for(report_file in cutadapt_log_files)
  {
    print(report_file)
    cell_id = gsub('.fastq.gz.log', '', report_file)
    lines = readLines(report_file)
    informative_lines = lines[grep(':', lines)]
    informative_lines = informative_lines[!grepl('strand', informative_lines)]
    informative_lines = gsub('\\(', ':', informative_lines)
    informative_lines = gsub('\\)', ':', informative_lines)

    df_report = read.delim(text = informative_lines, sep = '\t',  header = F)
    class(df_report)
    row_names = gsub(':', '', df_report$V1)
    result_list[[cell_id]] = as.character(df_report$V2)
  }

  row_names[1] = "Total_reads"
  row_names[2] = "Alignments"
  row_names[3] = "Alignment_rate"

  row_names[4] = "Sequences_with_no_alignments"
  row_names[6] = "Discarded_sequences"

  df_result = do.call(cbind, result_list)
  rownames(df_result) = row_names
  t_df_result =  t(df_result)
  final_table = data.frame(Cell_ID = rownames(t_df_result),  t_df_result, stringsAsFactors = F)

}

count_fastq_reads <- function(fastq_dir)
{
  print(paste('Reading fastq files from the directory:', fastq_dir))
  fastq_files_1 = list.files(fastq_dir, pattern = "*.fastq.gz")
  fastq_files_2 = list.files(fastq_dir, pattern = "*.fastq")

  fastq_files = union(fastq_files_1, fastq_files_2)
  fastq_files = fastq_files[!grepl('No_matching_index', fastq_files)]

  num_fastq = length(fastq_files)
  print(paste('I found:', num_fastq, ' fastq files: '))
  print(fastq_files)
  print(paste('Now counting reads inside fastqs '))


  cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
  registerDoSNOW(cl)
  clusterExport(cl, ls(), envir = environment())
  #clusterExport(cl, varlist = ls(), envir = environment())

  #read_counts = c()
  line_counts_list = foreach(i=1:length(fastq_files), .export = ls(globalenv())) %dopar%
  {
    library(ShortRead)
    fastq_file = fastq_files[i]
    print(fastq_file)
    countLines(fastq_dir, fastq_file)
  }

  line_counts = unlist(line_counts_list)
  read_counts = line_counts / 4
  names(read_counts) = gsub('.fastq.gz', '', names(read_counts))
  names(read_counts) = gsub('.fastq', '', names(read_counts))

  return(read_counts)

}

plot_preprocessing_results <- function(sample_name, demux_reports, demux_read_counts, trimmed_read_counts)
{
  #par(mfrow = c(2,2))

  layout(mat = matrix(c(1, 2, 1, 3),
                      2, 2, byrow = TRUE))

  par(mar = c(6,6,3,8))


  demux_reports_sorted = demux_reports[order(demux_reports$Total_reads), ]
  demux_matrix = t(as.matrix(demux_reports_sorted[,2:3])) / 1000000
  barplot(demux_matrix, col = c('aquamarine', 'brown1'), las = 2, horiz = T,
          main = 'Demultiplexing',  xlab = 'million reads', ylab = '')

  legend('bottomright', legend = c('With index', 'No index'), fill = c('aquamarine', 'brown1'), bty = "n")

  roof = max(demux_read_counts, trimmed_read_counts) / 1000000
  x_labels = c('Demux', 'Trimmed')
  boxplot(
    demux_read_counts/1000000,
    trimmed_read_counts/1000000,
    #ylim = c(0, roof),
    col = "#8dd3c7", names = x_labels,
    ylab = 'million reads',
    main = 'Read statistics'
    , outpch = 19, outcex = 0.2
    , las = 2
    #,cex.axis = 1.25, cex.names = 1.25,
   #cex.main = 1.25, cex.lab = 1.25
    #, outline =F
  )    #Plot number of reads


  loss_rates = (1 - trimmed_read_counts / demux_read_counts[names(trimmed_read_counts)] ) * 100

  #barplot(sort(loss_rates))
  plot(sort(loss_rates), type = 'l',
       xlab = 'Cells', ylab = '%Reads filtered out', col = 'blue',
       main = 'Adapter trimming'
       #, cex.main = 0.8
       )


  title(paste(sample_name, ': Preprocessing results'), line = -1, outer = TRUE)



}


