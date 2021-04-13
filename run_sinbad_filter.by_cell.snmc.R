#library(devtools)
#install_github("yasin-uzun/SINBAD.1.0")

#sample_name = '877476'
#sample_name = '878516'
#sample_name = '882304'
#working_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples//data/snmc//working_dir/LEUK/'
#read_type = 'R1'

library(SINBAD)
SINBAD::test()
packageVersion('SINBAD')

config_dir = '/mnt/isilon/tan_lab/uzuny/projects/sinbad/config_files/'
read_configs(config_dir)

print('**********************************')
print(paste('Duplicate remover: ', duplicate_remover))


args <- commandArgs(T)
working_dir = args[1]
sample_name = args[2]
read_type = args[3]
fastq_file = args[4]

#working_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples//data/snmc//working_dir/LEUK/'

source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/Main.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/preprocessing.R')

source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/alignment.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/methylation_calling.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/regional_quantification.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/dimensionality_reduction.R')




library(doSNOW)


sinbad_object = construct_sinbad_object(working_dir = working_dir,

                                        sample_name = sample_name)

print(sinbad_object)

log_dir= paste0(sinbad_object$main_log_dir, '/alignment/')
dir.create(log_dir, recursive = T)
setwd(alignment_dir)



  print(sinbad_object$alignment_dir)

  cat('** fastq_file: ', fastq_file, '\n')
  filter_cell(read_dir = sinbad_object$trimmed_fastq_dir,
             fastq_file = fastq_file,
             aligner = aligner,
             genomic_sequence_path = genomic_sequence_path,
             alignment_dir = sinbad_object$alignment_dir,
             log_dir = log_dir)

