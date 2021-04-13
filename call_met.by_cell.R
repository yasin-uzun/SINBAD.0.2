library(SINBAD)
SINBAD::test()
packageVersion('SINBAD')


source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/Main.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/preprocessing.R')

source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/alignment.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/methylation_calling.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/regional_quantification.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/dimensionality_reduction.R')




args <- commandArgs(T)
working_dir = args[1]
sample_name = args[2]
bam_file = args[3]
config_dir = args[4]

read_configs(config_dir)

#working_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples//data/snmc//working_dir/LEUK/'

library(doSNOW)



#sample_name <- 'Leuk_Batch_1'
sinbad_object = construct_sinbad_object(working_dir = working_dir,

                                        sample_name = sample_name)

print(sinbad_object)
#sinbad_object$alignment_dir = paste0(working_dir, sample_name, '/alignments_temp/')

sinbad_object$alignment_dir = gsub('alignments_r1_and_r2', 'alignments', sinbad_object$alignment_dir)

cell_id = gsub('.organism.bam', '', bam_file)

log_dir= paste0(sinbad_object$main_log_dir)
dir.create(log_dir, recursive = T)
setwd(alignment_dir)


call_methylation_sites_for_cell(alignment_dir = sinbad_object$alignment_dir,
                                methylation_calls_dir = sinbad_object$methylation_calls_dir,
                                cell_id = cell_id,
                                bme_param_settings = bme_param_settings,
                                log_dir = log_dir)

