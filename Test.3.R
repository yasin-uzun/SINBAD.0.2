#library(devtools)
#install_github("yasin-uzun/SINBAD.1.0")

library(SINBAD)
SINBAD::test()

config_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/config_files/'
read_configs(config_dir)

example_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/data/example/'
raw_fastq_dir = paste0(example_dir, '/reads/')
demux_index_file = paste0(example_dir, '/demux_index.txt')
working_dir = paste0(example_dir, '/working_dir/')

library(doSNOW)

source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/Main.R')
source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/preprocessing.R')

source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/alignment.R')
source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/methylation_calling.R')
source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/regional_quantification.R')
source('/mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/R/dimensionality_reduction.R')


sample_name <- 'Test'
sinbad_object = construct_sinbad_object(working_dir = working_dir,
                                        raw_fastq_dir = raw_fastq_dir,
                                        demux_index_file = demux_index_file,
                                        sample_name = sample_name)

#process_sample_wrapper(sinbad_object)

#Demux
sinbad_object = wrap_demux_fastq_files(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.01.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)

#Demux stats
sinbad_object = wrap_demux_stats(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.02.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)


#Trim
sinbad_object = wrap_trim_fastq_files(sinbad_object)
sinbad_object_file = paste0(working_dir, '/sinbad_object.03.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)

#Demux
sinbad_object = wrap_trim_stats(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.04.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)



#Align
sinbad_object = wrap_align_sample(sinbad_object)
sinbad_object_file = paste0(working_dir, '/sinbad_object.05.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)


#Align stats
sinbad_object = wrap_generate_alignment_stats(sinbad_object)
sinbad_object_file = paste0(working_dir, '/sinbad_object.06.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)



sinbad_object = wrap_call_methylation_sites(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.07.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)




sinbad_object = wrap_generate_methylation_stats(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.08.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)







sinbad_object = wrap_read_regions(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.09.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)





sinbad_object = wrap_quantify_regions(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.10.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)





sinbad_object = wrap_dim_red(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.11.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)




sinbad_object = wrap_dmr_analysis(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.12.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print(sinbad_object)

