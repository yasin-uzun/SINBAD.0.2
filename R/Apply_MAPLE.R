#library(devtools)
#install_github("yasin-uzun/SINBAD.1.0")

#sample_name = '877476'
#sample_name = '878516'
sample_name = '882304'
working_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples//data/snmc//working_dir/LEUK/'

#working_dir = args[1]
#sample_name = args[2]
#read_type = args[3]

#working_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples//data/snmc//working_dir/LEUK/'
library(SINBAD)
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/Main.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/preprocessing.R')

source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/alignment.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/methylation_calling.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/regional_quantification.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/dimensionality_reduction.R')

raw_fastq_dir = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/RawData/FirstBatch_0730/'

demux_index_file = paste0('/mnt/isilon/tan_lab/uzuny/projects/sinbad/data/example/', '/demux_index.txt')

library(doSNOW)


#sample_name <- 'Leuk_Batch_1'
sinbad_object = construct_sinbad_object(working_dir = working_dir,
                                        raw_fastq_dir = raw_fastq_dir,
                                        demux_index_file = demux_index_file,
                                        sample_name = sample_name)


library(MAPLE)
library(GenomicRanges)

#source('/mnt/isilon/tan_lab/uzuny/projects/MethylPredict/software_package/R_package/v15/MethylPredict/R/met_matrix_processing.R')
#source('/mnt/isilon/tan_lab/uzuny/projects/MethylPredict/software_package/R_package/v15/MethylPredict/R/set_ops.R')

intersect_bed

#Set directory names
annot_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/snmc/maple/annot/'
model_dir = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/snmc/maple/models/Clark/'
cov_dir = paste0(working_dir, sample_name, '/methylation_calls/')


#Set input files
annot_file = paste0(annot_dir,'/gencode.hg38.v29.genes.bed')
cpg_content_file =  paste0(annot_dir,'/regions.genes.tss_ud_5K.cpg_ratio.bin_size_500.hg38.rds')

#Compute binned data
binned_list = compute_binned_met_counts(cov_dir = cov_dir, annot_file = annot_file )


met_all = binned_list[["df_binned_met"]]
demet_all = binned_list[["df_binned_demet"]]

met_all$cell = gsub('.organism', '', met_all$cell)
demet_all$cell = gsub('.organism', '', demet_all$cell)

rownames(met_all) = gsub('.organism', '', rownames(met_all) )
rownames(demet_all) = gsub('.organism', '', rownames(demet_all) )

met_all[1:5, 1:5]
demet_all[1:5, 1:5]

sinbad_object_file = paste0(working_dir, sample_name, '/objects/sinbad_object.38.5.rds')
sinbad_object = readRDS(sinbad_object_file)
head(sinbad_object$df_alignment_stats)


failed_cells = rownames(sinbad_object$df_alignment_stats[sinbad_object$df_alignment_stats$pass == 0, ])
failed_cells = sub('_', '___', failed_cells)
binned_list[["df_binned_met"]] = met_all[! met_all$cell %in% failed_cells, ]
dim(binned_list[["df_binned_met"]] )

binned_list[["df_binned_demet"]] = demet_all[! demet_all$cell %in% failed_cells, ]



#Compute meta cells
meta_object = compute_meta_cells(df_met =  binned_list[["df_binned_met"]],
                                 df_demet =  binned_list[["df_binned_demet"]])

#Generate features
fr_list = get_fr_list(meta_data = meta_object, cpg_content_file = cpg_content_file)

#Load CNN model and predict
cnn_model_file = paste0(model_dir, '/cnn_model.hd5')
predict_cnn = cnn_predict(fr_list, cnn_model_file)

#Load Elastic model and predict
elastic_model_file = paste0(model_dir, '/elastic_model.rds')
predict_elastic = elastic_predict(fr_list, elastic_model_file)

#Load RF model and predict
rf_model_file = paste0(model_dir, '/rf_model.rds')
predict_rf = rf_predict(fr_list, rf_model_file)

#Compute Ensemble prediction
prediction_list = list(predict_cnn, predict_elastic, predict_rf)
predict_ensem = ensemble_predict(prediction_list)

print(head(predict_ensem))

#Convert prediction into matrix format (genesxcells)
gene_activity_matrix = convert_preds_to_matrix(predict_ensem)
gene_activity_matrix[1:5, 1:5]

gene_activity_matrix_file = paste0(working_dir, sample_name,  '/matrices/maple_matrix.rds')
saveRDS(gene_activity_matrix, file = gene_activity_matrix_file)


sinbad_object$met_mat_for_dim_red = gene_activity_matrix
sinbad_object$met_mat_for_features = gene_activity_matrix

sinbad_object$name_for_dim_red = 'MAPLE'
sinbad_object$name_for_features = 'MAPLE'


sinbad_object = wrap_dim_red(sinbad_object)




sinbad_object_file = paste0(working_dir, sample_name, '/objects/sinbad_object.40.rds')
#saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
gene_activity_matrix = sinbad_object$met_mat_for_dim_red
colnames(gene_activity_matrix) = gsub('.organism', '', colnames(gene_activity_matrix))
gene_activity_matrix[1:5, 1:5]


gene_activity_matrix_file = paste0(working_dir, sample_name, '/matrices/maple_matrix.rds')
saveRDS(gene_activity_matrix, file = gene_activity_matrix_file)



