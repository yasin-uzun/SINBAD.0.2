
library(Seurat)
library(SINBAD)

source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/dimensionality_reduction.R')
source('/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/p99/MethylProc/R/Main.R')

sinbad_file = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/snmc/LEUK/batch_1/working_dir/877476/objects/sinbad_object.40.rds'

sinbad_object = readRDS(sinbad_file)


sinbad_object$name_for_dim_red = 'MAPLE'
sinbad_object$name_for_features = 'MAPLE'

seurat_object = process_with_Seurat(sinbad_object)

seurat_object  <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.1)

seurat_object_maple = seurat_object

DimPlot(seurat_object)

VlnPlot(seurat_object, features = 'RUNX1', pt.size = 0)

VlnPlot(seurat_object, features = 'CD19', pt.size = 0)

VlnPlot(seurat_object, features = 'CD22', pt.size = 0)

VlnPlot(seurat_object, features = 'HOXA9', pt.size = 0)

VlnPlot(seurat_object, features = 'MEF2C', pt.size = 0)



sinbad_object$name_for_dim_red = 'Promoters'
sinbad_object$name_for_features = 'Promoters'

sinbad_object_pro = wrap_read_annots(sinbad_object = sinbad_object)

sinbad_object_pro = wrap_quantify_regions(sinbad_object = sinbad_object_pro)

sinbad_object_pro = wrap_dim_red(sinbad_object = sinbad_object, annot_type = 'promoters')

seurat_object = process_with_Seurat(sinbad_object)

seurat_object  <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.1)

seurat_object_prom = seurat_object

DimPlot(seurat_object)

