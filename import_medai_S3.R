library(aws.s3)

bucket <- "itx-acd-data-imm"
prefix <- "IBD/External/SHP647CD306/clinical/rawdata/shp647_306/"

arquivos <- get_bucket(bucket = bucket, prefix = prefix)

for (obj in arquivos) {
  save_object(
    object = obj$Key,
    bucket = bucket,
    file = file.path("/domino/datasets/local/clinical-trial-data/rawdata", basename(obj$Key))
  )
}