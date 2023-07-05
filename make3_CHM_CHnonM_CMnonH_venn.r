library(venn)
setwd('~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Overlap')

listInput_total <- list()
for(C_class in c('CHM', 'CHnonM', 'CMnonH')){
  earlyEmbryo <- as.character(read.table(paste0('EarlyEmbryogenesis.', C_class, '.txt'), sep = '\t', header = F)$V1)
  PGC <- as.character(read.table(paste0('PGCsDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
  spermatogenesis <- as.character(read.table(paste0('Spermatogenesis.', C_class, '.txt'), sep = '\t', header = F)$V1)
  retinal <- as.character(read.table(paste0('RetinalDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
  heart <- as.character(read.table(paste0('HeartDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
  liver <- as.character(read.table(paste0('LiverDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
  
  listInput <- list('Pre-implantation' = earlyEmbryo, 'PGC\ndevelopment' = PGC, Spermatogenesis = spermatogenesis, 
                    'Retinal\ndevelopment' = retinal, 'Heart\ndevelopment' = heart, 'Liver\ndevelopment' = liver)
  listInput_total[[C_class]] <- listInput
}


for(C_class in c('CHM', 'CHnonM', 'CMnonH')){
  pdf(file = paste0("~/CHMsInOtherContexts/figures/make3_venn6_R_", C_class, ".pdf"),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
  venn(
    listInput_total[[C_class]],
    zcolor = c('#9B1C3D', '#16557A', '#1B9E77', '#BF8B12', '#2E8BC0', '#666666'),
    col = c('#9B1C3D', '#16557A', '#1B9E77', '#BF8B12', '#2E8BC0', '#666666'),
    box = FALSE,
    ilcs = 0.8,
    sncs = 1
       )
  dev.off()
}
