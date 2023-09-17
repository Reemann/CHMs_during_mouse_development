#!/usr/bin/Rscript
setwd('~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Overlap')
library('UpSetR')
require(lattice)

listInput_total <- list()
for(C_class in c('CHM', 'CHnonM', 'CMnonH')){
	earlyEmbryo <- as.character(read.table(paste0('EarlyEmbryogenesis.', C_class, '.txt'), sep = '\t', header = F)$V1)
	PGC <- as.character(read.table(paste0('PGCsDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
	spermatogenesis <- as.character(read.table(paste0('Spermatogenesis.', C_class, '.txt'), sep = '\t', header = F)$V1)
	retinal <- as.character(read.table(paste0('RetinalDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
	heart <- as.character(read.table(paste0('HeartDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
	liver <- as.character(read.table(paste0('LiverDevelopment.', C_class, '.txt'), sep = '\t', header = F)$V1)
	
	listInput <- list(earlyEmbryo = earlyEmbryo, PGC = PGC, spermatogenesis = spermatogenesis, retinal = retinal, heart = heart, liver = liver)
  listInput_total[[C_class]] <- listInput
}



pdf(paste0('~/CHMsInOtherContexts/figures/make6_Overlap_UpSetR_CHM.pdf'), width = 8, height = 6, useDingbats = F, onefile=F)
p <- upset(fromList(listInput_total[['CHM']]), nsets = 6, point.size = 2, line.size = 0.75,
           mb.ratio = c(0.75, 0.25),
           order.by = 'freq',
           mainbar.y.label = paste0("Number of stable CHMs intersections"),
           sets.x.label = paste0("Number of stable CHMs"))
p
dev.off()

pdf(paste0('~/CHMsInOtherContexts/figures/make6_Overlap_UpSetR_CHnonM.pdf'), width = 8, height = 6, useDingbats = F, onefile=F)
p <- upset(fromList(listInput_total[['CHnonM']]), nsets = 6, point.size = 2, line.size = 0.75,
           mb.ratio = c(0.75, 0.25),
           order.by = 'freq',
           mainbar.y.label = paste0("Number of stable CH-non-Ms intersections"),
           sets.x.label = paste0("Number of stable CH-non-Ms"))
p
dev.off()

pdf(paste0('~/CHMsInOtherContexts/figures/make6_Overlap_UpSetR_CMnonH.pdf'), width = 8, height = 6, useDingbats = F, onefile=F)
p <- upset(fromList(listInput_total[['CMnonH']]), nsets = 6, point.size = 2, line.size = 0.75,
           mb.ratio = c(0.75, 0.25),
           order.by = 'freq',
           mainbar.y.label = paste0("Number of stable CM-non-Hs intersections"),
           sets.x.label = paste0("Number of stable CM-non-Hs"))
p
dev.off()


# Ref: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html