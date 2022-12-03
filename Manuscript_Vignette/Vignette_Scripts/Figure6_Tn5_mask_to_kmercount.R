library('doParallel')
masks = list.files('./')
masks = masks[grep('.tbl', masks)]
#
registerDoParallel(1)
foreach (i = 1:length(masks)) %dopar% {system(paste('seqOutBias table ',
        masks[i], ' ../C1_gDNA_rep1.bam > ', substr(masks[i], 1, nchar(masks[i])-4),
        '_fulldat_scale_factors.txt', sep = ''))}