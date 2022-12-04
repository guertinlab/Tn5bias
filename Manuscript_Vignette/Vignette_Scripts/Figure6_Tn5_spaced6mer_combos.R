source("https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_hcsearch.R")
#Input the span (46bp here) of the iterations
start_mask = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
#A function which iterates across the start_mask, creating all possible 5-mers
neighbors <- function(vecmask) {
  # generate list of neighboring masks that differ
  # by the addition of a single unmasked position
  result <- vector(mode="list", length=length(vecmask))
  n <- 0
  for (k in 1:length(vecmask)) {
    if (vecmask[k] == -1) { # X
      n <- n + 1
      vi = vecmask
      vi[k] = 1 # N
      vi[k+1] = 1 # 2nd N
      vi[k+2] = 1 # 3rd N
      result[[n]] <- vec.to.mask(vi)
    }
  }
  if (n == 0) return(NULL)
  #Remove last line so you don't get an extra N space (mappable bases = mer-1)
  zz <- (n-2)
  result[1:zz]
}
#Take the output from the first iteration of all possible 3-mers
initial_masks = unlist(neighbors(mask.to.vec(start_mask)))
#Make second iteration of 3-mers by making a permutation for each possible spacing
#of the original iteration across the 46-bp span
mask_store = vector(mode = 'list', length = length(initial_masks))
for (i in 1:length(initial_masks)) {
  mask_store[[i]] = unlist(neighbors(mask.to.vec(initial_masks[i])))
}
#Take each list object and turn it into a column of a dataframe
all_masks = data.frame(do.call(rbind, mask_store))
mask_col = NULL
#Take each column of the data frame and concatenate it
for (i in 1:ncol(all_masks)) {
  mask_col = c(mask_col, all_masks[,i])
}
#remove overlapping positions
mask_col = mask_col[-grep('NNNN' ,mask_col)]
#Find unique masks
mask_col = unique(mask_col)
write.table(mask_col, file = 'all_possible_spaced_6mers_46bp.txt', quote = FALSE,
            sep = '\n', row.names = FALSE)
