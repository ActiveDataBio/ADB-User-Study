test <- function (meta, group, null_string) {
  ## remove missing values
  mvidx = complete.cases(meta) & !(meta %in% null_string)
  rmgroup = group[mvidx]
  rmmeta = meta[mvidx]
  level_size = dim(table(unclass(rmmeta)))
  
  if (level_size == 1 || level_size == nrow(metadata)) {
    return (NULL) 
  } else {
    ## make a contingency table
    tempTable = table(rmmeta, rmgroup)
    if (dim(tempTable)[2] != 2) next
    tempRows = rownames(tempTable)
#     labels = paste(tempRows[!tempRows %in% null_string], collapse = ',')
#     gin[i] = paste(tempTable[!tempRows %in% null_string,"IN"], collapse = ',')
#     gout[i] = paste(tempTable[!tempRows %in% null_string,"OUT"], collapse = ',')
    
    # chi sq test
    # test = chisq.test(rmmeta, rmgroup)
    # fisher
    test = fisher.test(rmmeta, rmgroup)
#     testMethods[i] = gsub("\\'", "\\\\'", test$method)
#     pvalues[i] = test$p.value
#     types[i] = 'categorical'
    
    return (c(testMethods = gsub("\\'", "\\\\'", test$method),
              pvalues = test$p.value,
              types = 'categorical',
              labels = paste(tempRows[!tempRows %in% null_string], collapse = ','),
              gin = paste(tempTable[!tempRows %in% null_string,"IN"], collapse = ','),
              gout = paste(tempTable[!tempRows %in% null_string,"OUT"], collapse = ',')))
  }
}