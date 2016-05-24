#########################################################
# Pearson's Chi-Squared Test of Independence: Used with #
# categorical data with two or more possible outcomes.  #
# Used to determine if the outcome is independt of the  #
# chosen variable.                                      #
#########################################################

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
    ## rid of any extra rows
    tempTable = tempTable[tempTable[,1] != 0,]
    
        
    ## Pearson's chi-squared test of independence with continutiy correction
    test = chisq.test(rmmeta, rmgroup)
    
    return (c(testMethods = gsub("\\'", "\\\\'", test$method),
              pvalues = test$p.value,
              types = 'categorical',
              labels = paste(tempRows[!tempRows %in% null_string], collapse = ','),
              gin = paste(tempTable[!tempRows %in% null_string,"IN"], collapse = ','),
              gout = paste(tempTable[!tempRows %in% null_string,"OUT"], collapse = ',')))
    
       
  }
}