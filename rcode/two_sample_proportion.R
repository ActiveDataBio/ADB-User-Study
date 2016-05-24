#########################################################
# Two sample proportions: Used with categorical data    #   
# with two possible outcomes. Used to determine whether #
# two sample proportions are equal.                     #
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
    
    #reduce contingency table to yes/no groups
    selected = c("NO", "YES")
    tempTable = tempTable[selected,]
    
  #test of equality between two sample propotions  
  test = prop.test(tempTable, alternative = "two.sided")
  
  return (c(testMethods = gsub("\\'", "\\\\'", test$method),
            pvalues = test$p.value,
            types = 'categorical',
            labels = paste(tempRows[!tempRows %in% null_string], collapse = ','),
            gin = paste(tempTable[!tempRows %in% null_string,"IN"], collapse = ','),
            gout = paste(tempTable[!tempRows %in% null_string,"OUT"], collapse = ',')))
  
  }
}
