# Title: Two sample Wilcoxon rank sum Test
#
# Uses: The two sample Wilcoxon rank sum test is used as an
# alternative to the two-sample t-test when the data is not
# normally distributed, the sample size is small, or when one 
# variable is ordinal. This test assumes that that one variable
# is either continuous or ordinal and one variable is categorical
# with two groups, that the groups within the categorical
# varaible are independent, and that observations are independent.
# The null hypothesis is that the groups are identically
# distributed while the alternative hypothesis is that the groups
# are not identically distributed. If the distributions of the
# two groups are not similar then the Wilcoxon rank sum test
# determines if the distributions are significantly different. If
# the two distributions are similar then the test determines if
# the medians of the two groups are significantly different.
#
# Data format: Values except for null/missing characters are numeric (may be
# encoded as strings).
# 
# Author: Kaitlin Cornwell
#
# Date: June 3, 2016
#
# Notes: The two-sample Wilcoxon rank sum test is also known as the Mann-
# Whitney U test.

test <- function(meta, group, null_string) {
  
  ret = tryCatch({
    if (is.null(meta)) {
      stop("Custom: null")
    }
    
    if (!is.character(meta)) {
      meta = as.character(meta)
    }
    
    ## Remove missing values
    mvidx = !(meta %in% null_string)
    if (length(which(mvidx)) == 0) {
      stop("Custom: null")
    }
    rmgroup = group[mvidx]
    rmmeta = meta[mvidx]
    rmmeta = as.numeric(rmmeta)
    
    ## Remove NAs
    mvidx = (!is.na(rmmeta))
    rmgroup = rmgroup[mvidx]
    rmmeta = rmmeta[mvidx]
    
    ## Check data
    if (length(rmmeta) == 0) {
      stop("Custom: type")
    }
    if (length(which(rmgroup %in% "IN")) == 0 ||
        length(which(rmgroup %in% "OUT")) == 0) {
      stop("Custom: nullgroup")
    }
    
    ## Check data type
    table = table(rmmeta)
    if (dim(table) == 1) {
      stop("Custom: oneval")
    }
    for (i in 1:dim(table)) {
      if (table[[i]] > (.2 * length(rmmeta))) {
        stop("Custom: type")
      }
    }
    
    ## Separate into "in" and "out" groups
    group_index = (rmgroup %in% "IN")
    meta_in = rmmeta[group_index]
    meta_out = rmmeta[!group_index]
    
    test = wilcox.test(meta_in, meta_out, alternative = "two.sided")
  },
  ## error handler function
  error = function (e) {
    if (grepl("Custom:", e)) {
      if(grepl("nullgroup", e)) {
        return(c("No data in a group", 2))
      }
      if (grepl("null", e)) {
        return(c("No meta data", 2))
      }
      if (grepl("oneval", e)) {
        return(c("Same value for each observation", 3))
      }
      if (grepl("type", e)) {
        return(c("Incorrect data type: received character instead of continuous or ordinal", 3))
      }
    }
    return(c(e, 1))
})
  
  ## if length of return statement is 2 then an error occurred
  ## return the error message and code
  if (length(ret) == 2) {
    return(c(msg = ret[1], status = ret[2]))
  }
  
  return(c(testMethods = gsub("\\'", "\\\\'", test$method),
            pvalues = test$p.value,
            charts = paste(c("box", "scatter"), collapse = ','),
            labels = '',
            gin = paste(meta_in, collapse = ','),
            gout = paste(meta_out, collapse = ','),
            msg = '',
            code = 0))
}