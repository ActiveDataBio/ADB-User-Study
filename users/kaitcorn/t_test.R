# Title: Two-sample t-test
#
# Uses: The two-sample t-test is used when one variable is continuous and one
# variable is categorical with two groups. This test assumes that samples have
# been drawn from normally distributed popluations and that the popluations have
# equal variances. The two-sample t-test is used to test for differences in mean
# between the groups defined by the categorical variable. The null hypothesis is
# that the means are equal, while the alternative hypothesis is that the means
# are not equal. Since this test uses the mean, the two-sample t-test should not 
# be chosen when the data has a large range or is signifcanly skewed.
#
# Data format: numeric
#
# Author: Kaitlin Cornwell
# 
# Date: June 3, 2016
#
# Notes: 

error = function(meta, group, null_string) {
  ret = tryCatch({
    ## check data input
    if (is.null(meta)) {
      stop("Custom: null")
    }
    if(!is.character(meta)) {
      meta = as.character(meta)
    }
    
    ## remove missing values
    mvidx = !(meta %in% null_string)
    if (length(which(mvidx)) == 0) {
      stop("Custom: null")
    }
    rmgroup = group[mvidx]
    rmmeta = meta[mvidx]
    
    ## Change data to numeric form
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
    
    ## Split data into "IN" and "OUT" groups
    group_index = (rmgroup %in% "IN")
    meta_in = rmmeta[group_index]
    meta_out = rmmeta[!group_index]
    
    ## Check normality assumptions
    if (shapiro.test(meta_in)$p.value < .05 ||
        shapiro.test(meta_out)$p.value < .05) {
      stop("Custom: normality")
    }
    
    ## Test varaince and perform t-test
    var = var.test(meta_in, meta_out, alternative = "two.sided")
    if (var$p.value < .05) {
      t.test(meta_in, meta_out,
                    alternative = "two.sided", var.equal = FALSE)
    } else {
      t.test(meta_in, meta_out,
                    alternative = "two.sided", var.equal = TRUE)
    }
  }, 
  ## error handler function
  error = function (e) {
    if (grepl("Custom:", e)) {
      if (grepl("nullgroup", e)) {
        return(list("No data in a group", 2))
      }
      if (grepl("null", e)) {
        return(list("No meta data", 2))
      }
      if (grepl("oneval", e)) {
        return(list("Same value for each observation", 3))
      }
      if (grepl("type", e)) {
        return(list("Incorrect data type: received categorical instead of continuous", 3))
      }
      if (grepl("normality", e)) {
        return(list("Data is not normally distributed, use Wilcoxon test instead", 3))
      }
    }
    return(list(e$message, 1))
  })
  
  ## if length of return statement is 2 then an error occurred
  ## return the error message and code
  if (length(ret) == 2) {
    return(list(msg = ret[[1]], status = ret[[2]]))
  }
  
  return(list(
           msg = '',
           status = 0))
}

test = function(meta, group, null_string) {
  meta = as.character(meta)
  
  mvidx = !(meta %in% null_string)
  rmgroup = group[mvidx]
  rmmeta = meta[mvidx]
  rmmeta = as.numeric(rmmeta)
  
  ## Remove NAs
  mvidx = (!is.na(rmmeta))
  rmgroup = rmgroup[mvidx]
  rmmeta = rmmeta[mvidx]
  
  group_index = (rmgroup %in% "IN")
  meta_in = rmmeta[group_index]
  meta_out = rmmeta[!group_index]
  ## Test varaince and perform t-test
  var = var.test(meta_in, meta_out, alternative = "two.sided")
  if (var$p.value < .05) {
    test = t.test(meta_in, meta_out,
           alternative = "two.sided", var.equal = FALSE)
  } else {
    test = t.test(meta_in, meta_out,
           alternative = "two.sided", var.equal = TRUE)
  }
  
  return(list(method = gsub("\\'", "\\\\'", test$method),
              pvalue = test$p.value,
              charts = c('box','scatter'),
              labels = '',
              group_in = paste(meta_in, collapse = ','),
              group_out = paste(meta_out, collapse = ',')))
}