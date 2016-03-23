test <- function (meta, group, null_string) {
  ## remove missing values
  mvidx = complete.cases(meta) & !(meta %in% null_string)
  rmgroup = group[mvidx]
  rmmeta = meta[mvidx]
  # divide into two group
  rmgroup_in = rmmeta[which (rmgroup=='IN')]
  rmgroup_out = rmmeta[which (rmgroup=='OUT')]
  # the number of missing values
  if(length(rmgroup_in) == 0) return (NULL)
  if(length(rmgroup_out) == 0) return (NULL)
  
  if (!is.numeric(rmgroup_in)) rmgroup_in = as.numeric(as.character(rmgroup_in))
  if (!is.numeric(rmgroup_out)) rmgroup_out = as.numeric(as.character(rmgroup_out))
  if (length(rmgroup_in) == 0 || length(rmgroup_out) == 0) return (NULL)
  test = getTest(group_in=rmgroup_in, group_out=rmgroup_out)
  
  return (c(testMethods = gsub("\\'", "\\\\'", test$method),
            pvalues = test$p.value,
            types = 'continuous',
            labels = '',
            gin = paste(rmgroup_in, collapse = ','),
            gout = paste(rmgroup_out, collapse = ',')))
}

#' Run a statistical test whether two groups (with continous values) are significantly different or not.
#' 
#' @param group_in A set of values in a group.
#' @param group_out A set of values out of a group.
#' @return A test object.
#' @examples
#' getTest(group_in=rnorm(30, mean=0, sd=3), group_out=rnorm(40, mean=5, sd=10), in_num=30, out_num=40)
getTest<-function (group_in, group_out) {
  in_num = length(group_in)
  out_num = length(group_out)
  if (in_num < 10 || out_num < 10) {
    test = wilcox.test(group_in, group_out)
  } else if (in_num > 30 && out_num > 30) {
    test = var.test(group_in, group_out)
    if (test$p.value < 0.05) {
      test = t.test(group_in, group_out, var.equal=FALSE)
    } else {
      test = t.test(group_in, group_out, var.equal=TRUE)
    }
  } else {
    normality = TRUE
    if (in_num < 30) {
      in_norm_test = shapiro.test(group_in)
      normality = in_norm_test$p.value >= 0.1
    }
    
    if(normality && out_num < 30) {
      out_norm_test = shapiro.test(group_out)
      normality = out_norm_test$p.value >= 0.1
    }
    
    if(normality) {
      test = var.test(group_in, group_out)
      if (test$p.value < 0.05) {
        test = t.test(group_in, group_out, var.equal=FALSE)
      } else {
        test = t.test(group_in, group_out, var.equal=TRUE)
      }
    } else {
      test = wilcox.test(group_in, group_out)
    }
  }
  return(test)
}