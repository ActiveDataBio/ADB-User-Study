# Title: Fisher's exact test
#
# Uses: Fisher's exact test is used when both variables being tested
# are categorical. Fisher's exact t-test assumes that each observation
# only fits into one cell and that the row and column totals are fixed
# and not random. It is used to test for independence between the two 
# variables, when sample sizes are small. The null hypothesis is that
# the proportions of one variable are independent of the other variable,
# while the alternative hypothesis is that the proportions of one
# variable are not independent of the other variable.
#
# Data format: character or numeric
#
# Author: Kaitlin Cornwell
#
# Date: June 1, 2016
#
# Notes: 


Snippet <- setRefClass("Snippet", contains = "Data", fields = "datatable",
                       methods = list(
                         ## Missing value handler function
                         cleaning = function(null_string) {
                           tryCatch({
                            ## initial read in
                            read_check(meta)
                            meta <<- as.character(meta)
                           
                            ## remove missing values
                            index = !(meta %in% null_string)
                            missing_check(index)
                            group <<- group[index]
                            meta <<- meta[index]
                            group_check(group)
                         },
                         
                         error = function(e) {
                           e$message = gsub("\n", " ", e$message)
                           errors <<- list(e$message, 2)
                         },
                         
                         finally = {
                           return(result(error = errors))
                         }) 
                        },
                        
                        assumptions = function() {
                          tryCatch ({
                            value_check(meta, group)
                            datatable <<- table(meta, group)
                            if (!is.null(freq_check(datatable, length(meta)))) {
                              errors <<- c(errors, 
                                           list("At least 80% of the expected counts are >5", 1))
                              }
                            },
                          
                          error = function(e) {
                            e$message = gsub("\n", " ", e$message)
                            errors <<- list(e$message, 2)
                          },
                          
                          finally = {
                            return(result(error = errors))
                          })
                        },
                        
                        test = function() {
                          tryCatch({
                            print(datatable)
                            test = fisher.test(meta, group, workspace = 1e+7)
                            rows = rownames(datatable)
                            return(result(test, c("column", "stacked-column", "percent-column"),
                                          rows[!rows %in% null_string],
                                          datatable[!rows %in% null_string, "IN"],
                                          datatable[!rows %in% null_string, "OUT"],
                                          errors))
                          },
                          
                          error = function(e) {
                            e$message = gsub("\n", " ", e$message)
                            return(result(error = list(e$message, 2)))
                          })
                        }
                      ))


freq_check <- function(datatable, length) {
  ## find row and column sums
  row_sums = vector(mode = "numeric", length = 0)
  for (i in 1:dim(datatable)[1]) {
    row_sums = c(row_sums, sum(datatable[i,]))
  }
  col_sums = vector(mode = "numeric", length = 0)
  for (j in 1:dim(datatable)[2]) {
    col_sums = c(col_sums, sum(datatable[,j]))
  }
  
  ## find expected frequencies to test assumptions and count number of 
  ## expected frequences <5
  count = 0
  
  for (i in 1:dim(datatable)[1]) {
    for (j in 1:dim(datatable)[2]) {
      if (((row_sums[i] * col_sums[j])/length) > 5)
        count = count + 1
    }
  }
  
  ## return percentage of expected frequences >5
  if ((count/(dim(datatable)[1] * dim(datatable)[2])) > .8)
    return("chisq")
  return(NULL)
}