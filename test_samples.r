#' to run statistical tests based on the user preference
#' @arguments a working directory of a source file
#' @arguments node name
#' @examples
#' Rscript stat_col.r "/path/dir" node10
#' get arguments
#' args[1] path for working directory 
#' args[2] node, 
#' args[3] metaFile, 
#' args[4] dendroFile, 
#' args[5] result file, 
#' args[6] custom codes
#
args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])
library(rjson)

findNode <- function(parent, name) {
  if (!is.null(parent$name) && parent$name == name) {
    return (parent)
  } else {
    if (!is.null(parent$children)) {
      len <- length(parent$children)
      for (i in 1:len) {
        node <- findNode(parent$children[[i]], name)
        if (!is.null(node)) return (node)
      }
    }
    else return (NULL);
  }
}

findLeaves <- function(node) {
  if (!is.null(node$children)) {
    len <- length(node$children)
    for (i in 1:len) {
      findLeaves(node$children[[i]])
    }
  } else {
    count <<- count + 1
    leaf[count] <<- node$name
    #print (node$name)
  }
}

#' customized parameters
temp <- fromJSON(args[6])
custom_codes <- NULL
if (length(temp) > 0) {
  for (i in 1:length(temp)) {
    new <- data.frame(col=temp[[i]]$names,snippet=c(rep(paste0(temp[[i]]$snippet$sha,'_',temp[[i]]$snippet$name),length(temp[[i]]$names))))
    if (is.null(custom_codes))  custom_codes<-new
    else custom_codes<-rbind(custom_codes,new)
  }  
}

## fetch a list of children of this node
nodeName <- args[2]
data <- fromJSON(file=args[4], method='C')
count <- 0
leaf <- rep(NA,100)
findLeaves(findNode(data, nodeName))
sampleIds <- na.omit(leaf)

## read the metadata file
metadata <- read.delim(args[3])
metaconfig <- metadata[grep('#', metadata$id),]
metadata <- metadata[-grep('#', metadata$id),]

#' fix factor
#' @param x an input array
fixFactor <- function(x) {
  if(is.factor(x)) factor(x) else x
}

################################################################
## find a column name for unique ids
################################################################
unique_id <- 'id'

## to identify null strings
null_string <- c("","[Not Applicable]","[Not Available]","[Pending]","[]")

target_group <- sampleIds[sampleIds != ""]

## check the idx for target group
match_idx <- pmatch(target_group, eval(parse(text=paste0('metadata$',unique_id))), dup = FALSE)  
match_idx <- match_idx[complete.cases(match_idx)]

## make a group according to the each nodes (subset)
#' @NOTE: be careful whether there is 'temp_group' column in meta data...
# metadata$temp_group = "OUT"
# metadata[match_idx,]$temp_group = 'IN'
group <- c()
for (i in (1:nrow(metadata))) {
  if (i %in% match_idx) group[i] <- 'IN'
  else group[i] <- 'OUT'
}

testMethods <- c()
pvalues <- c()
types <- c()
gin <- c()
gout <- c()
labels <- c()

for (i in 2:ncol(metadata)) {
  testMethods[i] <- NA
  pvalues[i] <- NA
  types[i] <- NA
  labels[i] <- NA
  gin[i] <- NA
  gout[i] <- NA

  col_name <- names(metadata[i])
  ## user preference for this column
  eval(parse(text=paste0('config=metaconfig$',col_name)))
  if (sum(config == 'no') == 1) next
  config <- fixFactor(config)
  
  eval(parse(text=paste0('meta=metadata$',col_name)))
  meta <- fixFactor(meta)
  if (!is.null(custom_codes) && col_name %in% custom_codes$col) {
    snippet <- as.character(custom_codes$snippet[custom_codes$col==col_name])
    print(paste0("custom code for ", col_name, ": ", snippet))
  } else {
    snippet <- paste0('../rcode/',config[1],'.r')
  }
  source(snippet)

  test_result = test(meta, group, null_string)
  if (is.null(test_result)) next
  testMethods[i] = test_result['testMethods']
  pvalues[i] = test_result['pvalues']
  labels[i] = test_result['labels']
  types[i] = test_result['types']
  gin[i] = test_result['gin']
  gout[i] = test_result['gout']
}
statTestResult = data.frame(names=names(metadata), types = types, methods = testMethods, pvalues = pvalues, labels=labels, group_in = gin, group_out = gout)

################################################################
## generate json string
################################################################
str = '{"stats":['
count = 1
for (j in (1:dim(statTestResult)[1])) {
  if (is.na(statTestResult$methods[j])) next
  if (statTestResult$types[j] == 'continuous') {
    temp = paste0('{"attribute":"',statTestResult$names[j],'", "datatype":"', statTestResult$types[j],'", "group_in":"', statTestResult$group_in[j],'", "group_out":"', statTestResult$group_out[j], '", "method":"', statTestResult$methods[j], '", "pvalue":', statTestResult$pvalues[j], '}')  
  } else if (statTestResult$types[j] == 'categorical') {
    temp = paste0('{"attribute":"',statTestResult$names[j],'", "datatype":"', statTestResult$types[j],'", "labels":"', statTestResult$labels[j],'", "group_in":"', statTestResult$group_in[j],'", "group_out":"', statTestResult$group_out[j], '", "method":"', statTestResult$methods[j], '", "pvalue":', statTestResult$pvalues[j], '}')  
  }
  
  if (count == 1) {
    str = paste0(str, temp)
  } else {
    str = paste0(str, ',', temp)
  }
  count = count + 1
}
str = paste0(str, ']', ',"node":"', nodeName, '"}')

fileConn<-file(paste0(args[5]))
writeLines(str, fileConn)
close(fileConn)
