# get arguments (1: a working directory of a source file, 2: node name)
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
    sample[count] <<- node$name
    print (node$name)
  }
}

data <- fromJSON(file='dendro_col.json', method='C')
count <- 0
sample <- rep(NA,100)
findLeaves(findNode(data, args[2]))
sampleIds <- na.omit(sample)

nodeName <- args[2]

## read the metadata files (_data: data, _config: configurations for each column)
metadata <- read.delim("metadata.tsv")
metaconfig <- metadata[grep('#', metadata$id),]
metadata <- metadata[-grep('#', metadata$id),]

#' Run a statistical test whether two groups (with continous values) are significantly different or not.
#' 
#' @param group_in A set of values in a group.
#' @param group_out A set of values out of a group.
#' @param in_num The number of values in a group.
#' @param out_num The number of values out of a group.
#' @return A test object.
#' @examples
#' getTest(group_in=rnorm(30, mean=0, sd=3), group_out=rnorm(40, mean=5, sd=10), in_num=30, out_num=40)
getTest<-function (group_in, group_out, in_num, out_num) {
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
##############
fixFactor <- function(x) {
  if(is.factor(x)) factor(x) else x
}
##############
################################################################
## find a column name for unique ids
################################################################
unique_id <- 'id'
# for (i in 1:ncol(metadata)) {
#   eval(parse(text=paste0('config=metaconfig$',names(metadata[i]))))
#   if (sum(config == 'unique') == 1) {
#     unique_id <- names(metadata[i])
#     break;
#   }
# }

# null string
null_string <- c("","[Not Applicable]","[Not Available]","[Pending]","[]")

target_group <- sampleIds[sampleIds != ""]
## check the idx for target group
if (grepl("TCGA",target_group[1])) {
  match_idx <- pmatch(gsub("(\\w*)\\.(\\w*)\\.(\\d*)\\.(\\d*)", "\\2\\-\\3\\-\\4", target_group), eval(parse(text=paste0('metadata$',unique_id))), dup = FALSE)  
} else {
  match_idx <- pmatch(target_group, eval(parse(text=paste0('metadata$',unique_id))), dup = FALSE)  
}

match_idx <- match_idx[complete.cases(match_idx)]
# make a group according to the each nodes (subset)
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
  
  eval(parse(text=paste0('config=metaconfig$',names(metadata[i]))))
  # if (sum(config == 'no') == 1 || sum(config == 'unique') == 1) next
  if (sum(config == 'no') == 1) next
  config <- fixFactor(config)
  
  
  eval(parse(text=paste0('meta=metadata$',names(metadata[i]))))
  meta <- fixFactor(meta)
  # multivariate categorical data
  if (sum(config == 'categorical') == 1) {
    # remove missing values
    mvidx = complete.cases(meta)& !(meta %in% null_string)
    rmgroup = group[mvidx]
    rmmeta = meta[mvidx]
    level_size = dim(table(unclass(rmmeta)))
    
    if (level_size == 1 || level_size == nrow(metadata)) {
      #print(paste0(i, ',', colnames(metadata)[i], ', categorical, N/A'))
    } else {
      tempTable = table(rmmeta, rmgroup)
      if (dim(tempTable)[2] != 2) next
      tempRows = rownames(tempTable)
      labels[i] = paste(tempRows[!tempRows %in% null_string], collapse = ',')
      gin[i] = paste(tempTable[!tempRows %in% null_string,"IN"], collapse = ',')
      gout[i] = paste(tempTable[!tempRows %in% null_string,"OUT"], collapse = ',')
      
      # chi sq test
      test = chisq.test(rmmeta, rmgroup)
      # fisher
      #test = fisher.test(rmmeta, rmgroup)
      #print(paste0(i, ',', colnames(metadata)[i], ', categorical, ', test$method, ', ', test$p.value))
      testMethods[i] = gsub("\\'", "\\\\'", test$method)
      pvalues[i] = test$p.value
      types[i] = 'categorical'
    }
  } else if (sum(config == 'continuous') == 1) {
    # divide into two group
    group_in = meta[(1:length(meta)) %in% match_idx]
    group_out = meta[!(1:length(meta)) %in% match_idx]
    # the number of missing values
    in_idx = complete.cases(group_in)& !(group_in %in% null_string)
    in_num = sum(in_idx)
    if(in_num < 1) next
    out_idx = complete.cases(group_out)& !(group_out %in% null_string)
    out_num = sum(complete.cases(group_out))
    if(out_num < 1) next
    # remove missing values
    rmgroup_in = group_in[in_idx]
    rmgroup_out = group_out[out_idx]
    
    if (!is.numeric(rmgroup_in)) rmgroup_in = as.numeric(as.character(rmgroup_in))
    if (!is.numeric(rmgroup_out)) rmgroup_out = as.numeric(as.character(rmgroup_out))
    if (length(rmgroup_in) == 0 || length(rmgroup_out) == 0) next
    test = getTest(group_in=rmgroup_in, group_out=rmgroup_out, in_num=in_num, out_num=out_num)
    #print(paste0(i, ',', colnames(metadata)[i], ', continuous, ', test$method, ', ', test$p.value))
    testMethods[i] = gsub("\\'", "\\\\'", test$method)
    pvalues[i] = test$p.value
    types[i] = 'continuous'
    gin[i] = paste(rmgroup_in, collapse = ',')
    gout[i] = paste(rmgroup_out, collapse = ',')
  }
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

fileConn<-file(paste0("stat_",nodeName,".json"))
writeLines(str, fileConn)
close(fileConn)