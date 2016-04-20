# get arguments
# args[1] path for working directory 
# args[2] node, 
# args[3] genesetFile, 
# args[4] dendroFile, 
# args[5] backFile, 
# args[6] result file
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
    refseq[count] <<- node$name
    #print (node$name)
  }
}

data <- fromJSON(file=args[4], method='C')
count <- 0
refseq <- rep(NA,100)
findLeaves(findNode(data, args[2]))
refSeqIds <- na.omit(refseq)

if (!is.null(refSeqIds)) {
  # READS THE BACKGROUND FILE (ALL PROTEINS)
  NAMER=read.csv(args[5],sep=",",header=TRUE)
  
  temp = sub("(.*)\\.\\d+", "\\1", refSeqIds)
  if (sum(NAMER$id %in% temp)>0) {
    SOURCE = data.frame(id=temp)
  } else if (sum(NAMER$gene %in% temp)>0) {
    SOURCE = data.frame(gene=temp)
  } else {
    
  }
  
  SOURCE = merge(SOURCE, NAMER, all.x = TRUE)
  SOURCE = unique(SOURCE)
  TARGET = unique(SOURCE$gene[complete.cases(SOURCE$gene)])
  
  # READS THE DATABASE FILE
  load(args[3])
  
  # #########################################################
  # ### F-1) added preprocessing for gathering the number of genes commonly included in background and pathway (annotation).
  # #########################################################
  # USE THE GENE COLUMN AS INPUT
  BACK = unique(NAMER$gene)
  for (i in 1:length(k$names)) {
    thisname = k$names[i]
    thissize = k$sizes[i]
    thisdesc = k$desc[i]
    annot = k$matrix[i,1:thissize]
  
    back_annot = BACK[BACK %in% annot]
    thissize = k$sizes[i] = length(back_annot)
    if (thissize > 0) k$matrix[i,1:thissize] = paste(back_annot)
  }
  
  # #########################################################
  # ### G) Enrichment test
  # #########################################################
  enrichment_by_fishers <- function(group, background, annotation) {
    group_annot = sum(group %in% annotation)
    group_nonannot = length(group) - group_annot
    non_group_annot = length(annotation) - group_annot
    non_group_nonannot = length(background) - length(group) - non_group_annot
    
    test = matrix(c(group_annot, non_group_annot, group_nonannot, non_group_nonannot), nr=2,
            dimnames=list(c("Group", "NonGroup"), c("Annotated", "NonAnnotated")))
    
    per = c(test[1,1]/(test[1,1]+test[1,2]), test[2,1]/(test[2,1]+test[2,2]))
    
    ft = fisher.test(test, alternative = "greater")
    test = cbind(test, per)
    
    fold = per[1]/per[2]
    
    return(list(fisher=ft, mat=test, foldx=fold))
  }
  
  enrichment_in_groups <- function(geneset, targets, background, threshold=-1) {
    results = matrix(nrow=length(geneset$names), ncol=8)
    colnames(results) = c("geneset", "in_path", "out_path", "in_back", "out_back", "foldx", "pvalue", "bonferroni_pvalue")
    
    for (i in 1:length(geneset$names)) {
      thissize = geneset$sizes[i]
      if (thissize > 0) {
        thisname = geneset$names[i]
        thisdesc = geneset$desc[i]
        grouplist = geneset$matrix[i,1:thissize]
  
        # enr = enrichment_by_fishers(targets, background, grouplist, annot_size = thissize)
        enr = enrichment_by_fishers(targets, background, grouplist)
        p = enr$fisher$p.value
        f = enr$foldx
        mat = enr$mat
        
          if (mat[1,1] > threshold) {
            results[i,] = c(thisname, mat[1,1], mat[1,2], mat[2,1], mat[2,2], f, p, p*length(geneset$names))
          } else {
            
          }
      }
    }
    
    results[complete.cases(results), , drop=FALSE]
  }

  if(length(TARGET) > 0) {
    F=enrichment_in_groups(k, TARGET, BACK, threshold=0)
    if (dim(F)[1] == 0) {
      print("There is no pathway.")
      RESULT <- F
      new.p <- 1
    } else if (dim(F)[1] == 1) {
      print("There is one pathway.")
      RESULT <- F
      new.p <- as.numeric(RESULT[,'pvalue'])
    } else {
      print(paste0("There are ", dim(F)[1], " pathways."))
      RESULT <- F[order(as.numeric(F[,'pvalue'])),]
      # adjust the p-value for considering the multiple tests
      new.p <- p.adjust(as.numeric(RESULT[,'pvalue']), method ="BH")
    }

    RESULT <- cbind(RESULT, 'bh_pvalue'=new.p)
    idx <- which(k$names %in% RESULT[,'geneset'])
    desc <- data.frame('geneset'=k$names[idx],'desc'=k$desc[idx], 'link'=k$links[idx])
    RESULT <- merge(RESULT, desc, all.x = TRUE)
    write.table(RESULT, sep=",", file = paste0(args[6]), row.names=FALSE)
  }
}
