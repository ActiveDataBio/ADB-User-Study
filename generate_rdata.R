#########################################################
### A) Reading data and transform it into matrix format
#########################################################
data <- read.csv("dataset.csv", comment.char="#")
rnames1 <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames1                  # assign row names 

#########################################################
### B) reading meta data (group information)
#########################################################
metadata <- read.delim("metadata.tsv")
rnames2 <- metadata[3:nrow(metadata),1]                            # assign labels in column 1 to "rnames"
meta_data <- data.matrix(metadata[3:nrow(metadata),2:ncol(metadata)])  # transform column 2-5 into a matrix
rownames(meta_data) <- rnames2                  # assign row names 

# #########################################################
# ### C) generate cluster info for row
# #########################################################
d_row <- dist(mat_data)
d_row[is.na(d_row)] <- max(mat_data, na.rm = TRUE) ## adjust NA to maximum values
hc_row <- hclust(d_row)

# #########################################################
# ### D) generate cluster info for col
# #########################################################
if ('group' %in% colnames(meta_data)) {
  d_col <- dist(meta_data[,'group'])
} else {
  d_col <- dist(t(mat_data))
  d_col[is.na(d_col)] <- max(mat_data, na.rm = TRUE) ## adjust NA to maximum values
}
hc_col <- hclust(d_col)

# #########################################################
# ### E) create a list object
# #########################################################
hc.out <- list(rowDendrogram=hc_row, colDendrogram=hc_col, carpet=t(mat_data))
save(hc.out, file = 'data.RData')
