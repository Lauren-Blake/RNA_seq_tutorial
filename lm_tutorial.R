# In R, the lm(), or “linear model,” function can be used to create a simple regression model.
# We might be interested in creating a linear model with one or more independent variables and gene expression as the dependent variable. Similar to limma, we will create one linear model for each gene. Unlike limma, it does not do information sharing across the genes (variance shrinkage, etc.)
# Basically, for one gene, we will use the format: lm(dependent_variable ~ independent_variable1 + independent_variable2 + etc.)
# To look at the summary statistics for lm, use the command summary around the lm command e.g. summary lm(dependent_variable ~ independent_variable1 + independent_variable2 + etc.)

# This tutorial uses data from the RNA seq tutorial that we did on Thursday. Therefore, before you run this, run lines 42-206 of the RNA-seq analysis with R/Bioconductor tutorial. 

# Pull the normalized gene expression for the LP and ML cells
lp_ml <- c(1,2,5,6,8,9)
lp_ml_exp <- v$E[,lp_ml]

group <- as.data.frame(group)
lp_ml_exp_group <- as.character(t(group[lp_ml,]))

# Let's test the first gene
lm(lp_ml_exp[1,] ~ lp_ml_exp_group)

# Let's get the summary statistics

summary(lm(lp_ml_exp[1,] ~ lp_ml_exp_group))

# How do we interpret the coefficient lp_ml_exp_groupML? 

# We want to test all of the genes instead of just 1. To do, so we 

# Make an array to store the p values from the lm
mypvals <- c()
for (i in 1:dim(lp_ml_exp)[1]) {
  
  # Run the lm for each gene
  data.sum <- summary(lm(lp_ml_exp[i,] ~ lp_ml_exp_group))
  
  # Pull the p values
  pval <- pf(data.sum$fstatistic[1], data.sum$fstatistic[2], data.sum$fstatistic[3], 
             lower.tail = FALSE)
  mypvals <- c(mypvals, pval)
}

# Look at the p-values
head(mypvals)
length(mypvals)

# How many genes have P value lower than the Bonferroni correction?

  # Set the threshold
threshold <- 0.05 / length(mypvals)
threshold

  # FALSE = number of non-DE genes
summary(mypvals < threshold)

# Here, there are no genes that meet genomewide significance. 
# This could be due to many reasons e.g. no correction for mean-variance relationship (voom weights), no incorporation of lane effects, no variance shrinkage/sharing of information across genes, the Bonferroni correction is a very stringent threshold, etc.