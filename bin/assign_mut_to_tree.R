#!/usr/bin/env Rscript

options(stringsAsFactors = F)

# parsing arguments
{args = commandArgs(trailingOnly = TRUE)
  # Function to parse named arguments
  parse_arguments <- function(args) {
    parsed_args <- list()
    for (arg in args) {
      if (grepl("^--", arg)) {
        # Named argument found
        arg_split <- strsplit(arg, "=")[[1]]
        arg_name <- sub("^--", "", arg_split[1])
        arg_value <- arg_split[2]

        # Store the named argument in the list
        parsed_args[[arg_name]] <- arg_value
      }
    }
    return(parsed_args)
  }
  # Parse named arguments
  parsed_args <- parse_arguments(args)
  nr_path = parsed_args$nr_path
  nv_path = parsed_args$nv_path
  genotype_bin_path = parsed_args$genotype_bin_path
  tree_path = parsed_args$tree_path
  out_prefix = parsed_args$out_prefix
}

library(ape)
library(treemut)
library(ggplot2)
library(ggtree)


# load data
{
  NR_flt=read.table(nr_path, check.names = F)
  NV_flt=read.table(nv_path, check.names = F)
  genotype_bin = read.table(genotype_bin_path, check.names = F)
  tree=read.tree(tree_path)
  tree$edge.length=rep(1,nrow(tree$edge))
  tree=drop.tip(tree,"Ancestral")
  NR_flt = as.matrix(NR_flt[,tree$tip.label])
  NV_flt = as.matrix(NV_flt[,tree$tip.label])
}

# assigning mutations to tree
{
present_vars=rowSums(genotype_bin>0)>0
NR_flt=NR_flt[present_vars,]
NV_flt=NV_flt[present_vars,]
df=reconstruct_genotype_summary(tree)
res=assign_to_tree(tree = tree,
                   mtr=NV_flt,
                   dep=NR_flt)
filename = paste0(out_prefix, ".muts_assigned_to_tree.Rdata")
saveRDS(res,filename)
}

# getting the tree with branch length
{
tree_w_branchlength=tree
edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<0.01])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree_w_branchlength$edge.length=as.numeric(edge_length)
filename = paste0(out_prefix, ".tree_with_branch_length.tree")
write.tree(tree_w_branchlength, filename)

tree_w_branchlength_polytomised = as.polytomy(tree_w_branchlength, feature='branch.length', fun=function(x) as.numeric(x)==0)
filename = paste0(out_prefix, ".tree_with_branch_length_polytomised.tree")
write.tree(tree_w_branchlength_polytomised, filename)

p=ggtree(tree_w_branchlength_polytomised) +
  geom_tiplab(aes(x=branch),vjust=-0.3) +
  theme_tree2()+xlim(0,max(fortify(tree_w_branchlength_polytomised)$x)*1.3)
filename = paste0(out_prefix, ".tree_with_branch_length.pdf")
ggsave(filename = filename, plot = p)

tree_collapsed=tree_w_branchlength_polytomised
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
p = ggtree(tree_collapsed) + geom_tiplab(aes(x=branch),vjust=-0.3)+xlim(0,max(fortify(tree_collapsed)$x)*1.3)
filename = paste0(out_prefix, ".tree_with_equal_branch_length.pdf")
ggsave(filename = filename, plot = p)
}

# save VCF file to disk
{
Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_flt),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("CHROM","POS","REF","ALT")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<0.01,]
filename = paste0(out_prefix, ".muts_assigned_to_tree.txt")
write.table(Mutations_per_branch, filename, quote=F, row.names=F, sep="\t")
}

