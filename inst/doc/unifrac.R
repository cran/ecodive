## ----input_data, echo=FALSE---------------------------------------------------
library(ecodive)

counts <- matrix(
  data     = c(0, 0, 9, 3, 3, 1, 4, 2, 8, 0), 
  ncol     = 2, 
  dimnames = list(paste0('Species_', 1:5), c('Sample_A', 'Sample_B')) )

tree <- read_tree(
  underscores = TRUE,
  newick      = "
    (((Species_1:0.8,Species_2:0.5):0.4,Species_3:0.9):0.2,(Species_4:0.7,Species_5:0.3):0.6);" )

L <- tree$edge.length
A <- c(9,0,0,0,9,6,3,3)
B <- c(7,5,1,4,2,8,8,0)

# local({ # man/figures/unifrac-tree.png
#     
#   par(xpd = NA)
#   ape::plot.phylo(
#     x          = tree, 
#     direction  = 'downwards', 
#     srt        = 90, 
#     adj        = 0.5, 
#     no.margin  = TRUE,
#     underscore = TRUE,
#     x.lim      = c(0.5, 5.5) )
#   
#   ape::edgelabels(tree$edge.length, bg = 'white', frame = 'none', adj = -0.4)
# })
# 
# local({ # man/figures/unifrac-weights.png
#   
#   tree$edge.length <- c(1, 1, 1, 1, 2, 1, 2, 2)
#   
#   par(xpd = NA)
#   ape::plot.phylo(
#     x               = tree, 
#     direction       = 'downwards', 
#     srt             = 90, 
#     adj             = 0.5, 
#     no.margin       = TRUE,
#     underscore      = TRUE,
#     x.lim           = c(.8, 6) )
#   
#   ape::edgelabels(1:8, frame = 'circle')
#   
#   ape::edgelabels(paste('A =', A), bg = 'white', frame = 'none', adj = c(-0.4, -1.2))
#   ape::edgelabels(paste('B =', B), bg = 'white', frame = 'none', adj = c(-0.4,  0.0))
#   ape::edgelabels(paste('L =', L), bg = 'white', frame = 'none', adj = c(-0.3,  1.2))
# })

## ----input_data_tree, out.width = "100%", echo=FALSE--------------------------
  knitr::include_graphics('../man/figures/unifrac-tree.png')

## ----definitions, fig.align = 'center', out.width = "75%", echo=FALSE---------
knitr::include_graphics('../man/figures/unifrac-weights.png')

