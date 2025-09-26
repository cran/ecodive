# Copyright (c) 2025 ecodive authors
# Licensed under the MIT License: https://opensource.org/license/mit


#' documentation
#' 
#' @name documentation
#' @keywords internal
#' 
#' @param alpha   How much weight to give to relative abundances; a value 
#'        between 0 and 1, inclusive. Setting `alpha=1` is equivalent to 
#'        `normalized_unifrac()`.
#' 
#' @param counts   A numeric matrix of count data where each column is a 
#'        feature, and each row is a sample. Any object coercible with 
#'        `as.matrix()` can be given here, as well as `phyloseq`, `rbiom`, 
#'        `SummarizedExperiment`, and `TreeSummarizedExperiment` objects. 
#' 
#' @param cpus   How many parallel processing threads should be used. The
#'        default, `n_cpus()`, will use all logical CPU cores.
#' 
#' @param cutoff   The maximum number of observations to consider "rare".
#'        Default: `10`.
#' 
#' @param digits   Precision of the returned values, in number of decimal 
#'        places. E.g. the default `digits=3` could return `6.392`.
#' 
#' @param pairs   Which combinations of samples should distances be 
#'        calculated for? The default value (`NULL`) calculates all-vs-all. 
#'        Provide a numeric or logical vector specifying positions in the 
#'        distance matrix to calculate. See examples.
#' 
#' @param power   Scaling factor for the magnitude of differences between
#'        communities (\eqn{p}). Default: `1.5`
#' 
#' @param pseudocount   The value to add to all counts in `counts` to prevent 
#'        taking `log(0)` for unobserved features. The default, `NULL`, selects 
#'        the smallest non-zero value in `counts`.
#' 
#' @param norm   Normalize the incoming counts. Options are:
#'        \describe{
#'            \item{`norm = "percent"` - }{ Relative abundance (sample abundances sum to 1). }
#'            \item{`norm = "binary"`  - }{ Unweighted presence/absence (each count is either 0 or 1). }
#'            \item{`norm = "clr"`     - }{ Centered log ratio. }
#'            \item{`norm = "none"`    - }{ No transformation. }
#'        }
#'        Default: `'percent'`, which is the expected input for these formulas.
#' 
#' @param tree   A `phylo`-class object representing the phylogenetic tree for 
#'        the OTUs in `counts`. The OTU identifiers given by `colnames(counts)` 
#'        must be present in `tree`. Can be omitted if a tree is embedded with
#'        the `counts` object or as `attr(counts, 'tree')`.
#' 
NULL




#### alpha_div refs -----

# #' @references
# #' 
# #' Alroy, J. (2018). Limits to Species Richness in Terrestrial Communities.
# #' *Ecology Letters*, 21(12): 1781-1789. \doi{10.1111/ele.13152}
# #' 
# #' Chao A 1984.
# #' Non-parametric estimation of the number of classes in a population.
# #' Scandinavian Journal of Statistics, 11:265-270.
# #' 
# #' Chao, M. T., & Lee, S. M. (1992). Nonparametric estimation of the number of
# #' species in a given collection of species.
# #' *Journal of the American Statistical Association*, 87(419), 506-512.
# #' \doi{10.2307/2290711}
# #' 
# #' Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic foraminifera
# #' in deep-sea sediments. *Science*, 168(3937), 1345–1347.
# #' \doi{10.1126/science.168.3937.1345}
# #' 
# #' Brillouin, L. (1956). Science and information theory. Academic Press.
# #' 
# #' Fisher, R.A., Corbet, A.S., & Williams, C.B. (1943). The relation between the
# #' number of species and the number of individuals in a random sample of animal
# #' population. Journal of Animal Ecology, 12, 42-58.
# #' 
# #' McIntosh, R.P. (1967). An Index of Diversity and the Relation of Certain
# #' Concepts to Diversity. *Ecology*, 48(3), 392–404. \doi{10.2307/1932674}
# #' 
# #' Shannon CE, Weaver W 1949.
# #' The Mathematical Theory of Communication.
# #' University of Illinois Press.
# #' 
# #' Simpson EH 1949.
# #' Measurement of diversity.
# #' Nature, 163.
# #' \doi{10.1038/163688a0}
# #' 
# #' Faith DP (1992).
# #' Conservation evaluation and phylogenetic diversity.
# #' *Biological Conservation*, 61:1-10.
# #' \doi{10.1016/0006-3207(92)91201-3}
# #' 
# #' Gamito, S. (2010). "Caution is needed when applying Margalef diversity index".
# #' *Ecological Indicators*, 10(2), 550-551. \doi{10.1016/j.ecolind.2009.07.006}
# #' 
# #' Margalef, R. (1958). "Temporal succession of planktonic communities in a limited area".
# #' 
# #' Menhinick, E. F. (1964). A comparison of some species-individuals diversity
# #' indices applied to samples of field insects. *Ecology*, 45(4), 859-861.
# #' \doi{10.2307/1934933}
# #' 
# #' 
# #' 
# #' 



#### beta_div refs -----

# #' @references
# #' 
# #' 
# #' Aitchison, J. (1982), The Statistical Analysis of Compositional Data. *Journal
# #' of the Royal Statistical Society: Series B (Methodological)*, 44, 139-160.
# #' \doi{10.1111/j.2517-6161.1982.tb01195.x}
# #' 
# #' Aitchison, J. (1986) The Statistical Analysis of Compositional Data. Chapman
# #' and Hall, London, 416 p. \doi{10.1007/978-94-009-4109-0}
# #' 
# #' Aitchison, J. (1989) Measures of Location of Compositional Data Sets.
# #' *Mathematical Geology*, 21, 787-790. \doi{10.1007/BF00893322}
# #' 
# #' Aitchison, J (1992). On criteria for measures of compositional difference.
# #' *Mathematical Geology*, 24, 365–379. \doi{10.1007/BF00891269}
# #' 
# #' Aitchison, J., Barcelo-Vidal, C., Martín-Fernandez, J.A. et al (2000). Logratio
# #' Analysis and Compositional Distance. *Mathematical Geology*, 32, 271–275.
# #' \doi{10.1023/A:1007529726302}
# #' 
# #' Lozupone C, Knight R (2005). UniFrac: A new phylogenetic method for comparing
# #' microbial communities.  *Applied and Environmental Microbiology*, 71(12).
# #' \doi{10.1128/AEM.71.12.8228-8235.2005}
# #' 
# #' Lozupone CA, Hamady M, Kelley ST, Knight R (2007). Quantitative and
# #' Qualitative \eqn{\beta} Diversity Measures Lead to Different Insights into
# #' Factors That Structure Microbial Communities.
# #' *Applied and Environmental Microbiology*, 73(5). \doi{10.1128/AEM.01996-06}
# #' 
# #' Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD, Collman RG,
# #' Bushman FD, Li H (2012). Associating microbiome composition with
# #' environmental covariates using generalized UniFrac distances.
# #' *Bioinformatics*, 28(16). \doi{10.1093/bioinformatics/bts342}
# #' 
# #' Chang Q, Luan Y, Sun F (2011). Variance adjusted weighted UniFrac: a powerful
# #' beta diversity measure for comparing communities based on phylogeny.
# #' *BMC Bioinformatics*, 12. \doi{10.1186/1471-2105-12-118}
# #' 
# #' #' Bhattacharyya, A. (1943). On a measure of divergence between two statistical
# #' populations defined by their probability distributions. *Bulletin of the
# #' Calcutta Mathematical Society*, 35, 99–109.
# #' 
# #' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures
# #' between probability density functions. *International Journal of Mathematical
# #' Models and Methods in Applied Sciences*, 1(4), 300–307.
# #' 
# #' Clark, P. J., & Evans, F. C. (1954). Distance to nearest neighbor as a
# #' measure of spatial relationships in populations. *Ecology*, 35(4), 445–453.
# #' \doi{10.2307/1931034}
# #' 
# #' Hellinger, E. (1909). Neue Begründung der Theorie quadratischer Formen von
# #' unendlichvielen Veränderlichen. *Journal für die reine und angewandte
# #' Mathematik*, 1909(136), 210–271. \doi{10.1515/crll.1909.136.210}
# #' 
# #' Minkowski, H. (1896). *Geometrie der Zahlen*. Teubner,
# #' Leipzig.
# #' 
# #' Gower JC (1971). A general coefficient of similarity and some of its
# #' properties. *Biometrics*, 27(4). \doi{10.2307/2528823}
# #' 
# #' Gower JC, Legendre P (1986). Metric and Euclidean Properties of Dissimilarity
# #' Coefficients. *Journal of Classification*, 3. \doi{10.1007/BF01896809}
# #' 
# #' Legendre P, Caceres M (2013). Beta diversity as the variance of community
# #' data: dissimilarity coefficients and partitioning. *Ecology Letters*, 16(8).
# #' \doi{10.1111/ele.12141}
# #' 
# #' Levy, A., Shalom, B. R., & Chalamish, M. (2024). A guide to similarity
# #' measures. *arXiv*. \doi{10.48550/arXiv.2408.07706v1}
# #' 
# #' Cha, S.-H. (2007). Comprehensive survey on distance/similarity measures
# #' between probability density functions. *International Journal of Mathematical
# #' Models and Methods in Applied Sciences*, 1(4), 300–307.
# #' 
# #' Bray JR and Curtis JT (1957). An ordination of the upland forest communities
# #' of southern Wisconsin. *Ecological Monographs*, 27(4). \doi{10.2307/1942268}
# #' 
# #' Lance GN and Williams WT (1967). A general theory of classificatory sorting
# #' strategies II. Clustering systems. *The computer journal*, 10(3).
# #' \doi{10.1093/comjnl/10.3.271}
# #' 
# #' Horn, H. S. (1966). Measurement of "Overlap" in comparative ecological
# #' studies. The American Naturalist 100:419-424.
# #' 
# #' Morisita, M. (1959). "Measuring of the dispersion and analysis of
# #' distribution patterns". *Memoires of the Faculty of Science*, Kyushu
# #' University, Series E. Biology. 2: 215–235.
# #' 
# #' Jaccard P (1908). Nouvellesrecherches sur la distribution florale.
# #' *Bulletin de la Societe Vaudoise des Sciences Naturelles*, 44(163).
# #' \doi{10.5169/seals-268384}
# #' 
# #' Dice, Lee R. (1945). Measures of the Amount of Ecologic Association Between
# #' Species. *Ecology*. 26 (3): 297–302. \doi{10.2307/1932409}
# #' 
# #' Hamming, R. W. (1950). Error detecting and error correcting codes. *The Bell
# #' System Technical Journal*, 29(2), 147–160. \doi{j.1538-7305.1950.tb00463.x}
# #' 
# #' Ochiai, A. (1957). Zoogeographical studies on the soleoid fishes found in
# #' Japan and its neighbouring regions. *Bulletin of the Japanese Society of
# #' Scientific Fisheries*. 23(7-8): 458-466.
# #' 
# #' Sorensen, T. (1948). A method of establishing groups of equal amplitude in
# #' plant sociology based on similarity of species and its application to
# #' analyses of the vegetation on Danish commons.
# #' *Kongelige Danske Videnskabernes Selskab*. 5 (4): 1–34.
# #' 
