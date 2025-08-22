test_that("alpha diversity", {
  
  # Chao1 =====================================================================
  
  expect_equal(
    object   = chao1(counts), 
    expected = c(A = 3, B = NaN, C = NaN,  D = NaN) )
  
  expect_equal(chao1(1:10), 10.5)
  expect_equal(chao1(rep(1:10, 2)), 21)
  
  
  
  # Faith's Phylogenetic Diversity ============================================
 
  expect_equal( # abdiv::apply(counts, 2L, abdiv::faith_pd, tree)
    object   = faith(counts, tree), 
    expected = c(A = 2.319, B = 2.191, C = 2.191, D = 1.759) )
  
  
  
  # Inverse Simpson ===========================================================
 
  expect_equal( # vegan::diversity(t(counts), 'invsimpson')
    object   = inv_simpson(counts), 
    expected = c(A = 2.6, B = 2.84210526315789, 
                 C = 2.84516129032258, D = 1.8 ))
  
  
  
  # Shannon ===================================================================
 
  expect_equal( # vegan::diversity(t(counts), 'shannon')
    object   = shannon(counts), 
    expected = c(A = 1.01233083910317, B = 1.07204334357507,
                 C = 1.07101854240991, D = 0.636514168294813 ))
  
  
  
  # Simpson ===================================================================
 
  expect_equal( # vegan::diversity(t(counts), 'simpson')
    object   = simpson(counts), 
    expected = c(A = 0.615384615384615, B = 0.648148148148148, 
                 C = 0.648526077097506, D = 0.444444444444444 ))
  
  
  
  # Matrix with > 100 columns to trigger pthreading ===========================
  
  expect_silent(simpson(big_mtx))
  expect_silent(faith(big_mtx, tree))
  
})
