// Copyright (c) 2025 ecodive authors
// Licensed under the MIT License: https://opensource.org/license/mit


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // calloc, free
#include <string.h> // memset
#include "get.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  int    edge;
  int    parent;
  double length;
} node_t;


//======================================================
// Variables shared between main and worker threads.
//======================================================
static double *otu_mtx;
static int     n_otus;
static int     n_samples;
static int     n_edges;
static double *edge_lengths;
static node_t *nodes;
static int     n_threads;
static double *result_vec;
static int     mtx_len;



//======================================================
// Find the edges that an otu passes through.
//======================================================

static void *faith_result(void *arg) {
  
  int   thread_i     = *((int *) arg);
  char *has_edge_vec = calloc(n_edges, sizeof(char));
  
  if (has_edge_vec == NULL) { // # nocov start
    free(has_edge_vec);
    error("Unable to allocate memory for Faith PD calculation.");
    return NULL;
  } // # nocov end
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    memset(has_edge_vec, 0, n_edges * sizeof(char));
    
    // mark the tree edges that this OTU uses
    int otu = 0;
    for (int i = sample; i < mtx_len; i += n_samples) {
      if (otu_mtx[i] > 0) {           // OTU present in sample
        int node = otu;               // start at OTU tip/leaf in tree
        while (node > -1) {           // traverse until we hit the tree's root
          char *has_edge = has_edge_vec + nodes[node].edge;
          if (*has_edge) break;       // already traversed
          *has_edge = 1;
          node = nodes[node].parent;  // proceed on up the tree
        }
      }
      otu++;
    }
    
    double Faith = 0;
    for (int edge = 0; edge < n_edges; edge++) {
      Faith += has_edge_vec[edge] * edge_lengths[edge];
    }
    result_vec[sample] = Faith;
    
  }
  
  free(has_edge_vec);
  
  return NULL;
}



//======================================================
// R interface. Dispatches threads.
//======================================================
SEXP C_faith(
    SEXP sexp_otu_mtx,  SEXP sexp_phylo_tree,  
    SEXP sexp_n_threads ) {
  
  otu_mtx   = REAL( sexp_otu_mtx);
  n_otus    = ncols(sexp_otu_mtx);
  n_samples = nrows(sexp_otu_mtx);
  
  mtx_len = n_otus * n_samples;
  
  int *edge_mtx = INTEGER(get(sexp_phylo_tree, "edge"));
  n_edges       = nrows(  get(sexp_phylo_tree, "edge"));
  edge_lengths  = REAL(   get(sexp_phylo_tree, "edge.length"));
  
  n_threads = asInteger(sexp_n_threads);
  
  
  // intermediary values
  nodes = calloc(n_edges, sizeof(node_t));
  
  if (nodes == NULL) { // # nocov start
    free(nodes);
    error("Unable to allocate memory for Faith PD calculation.");
    return R_NilValue;
  } // # nocov end
  
  
  // sort edge data by child node
  for (int edge = 0; edge < n_edges; edge++) {
    
    int    parent = edge_mtx[0 * n_edges + edge] - 2;
    int    child  = edge_mtx[1 * n_edges + edge] - 1;
    double length = edge_lengths[edge];
    
    if (child  > n_otus) child--;
    if (parent < n_otus) parent = -1;
    
    nodes[child].edge   = edge;
    nodes[child].parent = parent;
    nodes[child].length = length;
  }
  
  
  // Create the diversity vector to return
  SEXP sexp_result_vec = PROTECT(allocVector(REALSXP, n_samples));
  result_vec           = REAL(sexp_result_vec);
  SEXP sexp_mtx_dimnames = getAttrib(sexp_otu_mtx, R_DimNamesSymbol);
  if (sexp_mtx_dimnames != R_NilValue) {
    SEXP sexp_mtx_rownames = VECTOR_ELT(sexp_mtx_dimnames, 0);
    if (sexp_mtx_rownames != R_NilValue) {
      setAttrib(sexp_result_vec, R_NamesSymbol, sexp_mtx_rownames);
    }
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1 && n_samples > 100) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      int       *args = calloc(n_threads, sizeof(int));
      
      int i, n = n_threads;
      for (i = 0; i < n; i++) args[i] = i;
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, faith_result, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      free(tids); free(args); free(nodes);
      
      UNPROTECT(1);
      return sexp_result_vec;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  n_threads    = 1;
  int thread_i = 0;
  faith_result(&thread_i);
  
  free(nodes);
  
  UNPROTECT(1);
  return sexp_result_vec;
}

