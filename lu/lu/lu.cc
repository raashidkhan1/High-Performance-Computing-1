#include <chrono>
#include <typeinfo>
#include <cstdio>

#include "matrix.h"

#include "lu.h"


using namespace std;

int
main(int argc, char **argv)
{
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  int nnz, n_rows, n_cols;
  bool ok(false);
    CRSFormatMatrix cfm;

  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, n_rows, n_cols,
                          values, col_ind, row_ptr_begin, row_ptr_end);
  if (!ok)
    {
      fprintf(stderr, "failed to load matrix.\n");
      return -1;
    }

  /* For debugging, can be removed when implementation is finished. */
//  dump_nonzeros(n_rows, cfm.values, cfm.col_ind, cfm.row_ptr_begin, cfm.row_ptr_end);
    
    // initialize CRS matrix
    
    copy(begin(values), end(values), begin(cfm.values));
    cfm.n_rows = n_rows;
    cfm.n_cols = n_cols;
    copy(begin(row_ptr_begin), end(row_ptr_begin), begin(cfm.row_ptr_begin));
    copy(begin(row_ptr_end), end(row_ptr_end), begin(cfm.row_ptr_end));
    cfm.init_memory_management();
    PMatrix pm;
    pm.identity(cfm.n_rows);
    
  auto factorization_start_time = std::chrono::high_resolution_clock::now();

  /* Perform LU factorization here */
    // perform PA = LU-- elements below diagonal on A matrix- update P matrix accordingly with pivoting
    lu_factorization(cfm, pm);
    
  auto factorization_end_time = std::chrono::high_resolution_clock::now();
  
    /* Construct x vectors here */
    double x1[n_rows], x2[n_rows], x3[n_rows],  x4[n_rows], x5[n_rows];
    computeXVector(n_rows, x1, vector_pattern[0]);
    computeXVector(n_rows, x2, vector_pattern[1]);
    computeXVector(n_rows, x3, vector_pattern[2]);
    computeXVector(n_rows, x4, vector_pattern[3]);
    computeXVector(n_rows, x5, vector_pattern[4]);
    
    /* Compute b vectors here */
    double b1[n_rows], b2[n_rows], b3[n_rows], b4[n_rows], b5[n_rows];
    computeBVector(cfm, x1, b1);
//    computeBVector(cfm, x2, b2);
//    computeBVector(cfm, x3, b3);
//    computeBVector(cfm, x4, b4);
//    computeBVector(cfm, x5, b5);
    
    
  auto solve_start_time = std::chrono::high_resolution_clock::now();
  
  /* Compute all 5 solution vectors here */
//    substitution(cfm, pm, b1, x1);
//    substitution(cfm, pm, b2, x2);
//    substitution(cfm, pm, b3, x3);
//    substitution(cfm, pm, b4, x4);
//    substitution(cfm, pm, b5, x5);

  auto solve_end_time = std::chrono::high_resolution_clock::now();
  
  
  double relative_errors[test_vector_count] = {0};
  
  /* Compute relative errors here */
  
  
  std::chrono::duration<double> factorization_elapsed_time = factorization_end_time - factorization_start_time;
  std::chrono::duration<double> solve_elapsed_time = solve_end_time - solve_start_time;
  
  
  /* Print results */
  fprintf(stdout, "factorization_elapsed_time %.20f\n", factorization_elapsed_time.count());
  fprintf(stdout, "solve_elapsed_time %.20f\n", solve_elapsed_time.count());
  for (size_t vector_idx = 0; vector_idx < test_vector_count; ++vector_idx)
    {
      fprintf(stdout, "%.20f\n", relative_errors[vector_idx]);
    }

  return 0;
}
