#include <chrono>
#include <typeinfo>
#include <cstdio>

#include "matrix.h"

#include "lu.h"

// provided solution vectors
static double solution_vector_x1[max_n_rows], solution_vector_x2[max_n_rows], solution_vector_x3[max_n_rows],  solution_vector_x4[max_n_rows], solution_vector_x5[max_n_rows];

// b vectors computed from b = A*solution_vector_x
static double b1[max_n_rows], b2[max_n_rows], b3[max_n_rows], b4[max_n_rows], b5[max_n_rows];

// x vectors for solving linear system
static double x1[max_n_rows], x2[max_n_rows], x3[max_n_rows], x4[max_n_rows], x5[max_n_rows];

// provided 5 solution vector patterns
double vector_pattern[5][2] = {
    {1, 1},
    {0.1, 0.1},
    {1, -1},
    {5, -5},
    {100, -100}
};


int main(int argc, char **argv)
{
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  bool ok(false);
    

  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, n_rows, n_cols,
                          values, col_ind, row_ptr_begin, row_ptr_end);
  if (!ok)
    {
      fprintf(stderr, "failed to load matrix.\n");
      return -1;
    }

  /* For debugging, can be removed when implementation is finished. */
//  dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);
    
    // initialize permutation matrix
    init_permutation_identity_matrix(n_rows);
    
    /* initialize solution vectors here */
    
    init_solution_vector(n_rows, solution_vector_x1, vector_pattern[0]);
    init_solution_vector(n_rows, solution_vector_x2, vector_pattern[1]);
    init_solution_vector(n_rows, solution_vector_x3, vector_pattern[2]);
    init_solution_vector(n_rows, solution_vector_x4, vector_pattern[3]);
    init_solution_vector(n_rows, solution_vector_x5, vector_pattern[4]);
    
    /* initialize b vectors here */
    
    init_b_vector(solution_vector_x1, b1);
    init_b_vector(solution_vector_x2, b2);
    init_b_vector(solution_vector_x3, b3);
    init_b_vector(solution_vector_x4, b4);
    init_b_vector(solution_vector_x5, b5);
    
  auto factorization_start_time = std::chrono::high_resolution_clock::now();

  /* Perform LU factorization here
   perform PA = LU-- elements below diagonal on A matrix- update P matrix accordingly with pivoting
   */
    lu_factorization();
    
  auto factorization_end_time = std::chrono::high_resolution_clock::now();
  
    

  auto solve_start_time = std::chrono::high_resolution_clock::now();
  
  /* Solve the linear system here for each b vector by forward and backward substitution */
    substitution(b1, x1);
    substitution(b2, x2);
    substitution(b3, x3);
    substitution(b4, x4);
    substitution(b5, x5);

  auto solve_end_time = std::chrono::high_resolution_clock::now();
  
  
  double relative_errors[test_vector_count] = {0};
  
  /* Compute relative errors here */
    relative_errors[0] = compute_relative_errors(solution_vector_x1, x1);
    relative_errors[1] = compute_relative_errors(solution_vector_x2, x2);
    relative_errors[2] = compute_relative_errors(solution_vector_x3, x3);
    relative_errors[3] = compute_relative_errors(solution_vector_x4, x4);
    relative_errors[4] = compute_relative_errors(solution_vector_x5, x5);
    
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
