#ifndef lu_h
#define lu_h

#include <cstdio>
#include <iostream>
#include <cmath>

/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

using namespace std;

// Test Execution parameters
//const int max_n_elements = 131072; //initial
const int max_n_elements = 1000000; //1M prevents segfault
const int max_n_rows = 16384;
const int test_vector_count = 5;

// Compressed Row Storage Format
static double values[max_n_elements];
static int col_ind[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];
int nnz, n_rows, n_cols;

// Utility matrices
static double pb[max_n_rows];
static double y[max_n_rows];
static int use_index[max_n_rows];
static double dense_values[max_n_rows];
static int permutation_index[max_n_rows];



/* Initialization functions **/
// compute x vectors
void init_solution_vector(int n_rows, double x[], double vector_pattern[]){
    for(int i=0;i<n_rows;i++){
        x[i] = vector_pattern[i%2];
    }
}

// multiply matrix with vector
void multiply_matrix_vector(double vector[], double product_vector[]){
    for (int row=0; row<n_rows; row++){
      product_vector[row] = 0;
      for (int i=row_ptr_begin[row]; i<=row_ptr_end[row]; i++){
        auto value = values[i];
        auto column = col_ind[i];
          product_vector[row] += value * vector[column];
      }
    }
}

// compute b vectors from A.x = b
void init_b_vector(double x[], double b[]){
    multiply_matrix_vector(x, b);
}

// initialize permutation matrix
void init_permutation_identity_matrix(int n_rows){
    for(int row_idx=0; row_idx<n_rows; row_idx++){
      permutation_index[row_idx] = row_idx;
    }
}


/* Utility functions**/

// compute relative errors
double compute_relative_errors(double solution_vector_x[], double x[]){
    
    double x_denominator = 0, x_numerator= 0, error=0;
    for(int i=0;i<n_rows;i++){
        if(solution_vector_x[i] == x[i]){
            continue;
        }
        x_denominator += solution_vector_x[i]*solution_vector_x[i];
        x_numerator += (x[i]-solution_vector_x[i])*(x[i]-solution_vector_x[i]);
    }
    // check for valid numerator and denominator to avoid nan error output
    if(x_denominator!=0 && x_numerator!=0){
        x_numerator = sqrt(x_numerator);
        x_denominator = sqrt(x_denominator);
        if(x_numerator>0 && x_denominator>0){
            error = x_numerator/x_denominator;
        }
    }
    return error;
}

// Swap rows
void swap_rows(int pivot_row_index, int swap_row){
    auto row_begin = row_ptr_begin[pivot_row_index];
    row_ptr_begin[pivot_row_index] = row_ptr_begin[swap_row];
    row_ptr_begin[swap_row] = row_begin;

    auto row_end = row_ptr_end[pivot_row_index];
    row_ptr_end[pivot_row_index] = row_ptr_end[swap_row];
    row_ptr_end[swap_row] = row_end;
}

// find column with the value
bool find_column(int search_row, int search_column, int actual_index){
    
  for(actual_index = row_ptr_begin[search_row];
      actual_index<=row_ptr_end[search_row];
      actual_index++){

    auto found_column = col_ind[actual_index];

    if (found_column == search_column){
      return true;
    }
  }
  return false;
}


/* Factorization functions**/

// find pivot and swap rows
double partial_pivoting(int pivot_row_index){
    int swap_row = pivot_row_index;
    double pivot_value = 0.0;
    
    int actual_index = 0;
    
    if (find_column(pivot_row_index, pivot_row_index, actual_index)){
      pivot_value = values[actual_index];
    } else {
      pivot_value = 0;
    }
    
    //find the row with the largest pivot value
    for (int change_row=pivot_row_index+1; change_row<n_rows; change_row++){
      for (int actual_index = row_ptr_begin[change_row];
          actual_index <= row_ptr_end[change_row];
          actual_index++){
        
        int column_index = col_ind[actual_index];
        double value = values[actual_index];
      
        if (column_index == pivot_row_index){
            // swap if value is greater than current pivot value
          if (abs(value) > abs(pivot_value)){
            pivot_value = value;
            swap_row = change_row;
          }
          break;
        }
      }
    }

    swap_rows(pivot_row_index, swap_row);
    
    // save the permutation in a separate matrix
    permutation_index[pivot_row_index] = swap_row;
    permutation_index[swap_row] = pivot_row_index;
    
    return pivot_value;
}

// Eliminate values below the pivot row
void elimination(int pivot_row, int add_row, double pivot){
    auto pivot_column = pivot_row;
    int pivot_column_in_other_row = 0;
    
    bool found = find_column(add_row, pivot_column, pivot_column_in_other_row);

    if (!found){
        return; // return if there is no non-zero column below the pivot
    }

    int last_index = 0;
    
    /* Apply scattering and gathering */
    //scatter into dense row
    for(int k=row_ptr_begin[add_row]; k<=row_ptr_end[add_row]; k++){
      auto value = values[k];
      auto column = col_ind[k];
      dense_values[column] = value;
      last_index = column;
    }
    
    auto mult = values[pivot_column_in_other_row]/pivot;
    
    for(int k=row_ptr_begin[pivot_row]+1; k<=row_ptr_end[pivot_row]; k++){
      auto column = col_ind[k];
      if (column > pivot_column){
          dense_values[column] -= mult*values[k];
      }
    }
    //set the element in the column below the pivot
    dense_values[pivot_column] = mult;

    //gather into compressed sparse row
    int actual_index = row_ptr_begin[add_row];
    for(int column=0; column<n_rows; column++){
      double value = dense_values[column];
      if (value==0){ //skip zero values
        continue;
      }
      values[actual_index] = value;
      col_ind[actual_index] = column;
      actual_index++;
    }
    //update row ptr in case row gets small
    row_ptr_end[add_row] = actual_index-1;
    
    //garbage collection: reset dense_values by filling with zeroes
    fill_n(dense_values, last_index+1, 0);
}

// PA = LU
void lu_factorization(){
    auto n_columns = n_rows; //square matrix
    for(int column = 0; column<n_columns; column++){
        double pivot = partial_pivoting(column);
        if (pivot == 0) {
            continue; //continue if current pivot is 0
        }
      //run elimination for each row below the current pivot
      for(int row = column+1; row<n_rows; row++){
        elimination(column, row, pivot);
      }
    }
}

// Utility function to permute b vector similar as LU
void permute_b(double b[], double pb[],
                         int length){
  //copy the array into pb
  for (int i=0; i<length; i++){
      pb[i] = b[i];
  }

  //create index
  for(int i=0; i<length; i++){
      use_index[i]=i;
  }

  for(int i=0; i<length; i++){
    auto v = permutation_index[i];
    auto temp = use_index[i];
      use_index[i] = use_index[v];
      use_index[v] = temp;
  }

  //permute vector b same as LU vector
  for(int i=0; i<length; i++){
    pb[i] = b[use_index[i]];
  }
}

// upon factorization
// Ly = Pb -- forward substitution
// Ux = y -- backward substitution
void substitution(double b[], double x[]){

    //permute vector b same as LU vector
    permute_b(b, pb, n_rows);

    // forward subsitution, gets y from Ly=Pb
    for(int row=0; row<n_rows; row++){
        y[row] = pb[row];
        for(int actual_index=row_ptr_begin[row]; actual_index<=row_ptr_end[row]; actual_index++){

            //only iterate till column<row-1
            auto column = col_ind[actual_index];
            if(column>row-1){
                break;
            }

            auto A = values[actual_index];
            y[row] -= A * y[column];
        }
    }

    //  back subsitution, gets x from Ux=y
    for(int row = n_rows-1; row>=0; row--){
        x[row] = y[row];
        for (int actual_index=row_ptr_begin[row]; actual_index<=row_ptr_end[row]; actual_index++){
          if (col_ind[actual_index]<row+1){
              continue; // continue if we are not in the next row
          }

          double A = values[actual_index];
          auto column = col_ind[actual_index];
          x[row] -= A * x[column];
        }
        double A;
        int actual_index = 0;
        if (find_column(row,row,actual_index)){
          A = values[actual_index];
        } else {
          A = 1.0;
        }
        x[row] = x[row] / A;
    }
}


#endif /* lu_h */
