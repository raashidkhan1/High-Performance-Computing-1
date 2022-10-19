//
//  lu.h
//  lu
#ifndef lu_h
#define lu_h

#include <cstdio>
#include <ctime>
#include <iostream>
#include <cstring>
#include <cmath>

/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

//const int max_n_elements = 131072; //initial
const int max_n_elements = 12000000; //initial
const int max_n_rows = 16384;
const int test_vector_count = 5;

const int max_n_cols = max_n_rows; //max colums = max rows ---Square matrix

static double values[max_n_elements];
static int col_ind[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];

double vector_pattern[test_vector_count][2] = {
    {1, 1},
    {0.1, 0.1},
    {1, -1},
    {5, -5},
    {100, -100}
};

// a template to store common functionalities for Compressed Row Storage Format Matrix
struct CRSFormatMatrix{
    int n_rows = 0, n_cols = 0;
    int col_ind[max_n_elements];
    double values[max_n_elements];
    int row_ptr_begin[max_n_rows];
    int row_ptr_end[max_n_rows];
    int numberOfElementsInRow(int row){
        return row_ptr_end[row] - row_ptr_begin[row]+1;
    };
    int row_ptr_reserved[max_n_cols];
    int free = 0; //points to the first free element in the memory array
    int old_space_end = max_n_elements; //end of fromspace
    // Swap rows
    void swap_rows(int pivot_row_index, int swap_row){
        auto row_begin = row_ptr_begin[pivot_row_index];
        row_ptr_begin[pivot_row_index] = row_ptr_begin[swap_row];
        row_ptr_begin[swap_row] = row_begin;

        auto row_end = row_ptr_end[pivot_row_index];
        row_ptr_end[pivot_row_index] = row_ptr_end[swap_row];
        row_ptr_end[swap_row] = row_end;

        auto row_reserved = row_ptr_reserved[pivot_row_index];
        row_ptr_reserved[pivot_row_index] = row_ptr_reserved[swap_row];
        row_ptr_reserved[swap_row] = row_reserved;
    }
    
    void allocate(int numb_elements, int ptr_begin, int ptr_reserved){
        if(old_space_end-free<(numb_elements)){
          //do stop and copy to compact
          //stop_and_copy(); //TODO re-enable for larger matrices
        }
        
        ptr_begin = free;
        ptr_reserved = free + numb_elements;
        free += numb_elements +1;
      }
    void init_memory_management(){
      free = row_ptr_end[n_rows-1]+1;
      for(int i=0; i<n_rows; i++){
        row_ptr_reserved[i] = row_ptr_end[i];
      }
    }
};

// a template to store common functionalities for Permutation Matrix
struct PMatrix{
    long permuted_to_original_index[max_n_rows];
    void identity(size_t n_rows){
        for(size_t row_idx=0; row_idx<n_rows; row_idx++){
          permuted_to_original_index[row_idx] = row_idx;
        }
    }
    void mark_swap(const int row, const int swap_row){
        permuted_to_original_index[row] = swap_row;
        permuted_to_original_index[swap_row] = row;
    }
};

struct DenseIndexedRow {
  double values[max_n_cols] = {0};
};


bool find_column(CRSFormatMatrix m, int haystack_row,
                 int needle_column, int flat_index){
  for(flat_index = m.row_ptr_begin[haystack_row];
      flat_index<=m.row_ptr_end[haystack_row];
      flat_index++){

    auto hay_column = m.col_ind[flat_index];

    if (hay_column == needle_column){
      return true;
    }
  }
  return false;
}

// find pivot on the CSR matrix and swap rows
double partial_pivoting(CRSFormatMatrix cfm, PMatrix pm, int pivot_row_index){
    int swap_row = pivot_row_index;
    double largest_val = 0.0;
    
    int flat_index = 0;
    if (find_column(cfm, pivot_row_index, pivot_row_index, flat_index)){
      largest_val = cfm.values[flat_index];
    } else {
      largest_val = 0;
    }
//
        //find the row with the best pivot value
    for (int replacement_row=pivot_row_index+1; replacement_row<cfm.n_rows; replacement_row++){
      for (int flat_index = cfm.row_ptr_begin[replacement_row];
          flat_index <= cfm.row_ptr_end[replacement_row];
          flat_index++){
        
        int column_index = cfm.col_ind[flat_index];
        double value = cfm.values[flat_index];
      
        if (column_index == pivot_row_index){

          if (std::abs(value) > std::abs(largest_val)){
            largest_val = value;
            swap_row = replacement_row;
          }
          break;
        }
      }
    }

    cfm.swap_rows(pivot_row_index, swap_row);

    return largest_val;
}


void overwrite_sparse_row_with_dense(DenseIndexedRow dense_row,
                                     int row, CRSFormatMatrix lu){

  //count number of non zeros
  int nnz = 0;
  for(int i=0; i<lu.n_rows; i++){
    if (dense_row.values[i] != 0.0){
      nnz++;
    }
  }

  //+1 as lu.row_ptr_reserved points to the last element still reserved for
  //this row. thus reserved-row is one smaller then the reserved number of elements
  if(lu.row_ptr_reserved[row] - lu.row_ptr_begin[row] +1 < nnz){
    //allocate more space
    //opt, find way to just extend reserved space if there is free space
    lu.allocate(nnz, lu.row_ptr_begin[row],
                lu.row_ptr_reserved[row]);

    //copy values not in dense_row not needed as dense row now contains everything
  }

  //gather into sparse row
  auto flat_index = lu.row_ptr_begin[row];
  for(int column=0; column<lu.n_rows; column++){
    auto value = dense_row.values[column];
    if (value==0){ //skip zero values
      continue;
    }

    lu.values[flat_index] = value;
    lu.col_ind[flat_index] = column;
    flat_index++;
  }
  //update row ptr in case row shrunk
  lu.row_ptr_end[row] = flat_index-1;
}

void elimination(CRSFormatMatrix cfm, int pivot_row, int add_row, double pivot){
    //walk until we get the flat index of the pivot column in the other row
    auto pivot_column = pivot_row;
    int pivot_column_in_other_row = 0; //set by find column
    bool found = find_column(cfm, add_row, pivot_column, pivot_column_in_other_row);
    if (!found){return;} // no pivot in this column => we are done

    DenseIndexedRow dense_row; //opt make static, and reset each run

    //scatter elements after mult of other_row to a temp row in dense form
    for(int k=cfm.row_ptr_begin[add_row]; k<=cfm.row_ptr_end[add_row]; k++){
      auto value = cfm.values[k];
      auto column = cfm.col_ind[k];
      dense_row.values[column] = value;
    }
    
    auto mult = cfm.values[pivot_column_in_other_row]/pivot;
    //add scaled pivot_row to scatterd other row in dense form
    for(int k=cfm.row_ptr_begin[pivot_row]+1; k<=cfm.row_ptr_end[pivot_row]; k++){
      auto column = cfm.col_ind[k];
      if (column > pivot_column){
        dense_row.values[column] -= mult*cfm.values[k];
      }
    }
    //set the element in the column below the pivot
    dense_row.values[pivot_column] = mult;

    //printf("mult: %f, dense row [%i]: %f %f %f %f\n", mult, other_row, dense_row.values[0], dense_row.values[1], dense_row.values[2], dense_row.values[3]);
//    overwrite_sparse_row_with_dense(dense_row, add_row, cfm);
}



// PA = LU
void lu_factorization(CRSFormatMatrix cfm, PMatrix pm){
    auto n_columns = cfm.n_rows;
    for(int column=0; column<n_columns; column++){
        double pivot = partial_pivoting(cfm, pm, column);
      if (pivot==0) {continue;} //skip lines with pivot value 0
      //for each row below the current pivot
      for(int row = column+1; row<cfm.n_rows; row++){
        elimination(cfm, column, row, pivot);
      }
    }
}

// print iterator
void print_matrix();

// multiply CRS matrix with vector
void multiply_matrix_vector(CRSFormatMatrix cfm, double vector[], double product_vector[]){
    for (int row=0; row<cfm.n_rows; row++){
      product_vector[row] = 0;
      for (int i=cfm.row_ptr_begin[row]; i<=cfm.row_ptr_end[row]; i++){
        auto value = cfm.values[i];
        auto column = cfm.col_ind[i];
          product_vector[row] += value * vector[column];
      }
    }
}

// regular sparse matrix-matrix multiplication
void multiply_matrix(double a[][max_n_cols], double b[][max_n_cols], double mul[][max_n_cols]){
    for(int i=0;i<max_n_cols;i++){
        for(int j=0;j<max_n_cols;j++){
            mul[i][j] = 0;
            for(int k=0;k<max_n_cols;k++){
                mul[i][j] += a[i][k]*b[k][i];
            }
        }
    }
}

/* Initialization functions**/
// compute x vectors
void computeXVector(int n_rows, double x[], double vector_pattern[]){
    for(int i=0;i<n_rows;i++){
        x[i] = vector_pattern[i%2];
    }
}

// compute b vectors from A.x = b
// here cfm.values[]*x = b
void computeBVector(CRSFormatMatrix cfm, double x[], double b[]){
    multiply_matrix_vector(cfm, x, b);
}

void get_permuted_vector(const double array[], double pb[],
                         int length, PMatrix p){
  //copy the array into pb
  for (int i=0; i<length; i++){pb[i] = array[i];}

  //create index from permutationMatrix
  int index[max_n_rows];
  for(int i=0; i<length; i++){index[i]=i;}

  for(int i=0; i<length; i++){
    auto v = p.permuted_to_original_index[i];
    auto temp = index[i];
    index[i] = index[v];
    index[v] = temp;
  }

  //permute vector b to match up LU vect.
  for(int i=0; i<length; i++){
    pb[i] = array[index[i]];
  }
}

// upon factorisation
// Ly = Pb -- forward substitution
// Ux = y -- backward substitution
void substitution(CRSFormatMatrix lu, PMatrix p,
                         double org_b[], double x[]){
    double pb[max_n_rows];
    double c[max_n_rows];
    //permute vector b to match up LU vect.
    get_permuted_vector(org_b, pb, lu.n_rows, p);

      //forward subsitution, solves y from Ly=Pb
      for(int row=0; row<lu.n_rows; row++){
          c[row] = pb[row];

          //iterate over all columns of A
          //if A[] is zero nothing happens to c =>
          //we only need to iterate the nonzero elements
          //for(int column=0; column<row-1; column++){
          for(int flat_idx=lu.row_ptr_begin[row];
              flat_idx<=lu.row_ptr_end[row]; flat_idx++){

              //only iterate till column<row-1
              auto column = lu.col_ind[flat_idx];
              if(column>row-1){break;}

              auto A = lu.values[flat_idx];
              c[row] -= A * c[column];
          }
      }

      //back subsitution, solves x from Ux=y
      for(int row = lu.n_rows-1; row>=0; row--){
          x[row] = c[row];
          //iterate over all columns of A
          //if A[] is zero nothing happens to c =>
          //we only need to iterate the nonzero elements
          //for(int column=row+1; column<row; column++){
          
          //find first non zero element in row of A starting at column row+1
          for (int flat_idx=lu.row_ptr_begin[row]; flat_idx<=lu.row_ptr_end[row]; flat_idx++){
            //only truely start iterating over columns once we are past column row+1
            if (lu.col_ind[flat_idx]<row+1){ continue; }

            double A = lu.values[flat_idx];
            auto column = lu.col_ind[flat_idx];
            x[row] -= A * x[column];
          }
          double A;
          int flat_idx = 0;
          if (find_column(lu,row,row,flat_idx)){
            A = lu.values[flat_idx];
          } else {
  //          dbg(row);
            A = 1.;
          }
          x[row] = x[row] / A;
      }
}


#endif /* lu_h */
