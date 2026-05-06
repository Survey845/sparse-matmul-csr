# Sparse Matrix-Vector Multiplication using Compressed Sparse Row

This repository contains the implementation of the [RV-Sparse coding challenge](https://docs.google.com/document/d/15LcAv0bXG6-J-n7p6rQUkhout0UnwftpRGL43cbOrUg/edit?tab=t.0), developed as a foundational stepping stone for the [RISC-V Vectorized Sparse Linear Algebra library](https://github.com/merledu/rv-sparse). 

The core objective is to implement a highly efficient Sparse Matrix-Vector Multiplication routine with a strict zero dynamic memory allocation constraint. This mimics the stringent memory environments of bare-metal embedded systems and prepares the algorithmic logic for hardware-level vectorization using the RISC-V Vector Extension (RVV).

## Overview

The challenge asks us to write a C program that multiplies a sparse matrix by a vector. A sparse matrix is a large grid of numbers where almost all of the values are exactly zero. Doing math with all those zeros wastes a lot of time and computer memory. To solve this problem, the challenge requires the code to do two specific things:

1. **Compress the data**: The program scans the grid and ignores all the zeros. It saves only the useful numbers using a storage method called Compressed Sparse Row (CSR).

2. **Calculate the result**: The program multiplies only these compressed numbers by an input vector to get the final answer.

## Constraints

1. **Zero Dynamic Allocation**: The program cannot request new memory while running (no malloc, calloc, or VLAs).

2. **Pre-Allocated Buffers**: All data compression and math must use only the empty memory arrays provided by the caller at the start.

## Prerequisites

- **C Compiler**: GCC or Clang (C99 standard or newer)
- **Math Library**: Standard C math library (linked via -lm)
- **Version Control**: Git (to clone repository)

## Build & Run

1. **Clone**
```
git clone https://github.com/Survey845/sparse-matmul-csr.git
cd sparse-matmul-csr
```
2. **Compile**
```
gcc -lm -o run challenge.c
```
3. **Execute**
```
./run
```

## Implementation
1. **Data Compression**
- **Scan the matrix**: Iterate through the flat, row-major dense input matrix ```A```.
- **Extract non-zeros**: When a non-zero value is found, it is appended to the `values` array. Record its column position in `col_indices` array and increment the non-zero counter `nnz`.
- **Map row boundaries**:  Record the starting index of each row `i` by assigning `row_ptrs[i] = nnz`.
- **Cap the pointers**:  Set the final element `row_ptrs[rows] = nnz`. With this, we can calculate the length of any row `i` can be calculated as `row_ptrs[i+1] - row_ptrs[i]`.

2. **Multiplication**
- Iterate through the rows `i` using row_ptrs to determine the exact bounds of the non-zero element. 
- Calculate the dot product using the CSR and input vector `x`. The sum is written to the pre-allocated output vector `y`. The mathematical operations executed for each row is:

```
y[i] = sum( values[k] * x[col_indices[k]] )
for k = row_ptrs[i] to row_ptrs[i+1] - 1
```

## Design
- **Format Selection**: Compressed Sparse Row format was chosen over the Coordinate (COO) format because CSR is far more memory-efficient. It is also better for row slicing and faster matrix-vector multiplication.
- **Memory Decoupling**: The memory is provided by the caller for `values`, `col_indices`, `row_ptrs`, `y` which ensures safety in embedded environments.
- **Cache locality** - The CPU cache can efficiently pre-fetch data from the CSR contiguous memory. This is an important Single Instruction Multiple Data principle for Vector processors.

## Alternatives
Two other sparse formats which were considered
1. **Coordinate Format (COO)**: Stores each non-zero value as a triple `(row, column, value)`. It is the simplest format to construct but less efficient for matrix-vector multiplications it requires sorting and does not support fast row slicing. It also used more memory since it stores a full row index for every non-zero element.

2. **Compressed Sparse Columns (CSC)**: The column major equivalent of CSR. It is well-suited for column-wise access patterns. However, since matrix-vector multiplication is primarily a row wise operation, CSC would require more complex indexing and lose the cache locality advantage provided naturally by CSR.

## Complexity
- **Time Complexity**
1. **Dense to CSR**: `O(MxN)` where `M` is number of rows and `N` is number of columns
2. **Matrix-Vector Multiplication**: `O(M+NNZ)` where `M` is the number of rows and `NNZ` is the number of non-zero values

**Total Time Complexity**: `O(MxN)` as `NNZ` can never be larger than `M x N`, so the overall execution time is dominated by the Dense to CSR conversion.

- **Space**
1. **Auxillary space**: `O(1)` as only a few local variables are declared inside the function
2. **Data structure space**: `O(NNZ+M)` because `values` and `col_indices` scales with `NNZ` and `row_ptrs` and `y` scales with `M`.

**Total Space Complexity**: `O(1)` is the total space complexity as the caller is responsible for allocating memory passed to the function.
## Testing & Validation
The main function inside `challenge.c` acts as a test harness, proving the mathematical accuracy and stability of the sparse_multiply function.
1. **Iterations**: 100 independent test cycles
2. **Randomized Values**: Generate a matrix between `5 x 5` and `45 x 45` and a random sparsity density (5% - 40%)
3. **True Answer**:  Calculates a reference vector output using a standard dense matrix multiplication
4. **Error Tolerance**:  To account for floating-point drift, it used mixed absolute and relative tolerance formula: `1e-7 + 1e-7 * fabs(y_ref[i])`
5. **Exit Conditions**: The program logs the pass/fail status of every iteration and will only exit with a success code (0) if all 100 randomized matrices compute perfectly.

## Results
This code passes all the automated tests run inside the main function block

- **Iteration**

```
Iter  0 [ 19x 45, density=0.22, nnz= 199]: PASS (Max error: 0.00e+00)
```
- **Final**

```
All tests passed! (100/100 iterations passed)
```

## Usage
The sparse multiply function has the following signature:
```
void sparse_multiply(int rows, int cols, const double *A, const double *x,
                     int *out_nnz, double *values, int *col_indices,
                     int *row_ptrs, double *y);
```
Where the user is responsible for allocating all buffers before calling the function
- `values` and `col_indices` are at least `rows x columns`
- `row_ptrs` is at least `rows + 1`
- `y` is at least `rows`

**Minimal working example**:
```
int rows = 3, cols = 4;

double A[] = {
    1.0, 0.0, 0.0, 2.0,
    0.0, 0.0, 3.0, 0.0,
    0.0, 4.0, 0.0, 5.0
};
double x[] = { 1.0, 2.0, 3.0, 4.0 };

double values[12];
int    col_indices[12];
int    row_ptrs[4];
double y[3];
int    out_nnz = 0;

sparse_multiply(rows, cols, A, x, &out_nnz, values, col_indices, row_ptrs, y);
```

## Future Applications
- **RISC-V Vector Extension Acceleration**: The CSR traversal loop is a natural candidate for vectorization using RVV intrinsics. The contiguous memory layout allows vector load instructions to process multiple non-zeros per clock cycle for data level parallelism on RVV enabled processors.

- **Scientific computing and simulation**: It is a core operation in iterative solvers like Conjugate Gradient and GMRES, used in large-scale physics simulation and computational fluid dynamics

- **Lightweight sparse BLAS library**:u Part of a clean, portable and lightweight sparse linear algebra library for the RISC-V ecosystem, filling the gap between dense-optimized OpenBLAS libraries and equivalent sparse support for RVV enabled processors.

## Reading Material
- [RISC-V Vector ISA](https://riscv.org/wp-content/uploads/2024/12/15.20-15.55-18.05.06.VEXT-bcn-v1.pdf)
- [RV-Sparse library](https://github.com/merledu/rv-sparse)
- [Sparse Matrices: COO v/s CSR](https://medium.com/@oz3dprinter/mastering-sparse-matrices-a-deep-dive-into-coo-and-csr-formats-for-python-practitioners-09cfda501ffa)
- [Vector Processors](https://enesharman.medium.com/vector-processors-f1b082ee9cdb)
- [Conjugate Gradient using CSR](https://www.sciencedirect.com/science/article/pii/S0377042711002196#s000030)
- [GMRES using CSR](https://jomardpublishing.com/UploadFiles/Files/journals/AMMAV1N1/V7N2/Boumzough_et_al.pdf)
