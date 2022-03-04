# about matrices
operations on matrices,
useful in linear algebra


 
  This jar file contains class "matrices",
  which does various kinds of operations on matrices in a form of double[][]
  this class contains some public and private methods,
  public methods are key methods they do the main operations about matrices
  private methods are some helpful methods for public methods to provide decomposed and clean code.
  all the public methods which provide you with the ability to easily do operations on matrices are:
 
 * printMatrix() --> prints matrix in a console program of ACM library
 * identityMatrix(rank) --> crates and returns n*n identity matrix of n rank
 * multiplyTwo(array1, array2) --> multiplies two matrices from left to right and returns it { column number of the left matrix must be equal to row number of right matrix}
 * multiply(array...) --> this method is generalized version of multiplyTwo { column number of the left matrix must be equal to row number of right matrix}
 * multiplyByConstant(int K, array) --> this method multiplies array by constant and returns result
 * addTwo(arr1, arr2) --> adds two matrices and returns it { sizes of the matrices must be equal, this method does not change matrices passed by argument}
 * add(array...) --> this method is generalized version of addTwo { sizes of the matrices must be equal}
 * subtract(arr1, arr2) --> subtracts second matrix from first 
 * transposeMatrix(array) --> returns transposes matrix passed by argument and returns it
 * inverse(array) --> returns inverse of a matrix 
 * determinant(array) --> returns determinant of a matrix as a double 
 * cofactorMatrix(array) --> returns co-factor matrix 
 * upperTriangularMatrix(array) --> returns upper triangular matrix { and also saves the matrix needed to get this transformation }
 * lowerTriangularMatrixFromTranspose(array) --> finds upper triangular matrix and then transposes it
 * upperTriangularMatrixForRREF(array) --> returns upper triangular matrix that has all the zero rows (if any) at the bottom
 * RREF(array) --> returns row reduced echelon form of a matrix and also saves 
 * nullSpaceBasisMatrix(array) --> finds basis vectors for nullSpace and puts it in a matrix { and also saves nonPivotColumns, pivotColumns and rank as global variables }
 * projectionMatrix(array) --> finds n*n projection matrix to project on the column space of n*m matrix
 * projection(array, arrayOfVectors) --> finds projection of the K vectors (from n*K arrayOfVectors) on the column space of n*m array and returns (n*K arrayOfProjections which contains) K projections for each vector { also saves the error vector }
 * error(array, arrayOfVectors) --> returns arrayOfVectors - projection(array, arrayOfVectors)
 
 
 here are some not so famous operations but anyway public:
 
 * removeColumn(array, int columnToRemove) --> removes column from an array and returns sub matrix { this method trusts the user that passed arguments will be correct: 0 <= columnToRemove < total numbers of columns}  
 * findSubMatrix(array, int rowToRemove, int columnToRemove) --> if n*m matrix is passed as an argument return result is (n-1)*(m-1) with removed row and column  { this method trusts the user that passed arguments will be correct: 0 <= rows < total number of rows, 0 <= column < total numbers of columns} 
 * subtractRowsToGetZeroBelow(array, int mainRow, int column, int rowToSubtractFrom) --> subtracts K * mainRow from rowToSubtract and returns matrix which has zero below array[mainRow][column} entry
 * switchRows(array, int row1, int row2) --> switches row1 and row2 and returns matrix with switched rows
 * isNonZeroInRow(array, int row) --> finds if there is nonzero number in a row
 * findNonZeroCoefficient(array, int row) --> finds nonzero number in a row and returns its coefficient
 
