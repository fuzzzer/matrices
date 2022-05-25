
/* ABOUT 
 * This program contains class matrices,
 * which does various kinds of operations on matrices in a form of a double[][]
 * this class contains some public and private methods,
 * public methods are key methods they do the main operations about matrices
 * private methods are some helpful methods for public methods to provide decomposed and clean code.
 * all the public methods which provide you with the ability to easily do operations on matrices are:
 * 
 * printMatrix() --> prints matrix in a console program of ACM library

 * identityMatrix(rank) --> crates and returns n*n identity matrix of n rank

 * multiplyTwo(array1, array2) --> multiplies two matrices from left to right and returns it { column number of the left matrix must be equal to row number of right matrix}

 * multiply(array...) --> this method is generalized version of multiplyTwo { column number of the left matrix must be equal to row number of right matrix}

 * multiplyByConstant(int K, array) --> this method multiplies array by constant and returns result

 * addTwo(arr1, arr2) --> adds two matrices returns it { sizes of the matrices must be equal, this method does not change matrices passed by argument}

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

 * nullSpaceBasisMatrix(array) --> finds basis vectors for nullSpace and puts it in a matrix { and also saves nonPivotColumns, pivotColumns and rank as a global variables }

 * projectionMatrix(array) --> finds n*n projection matrix to project on the column space of n*m matrix

 * projection(array, arrayOfVectors) --> finds projection of the K vectors (from n*K arrayOfVectors) on the column space of n*m array and returns (n*K arrayOfProjections which contains) K projections for each vector { also saves the error vector }

 * error(array, arrayOfVectors) --> returns arrayOfVectors - projection(array, arrayOfVectors)

 * 
 * 
 * 
 * some not so popular, but anyway public methods:

 * removeColumn(array, int columnToRemove) --> removes column from an array and returns sub matrix { this method trusts the user that passed arguments will be correct: 0 <= columnToRemove < total numbers of columns}  
 * findSubMatrix(array, int rowToRemove, int columnToRemove) --> if n*m matrix is passed as an argument return result is (n-1)*(m-1) with removed row and column  { this method trusts the user that passed arguments will be correct: 0 <= rows < total number of rows, 0 <= column < total numbers of columns} 
 * subtractRowsToGetZeroBelow(array, int mainRow, int column, int rowToSubtractFrom) --> subtracts K * mainRow from rowToSubtract and returns matrix which has zero below array[mainRow][column} entry
 * switchRows(array, int row1, int row2) --> switches row1 and row2 and returns matrix with switched rows
 * isNonZeroInRow(array, int row) --> finds if there is nonzero number in a row
 * findNonZeroCoefficient(array, int row) --> finds nonzero number in a row and returns its coefficient
 * 
 */



import java.util.*;
import java.lang.Math;

public class matrices  {
	
	//-------------------MAIN METHODS (0)--------------------------------
	
	private void printMatrix(double[][] arr) {
		System.out.println("\n");
		for(int i=0; i<arr.length ; i++) {
			String row = "";
			for(int j=0 ; j<arr[0].length ; j++) {
				row = row + Math.round(arr[i][j] * 100) / 100.0 + " ";
			}
			
			System.out.println("\n" + row);
		}
	}
	
	public double[][] identityMatrix(int rank) {
		double[][] result = new double[rank][rank];
		for(int i=0 ; i<rank ; i++) {
			result[i][i] = 1;
		}
		return result;
	}

	public double[][] multiplyTwo(double[][] arr1, double[][] arr2){
		if(arr1 == null || arr2 == null) return null;
		if(arr1[0].length != arr2.length) return null;
		
		double[][] result = new double[arr1.length][arr2[0].length];
		for(int i=0 ; i<arr1.length ; i++) {
			for(int k=0 ; k<arr2[0].length ; k++) {
				double sum = 0;
				for(int j=0 ; j<arr1[0].length ;j++) {
					sum = sum + arr1[i][j] * arr2[j][k];
				}
				result[i][k] = sum;
			}
		}
		return result;
	}
	
	public double[][] multiply(double[][]... A) {
		double[][] result = A[0];
		for(int i=1 ; i<A.length ; i++) {
			result = multiplyTwo(result, A[i]);
		}
		return result;
	}
	
	public double[][] multiplyByConstant(double K, double[][] arr){
		double[][] result = new double[arr.length][arr[0].length];
		
		for (int i=0 ; i<arr.length ; i++) {
			for(int j=0 ; j<arr[i].length ; j++) {
				result[i][j] = arr[i][j] * K;
			}
		}
		return result;
	}
	
	public double[][] addTwo(double[][] arr1, double[][] arr2){
		if(arr1 == null || arr2 == null) return null;
		if(arr1.length != arr2.length || arr1[0].length != arr2[0].length) return null;
		
		double[][] result = new double[arr1.length][arr1[0].length];
		for(int i=0 ; i<arr1.length ; i++) {
			for(int j=0 ; j<arr1[0].length ; j++) {
				result[i][j] = arr1[i][j] + arr2[i][j];
			}
		}
		return result;
	}
	
	public double[][] add(double[][]... A) {
		double[][] result = A[0];
		for(int i=1 ; i<A.length ; i++) {
			result = addTwo(result, A[i]);
		}
		return result;
	}
	
	public double[][] subtract(double[][] arr1, double[][] arr2) {
		double[][] result = add(arr1, multiplyByConstant(-1, arr2));
		return result;
	}
	
	public double[][] transposeMatrix(double[][] arr){
		double[][] result = new double[arr[0].length][arr.length];
		for (int i=0 ; i<arr.length ; i++) {
			for(int j=0 ; j<arr[i].length ; j++) {
				result[j][i] = arr[i][j];
			}
		}
		return result;
	}
	
	public double[][] inverse(double[][] arr){
		if(arr.length != arr[0].length) return null;
		
		double[][] rref = RREF(arr); // RREF method will modify and set up matrixColumnRank and rowCombiantins
		
		if(matrixCollumnSpaceRank == arr.length) {
			return rowCombinations;
		} else {
			return null;
		}
	}
	
	public double[][] inverseFromCofactors(double[][] arr) {
		if(arr.length != arr[0].length) return null;
		if(arr.length == 1) return new double[][]{{1/arr[0][0]}};
		
		double[][] cofactorMatrix = cofactorMatrix(arr);
		double[][] transposedCofactorMatrix = transposeMatrix(cofactorMatrix);
		double[][] diagonal = multiply(arr, transposedCofactorMatrix);
		double determinant = diagonal[0][0];
		double[][] result = multiplyByConstant(1 / determinant , transposedCofactorMatrix);
		
		return result;
	}
	
	public double determinant(double[][] arr){
		if(arr.length != arr[0].length) return 0;
		
		double result = 1;
		
		double[][] upper = upperTriangularMatrix(arr);
		for(int i=0; i<upper.length; i++) {
			result = result * upper[i][i];
		}
		
		return result;
	}
	
	public double factorialDeterminant(double[][] arr) {
		if(arr.length != arr[0].length) return 0;
		int N = arr.length;
		int[][] allCombinationOfEntries = findAllCombinations(N);
		double result = 0;
		for(int i=0 ; i<allCombinationOfEntries.length ; i++) {
			double mult = 1;
			for(int j=0 ; j<allCombinationOfEntries[i].length ; j++) {
				int column = allCombinationOfEntries[i][j];
				mult = mult * arr[j][column];
			}
			if(isOddCombination(allCombinationOfEntries[i])) mult = mult * (-1);
			result = result + mult;
		}
		return result;
	}
	
	public double[][] cofactorMatrix(double[][] arr){
		if(arr.length != arr[0].length) return null;
		if(arr.length < 2) return null;
		
		double[][] result = new double[arr.length][arr[0].length];
		for(int i=0 ; i<arr.length ; i++) {
			for(int j=0 ; j<arr[i].length ; j++) {
				double[][] subMatrix = findSubMatrix(arr, i, j);
				result[i][j] = determinant(subMatrix) ;
				if(((i+j) % 2) == 1) result[i][j] = result[i][j] *(-1);
			}
		}
		return result;
	}
	
	public double[][] findSubMatrix(double[][] arr, int row, int column) {
		
		double[][] result = new double[arr.length-1][arr[0].length-1];
		for(int i=0 ; i<arr.length ; i++) {
			for(int j=0 ; j<arr[i].length ; j++) {
				if(row != i && column != j) {
					int rowEntry = i;
					int columnEntry = j;
					if(row < i) rowEntry--;
					if(column < j) columnEntry--;
					result[rowEntry][columnEntry] = arr[i][j];
				}
			}
		}
		
		return result;
	}
	
	public double[][] upperTriangularMatrix(double[][] arr){
		int rowsNumber = arr.length;
		int columnsNumber = arr[0].length;
		int iterationsNumber = 0;
		double[][] result = makeMatrixCopy(arr);
		
		rowCombinations = identityMatrix(rowsNumber);// to save matrix on which the starting matrix was multiplied from left, later this matrix will be modified as a process computes upper triangular matrix
		
		if(rowsNumber > columnsNumber) {
			iterationsNumber = columnsNumber;
		} else {
			iterationsNumber = rowsNumber-1;
		}
		
		for(int i=0 ; i<iterationsNumber ; i++) {
			for(int j=i+1 ; j<rowsNumber ; j++) {
				if(result[j][i] != 0) {
					if(result[i][i] == 0) {
						result = switchRows(result, i, j);
						rowCombinations = switchRows(rowCombinations, i, j);
					} else { 
						result = subtractRowsToGetZeroBelow(result, i, i, j);
					}
				}
			}
		}
			
		return result;
	}
	
	public double[][] lowerTriangularMatrixFromTranspose(double[][] arr){
		double[][] result = makeMatrixCopy(arr);
		
		result = transposeMatrix(result);
		result = upperTriangularMatrix(result);
		result = transposeMatrix(result);
		return result;
	}
	
	public double[][] RREF(double[][] arr){
		
		double[][] result = upperTriangularMatrixForRREF(arr);
		
		for(int i=result.length-1 ; i>0 ; i--) {
			if(isNonZeroInRow(result, i)) {
				int firstNonZeroInRow = findNonZeroCoefficient(result, i);
				for(int j=i-1 ; j>=0 ; j--) {
					result = subtractRowsToGetZeroBelow(result , i, firstNonZeroInRow, j);
				}
				
				double divisionCoefficient = result[i][firstNonZeroInRow];
				result = devideRowByCoefficient(result, i, divisionCoefficient);
		
				rowCombinations = devideRowByCoefficient(rowCombinations, i, divisionCoefficient);
				
			}
		}
		
		if(isNonZeroInRow(result, 0)) {
			double divisionCoefficient = result[0][findNonZeroCoefficient(result, 0)];
			result = devideRowByCoefficient(result, 0, divisionCoefficient);
			rowCombinations = devideRowByCoefficient(rowCombinations, 0, divisionCoefficient);
		}
		
		return result;
	}
	
	public double[][] nullSpaceBasisMatrix(double[][] arr){
		double[][] rref = RREF(arr);
		
		int nullSpaceDimension = rref[0].length;
		
		int nullSpaceRank = nullSpaceDimension - matrixCollumnSpaceRank;
		
		
		if(nullSpaceDimension == nullSpaceRank) return identityMatrix(nullSpaceRank);
		if(nullSpaceRank == 0) return new double[][]{{0}};
		
		
		double[][] nullSpaceMatrix = new double[nullSpaceDimension][nullSpaceRank];
		
		for(int j=0 ; j<nullSpaceRank ; j++) {
			int currentNPC = nonPivotColumns.get(j);// currentNonPivotColumn
			nullSpaceMatrix[currentNPC][j] = 1;
	
			for(int i=0 ; i<pivotColumns.size() && pivotColumns.get(i)<currentNPC ; i++) {
				nullSpaceMatrix[pivotColumns.get(i)][j] = -rref[i][currentNPC];
			}
		}
		
		return nullSpaceMatrix;
	}
	
	public double[][] projectionMatrix(double[][] arr) {// this method finds column space projection matrix
		double[][] A = makeMatrixCopy(arr);
		double[][] rrefA = RREF(A);
		A = removeColumns(A, nonPivotColumns);
		double[][] transposedA = transposeMatrix(A);
		double[][] inverseOfATotransposedA = inverse(multiply(transposedA, A));
		double[][] P = multiply(A, inverseOfATotransposedA, transposedA);
		return P;
	}
	
	public double[][] projection(double[][] matrixToProjectOn, double[][] toBeProjected) {// finds projection on the column space
		if(matrixToProjectOn.length != toBeProjected.length) return null;
		double[][] P =  projectionMatrix(matrixToProjectOn);
		double[][] result = multiply(P, toBeProjected);
		
		errorOfProjection = subtract(toBeProjected, result);
		
		return result;
	}
	
	public double[][] error(double[][] matrixToProjectOn, double[][] toBeProjected){
		projection(matrixToProjectOn, toBeProjected);
		return errorOfProjection;
	}
	
	//----------------------------AUXILIARY METHODS (1)----------------------------------------------------------
	
	public double[][] upperTriangularMatrixForRREF(double[][] arr){
		int rowsNumber = arr.length;
		int columnsNumber = arr[0].length;
		int iterationsNumber = 0;
		double[][] result = makeMatrixCopy(arr);
		
		
		if(rowsNumber > columnsNumber) {
			iterationsNumber = columnsNumber;
		} else {
			iterationsNumber = rowsNumber-1;
		}
		
		int nonPivotColumnsCounter = 0;
		matrixCollumnSpaceRank = 0;
		nonPivotColumns.clear();
		pivotColumns.clear();
		rowCombinations = identityMatrix(rowsNumber);// to save matrix on which the starting matrix was multiplied from left
		
		for(int i=0 ; i<iterationsNumber ; i++) {
			int zerosBelowDiagonalEntry = 0;
			
			for(int j=i - nonPivotColumnsCounter + 1 ; j<rowsNumber ; j++) {
				if(result[j][i] != 0) {
					if(result[i - nonPivotColumnsCounter][i] == 0) {
						result = switchRows(result, i - nonPivotColumnsCounter, j);
						rowCombinations = switchRows(rowCombinations, i - nonPivotColumnsCounter, j);
					} else {
						result = subtractRowsToGetZeroBelow(result, i - nonPivotColumnsCounter, i, j);
					}
				} else {
					zerosBelowDiagonalEntry++;
				}
			}
			if(result[i - nonPivotColumnsCounter][i] == 0) {
				zerosBelowDiagonalEntry++;
			} else {
				pivotColumns.add(i);
			}
			
			if(zerosBelowDiagonalEntry == rowsNumber - (i - nonPivotColumnsCounter)) {
				
				nonPivotColumns.add(i);
				nonPivotColumnsCounter++;
				if(iterationsNumber < columnsNumber) iterationsNumber++;
				
			}
		}
		
		if(columnNumberOfLastRowsPivot(result) != -1) pivotColumns.add(columnNumberOfLastRowsPivot(result));
		if(pivotColumns.size() + nonPivotColumns.size() < columnsNumber) addLastNonPivotColumns(columnsNumber);
		matrixCollumnSpaceRank = pivotColumns.size();
		
			
		return result;
	}
	
	public double[][] removeColumn(double[][] arr, int columnToRemove){
		double[][] result = new double[arr.length][arr[0].length - 1];
		
		int columnsRemoved = 0;
		
		for(int j=0 ; j<arr[0].length ; j++) {
			if(j != columnToRemove) {
				for(int i=0 ; i<arr.length ; i++) {
					result[i][j-columnsRemoved] = arr[i][j];
				}
			} else {
				columnsRemoved++;
			}
		}
		return result;
	}
	
	public double[][] addColumns(double[][] arr, HashMap<Integer, double[]> columnsList){
		double[][] result = new double[arr.length][arr[0].length + columnsList.keySet().size()];
		
		int columnsAdded = 0;
		
		for(int j=0 ; j<result[0].length ; j++) {
			if(columnsList.keySet().contains(j)) {
				for(int i=0; i<result.length ; i++) {
					double[] currentColumn = columnsList.get(j);
					result[i][j] = currentColumn[i];
				}
				columnsAdded++;
			} else {
				for(int i=0; i<result.length ; i++) {
					result[i][j] = arr[i][j-columnsAdded];
				}
			}
		}
		return result;
	}
	
	public double[][] subtractRowsToGetZeroBelow(double[][] arr, int mainRow, int column, int rowToSubtractFrom){
		double result[][] = makeMatrixCopy(arr);
		
		double substractionCoefficient = result[rowToSubtractFrom][column] / result[mainRow][column];
		for(int i=column ; i<result[0].length ; i++) {
			double input = result[rowToSubtractFrom][i] - substractionCoefficient * result[mainRow][i];
			result[rowToSubtractFrom][i] = input;
			
		}
		
		if(rowCombinations != null && arr.length == rowCombinations.length) {
			for(int i=0 ; i<rowCombinations[0].length ; i++) {
				double input = rowCombinations[rowToSubtractFrom][i] - substractionCoefficient * rowCombinations[mainRow][i];
				rowCombinations[rowToSubtractFrom][i] = input;
			}
		}
		
		return result;
	}
	
	public double[][] switchRows(double[][] arr, int row1, int row2){// this method changes initial matrix
		double[][] result = makeMatrixCopy(arr);
		
		for(int i=0; i<result[0].length ; i++) {
			double row2Input = arr[row1][i];
			double row1Input = arr[row2][i];
			result[row2][i] = row2Input;
			result[row1][i] = row1Input;
		}
		
		return result;
	}
	
	public boolean isNonZeroInRow(double[][] arr, int row) {
		boolean result = false;
		for(int i=0 ; i<arr[row].length ; i++) {
			if(arr[row][i] != 0) {
				result = true;
			}
		}
		return result;
	}
	
	public int findNonZeroCoefficient(double[][] arr, int row) {
		for(int i=0 ; i<arr[row].length ; i++) {
			if(arr[row][i] != 0) {
				return i;
			}
		}
		return -1;
	}
	
	//----------------------------AUXILIARY METHODS (2)----------------------------------------------------------
	
	private double[][] makeMatrixCopy(double[][] arr){
		double[][] result = new double[arr.length][arr[0].length];
		for(int i=0 ; i<arr.length ; i++) {
			for(int j=0 ; j<arr[i].length ; j++) {
				result[i][j] = arr[i][j];
			}
		}
		return result;
	}
	
	private int[][] findAllCombinations(int N){
		int[][] result;
		
		totalCombinations = 1;
		for(int i=1; i<=N ; i++) {
			totalCombinations = totalCombinations * i;
		}
		
		recursionCounter = 0;
		countedCombinations = 0;
		combinationsArray = new int[totalCombinations][N];
		check = new ArrayList<Integer>();
		
		
		forRecursion(N);
		result = combinationsArray;
		return result;
	}
	
	private void forRecursion(int N) { // this recursion is used to generate combinationsArray { which will contain all the combinations which can be constructed from N numbers. for example if N = 3, combiantionsArray will contain 6 different enumerated combinations (0,1,2) (0,2,1) (1,0,2) (1,2,0) (2,0,1) (2,1,0)
		for(int i=0 ; i<N ; i++) {
			if(!check.contains((Integer)i)) {
				check.add(i);
				recursionCounter++;
				if(recursionCounter < N) {
					forRecursion(N);
				} else {
					for(int k=0 ; k<N ; k++) {
						combinationsArray[countedCombinations][k] = check.get(k);
					}
					recursionCounter--;
					countedCombinations++;
				}
			}
			if(i==N-1) {
				if(recursionCounter == i)  check.remove(check.size()-1);
					
				if(check.size() != 0) check.remove(check.size()-1);
				recursionCounter--;
			}
		}
	}
	
	//-----------------------
	
	private boolean isOddCombination(int[] arr) {// true means odd
		int control = 1;
		for(int i=0 ; i<arr.length ; i++) {
			if(arr[i] != i) {
				setOnItsPlace(arr, i);
				control = control * (-1);
			}
		}
		if(control == -1) {
			return true;
		} else {
			return false;
		}
		
	}
	
	private void setOnItsPlace(int[] arr, int n) {
		int coefficient = findCoefficient(arr, n);
		int save = arr[n];
		arr[n] = arr[coefficient];
		arr[coefficient] = save;
	}
	
	private int findCoefficient(int[] arr, int n) {
		int result = -1;
		for(int i=0 ; i<arr.length ; i++) {
			if(arr[i] == n) result = i;
		}
		return result;
	}
	
	//------------------------------------------------
	
	private ArrayList<Integer> findZeroColums(double[][] arr) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		for(int i=0; i<arr[0].length ; i++) {
			if(isColumnOfZeros(arr, i)) {
				result.add(i);
			}
		}
		return result;
	}
	
	private boolean isColumnOfZeros(double[][] arr, int column) {
		for(int i=0 ; i<arr.length ; i++) {
			if(arr[i][column] != 0) return false;
		}
		return true;
	}
	
	private double[][] removeColumns(double[][] arr, ArrayList<Integer> columnsList){
		double[][] result = new double[arr.length][arr[0].length - columnsList.size()];
		int columnsSize = arr[0].length-columnsList.size();
		int columnsRemoved = 0;
		for(int j=0 ; j<arr[0].length ; j++) {
			if(!columnsList.contains(j)) {
				for(int i=0 ; i<arr.length ; i++) {
					result[i][j-columnsRemoved] = arr[i][j];
				}
			} else {
				columnsRemoved++;
			}
		}
		return result;
	}
	
	private double[][] devideRowByCoefficient(double[][] arr, int row, double coefficient){
		double[][] result = makeMatrixCopy(arr);
		for(int i=0 ; i<arr[row].length ; i++) {
			result[row][i] = result[row][i] / coefficient;
		}
		return result;
	}
	
	private void addLastNonPivotColumns(int columnsNumber) {
		int lastAddedColumn = pivotColumns.get(pivotColumns.size() - 1);
		while(lastAddedColumn < columnsNumber-1) {
			lastAddedColumn++;
			nonPivotColumns.add(lastAddedColumn);
		}
	}
	
	private int columnNumberOfLastRowsPivot(double[][] arr) {
		int rows = arr.length;
		int columns = arr[0].length;
		for(int j=0 ; j<columns ; j++) {
			if(arr[rows-1][j] != 0) return j;
		}
		return -1;
	}
	
	
	
	//instance variables
	ArrayList<Integer> check;
	int totalCombinations;
	int[][] combinationsArray;
	private int recursionCounter = 0;
	private int countedCombinations = 0;
	private int matrixCollumnSpaceRank = 0; // this will show rank of the last matrix on which findUpperTriangularMatrixForRREF was called (or a method that in itself calls this method)
	private ArrayList<Integer> nonPivotColumns = new ArrayList<Integer>(); // this ArrayList will contain all the non pivot column numbers of the last matrix on which findUpperTriangularMatrixForRREF was called (or a method that in itself calls this method)
	private ArrayList<Integer> pivotColumns = new ArrayList<Integer>(); // this ArrayList will contain all the non pivot column numbers of the last matrix on which findUpperTriangularMatrixForRREF was called (or a method that in itself calls this method)
	private double[][] errorOfProjection; // this variable is initialized when "projection" method is called (or a method that in itself calls this method)
	private double[][] rowCombinations;
	
	
}

