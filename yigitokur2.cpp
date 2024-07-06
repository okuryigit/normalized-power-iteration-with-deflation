#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include<string>
#include <math.h>
#include <random>

using namespace std;

class Matrix {
private:
	int columns, rows;//member vatiables of Matrix class are defined.
	double** matarray;
public:
	Matrix(int row_size, int col_size) : rows(row_size), columns(col_size) { //contstructor
		matarray = new double* [rows];
		for (int i = 0; i < rows; i++) {
			matarray[i] = new double[columns];
		}
	}
	~Matrix() {
		for (int i = 0; i < rows; i++) {
			delete[] matarray[i]; // Delete each sub-array
		}
		delete[] matarray; // Delete the main array
	}

	Matrix(const Matrix& other) : rows(other.rows), columns(other.columns) {//copy constructor helps us to work with dynamically
		matarray = new double* [rows];                                      //allocated memory properly.
		for (int i = 0; i < rows; i++) {
			matarray[i] = new double[columns];
			for (int j = 0; j < columns; j++) {
				matarray[i][j] = other.matarray[i][j]; // each element is copied.
			}
		}
	}

	Matrix& operator=(const Matrix& mat_equal) {
		if (this != &mat_equal) {// delete the memory if they are not same.
			for (int i = 0; i < rows; i++) {
				delete[] matarray[i];
			}
			delete[] matarray;

			rows = mat_equal.rows;
			columns = mat_equal.columns;
			matarray = new double* [rows];
			for (int i = 0; i < rows; i++) {// dynamic memory of the new equaled matrix allocated.
				matarray[i] = new double[columns];
				for (int j = 0; j < columns; j++) {
					matarray[i][j] = mat_equal.matarray[i][j]; // each elementis copied to allocate dynamically.
				}
			}
		}
		return *this;
	}

	double get_value(int rows, int columns) const {   // returns the desired element of the matrix without changing anyhting.
		return matarray[rows][columns];
	}

	double set_value(int rows, int columns, double value) {  //setting a value for a speficic element of the matrix.
		return matarray[rows][columns] = value;
	}

	void readFromFile(const std::string& name_f) { // member function that reads the input matrix file and fills the values of the object. 
		ifstream file(name_f);
		if (!file.is_open()) {
			cout << "Error opening file: " << name_f << endl;
			return;
		}
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				file >> matarray[i][j];
			}
		}
		file.close();
	}
	Matrix operator*(const double& factor) const {//overloading function for matrix and scalar multiplication.
		Matrix multiplied(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {

				multiplied.set_value(i, j, matarray[i][j] * factor);
			}
		}
		return multiplied;
	}

	Matrix operator*(const Matrix& mat2) const {//overloading function for matrix multiplication.
		Matrix multiplied(rows, mat2.columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < mat2.columns; j++) {
				double sum = 0;
				for (int k = 0; k < columns; k++) {
					sum += matarray[i][k] * mat2.matarray[k][j];
				}
				multiplied.set_value(i, j, sum);
			}
		}
		return multiplied;
	}
	Matrix operator-(const Matrix& mat2) const {//overloading function that subtracts matrix in the argument from the object matrix. 
		Matrix subtracted(rows, columns);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < mat2.columns; j++) {
				double result = 0;
				result = matarray[i][j] - mat2.matarray[i][j];
				subtracted.set_value(i, j, result);
			}
		}
		return subtracted;
	}

	Matrix transpose() const {//function that takes the transpose of the matrix object.
		Matrix tmat(columns, rows);
		for (int i = 0; i < rows; i++) {
			for (int k = 0; k < columns; k++) {
				tmat.set_value(k, i, matarray[i][k]);
			}
		}
		return tmat;
	}

	double dot_product(const Matrix& v) const {//dot product of the two vectors.

		double result = 0.0;                     // needed in the deflation part.
		for (int i = 0; i < rows; ++i) {
			result += (*this).matarray[i][0] * v.matarray[i][0];// multiplying corresponding elements and adding them all.
		}
		return result;
	}

	double two_norm() const { //member function to find two norm of a vector.
		double norm_squared = 0.0;
		for (int i = 0; i < rows; i++) {
			norm_squared += matarray[i][0] * matarray[i][0];//adding the squares of the element in each row.
		}
		return sqrt(norm_squared);//taking the square of the summation obtained above.
	}

	double inf_norm()const {  //infinite norm of a matrix means absolute max row sum.
		double max_row_sum = 0;
		for (int i = 0; i < rows; i++) {
			double temp_sum = 0;
			for (int j = 0; j < columns; j++) {
				temp_sum += abs(matarray[i][j]); //absolute row sum found.
			}
			if (temp_sum >= max_row_sum) { max_row_sum = temp_sum; }
		}
		return max_row_sum;

	}

	void normalizer() {//normalizing the matrix object using previous infinite norm function.
		double norm = inf_norm();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matarray[i][j] /= norm;
			}
		}
	}
	void power_iteration(const double tol, double& eigenvalue, Matrix& eigenvector) {// member function that performs power iteration.
		Matrix vector_x(rows, 1);
		double diff = tol + 1;// just to ensure it is higher than the tolerance at the start.
		for (int i = 0; i < rows; i++) {// initializing vector_x with random values between 0 and 1
			vector_x.set_value(i, 0, static_cast<double>(rand()) / RAND_MAX);//random value found converted to double and divided to max
		}                                                                    //to ensure it is between 0 and 1.     

		int max_iterations = 1000; // limit for max iteration.
		int iterations = 0; //current value for iteration
		double pre_eigen = 0.0;  //eigenvalue initialized.

		while (diff > tol && iterations < max_iterations) {//condtions for tolerance and max iterations.
			Matrix Ax = (*this) * vector_x;  // calculating the next vector using matrix-vector multiplication using the formula.
			double new_eigen = Ax.inf_norm();// infinite norm of the product that is expected to converge to the eigenvalue after iterations.
			Ax.normalizer();   //normalized the product to use in the next iteration as a vector_x.
			if (Ax.dot_product(vector_x) < 0) {
				new_eigen *= -1;//considering negative dominant eigenvalue, checking its sign 
				// and correcting if necessary.
				Ax = Ax * (-1);// adjust also the sign of the eigenvector accordingly.
			}
			vector_x = Ax;   //updating the vector_x with normalized Ax.                             
			diff = fabs(new_eigen - pre_eigen);// checking the difference of consequent eigenvalues
			//if it is less than tolerance, while lopp will end.
			pre_eigen = new_eigen;// updating the eigenvalue for the next iteration.

			iterations++; // increasing iterration number.
		}
		eigenvalue = pre_eigen;//resulting eigenvalue and eigenvector are assigned to the argument to use in the main.
		eigenvector = vector_x;
	}

	Matrix deflation_eigen2(Matrix& eigenvector) {//the deflation part to find second eigenvalue, argument is eigenvector found after iterations.
		Matrix identity(rows, columns);           //identity matrix is constructed.
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (i == j) {
					identity.set_value(i, j, 1);
				}
				else identity.set_value(i, j, 0);
			}
		}

		Matrix v = eigenvector;//v is initialized as eigenvector which is found from power iteration func.
		double norm2 = v.two_norm();//v is modified to perform H=I-2vvT/vTv householder transformation.
		if (v.get_value(0, 0) < 0) {
			v.set_value(0, 0, (v.get_value(0, 0) - norm2));
		}
		else v.set_value(0, 0, v.get_value(0, 0) + norm2);

		Matrix b = v * v.transpose();//in order to calculate H=I-2vvT/vTv, matrix multiplication and dot product calculated.
		double value = v.dot_product(v);
		double val2 = 2.0 / value;
		Matrix H = identity - (b * val2);
		Matrix G = (H * (*this)) * H.transpose();//since (H A H^-1 = H A H^T^) holds,matrix G calculated.
		Matrix A_new(rows - 1, columns - 1);//A_new (n-1 x n-1) constructed by eliminating first row and column of G.
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < columns; j++) {
				A_new.set_value(i - 1, j - 1, G.get_value(i, j));
			}
		}
		return A_new;
	}
};

int main(int argc, char* argv[]) {
	// Control the input number.
	if (argc != 4) {  //argument count argc must be 4 since argv[0] will always contain the name of the program.
		cout << "Wrong number of inputs,please check your inputs.";
		return 0;
	}

	ifstream mfile;
	mfile.open(argv[1]);
	string line;
	int n = 0;
	if (mfile.is_open()) {// finding the column and row number (n) of the input matrix A.
		while (getline(mfile, line)) {
			n = n + 1;
		}
	}
	else cout << "Unable to open the file";
	mfile.clear();
	mfile.seekg(0);//resetting and closing the file.
	mfile.close();

	double tolerance = atof(argv[2]);// tolerance value inserted as string in command line ,converted to double.
	Matrix A(n, n);   // nxn matrix object A constructed.
	double eigenvalue; // eigenvalue is defined.
	Matrix eigenvector(n, 1);//nx1 matrix object eigenvector constructed.
	A.readFromFile(argv[1]);// matrix A is filled with the values read from A.txt file.
	A.power_iteration(tolerance, eigenvalue, eigenvector);// power iteration performed to the matrix A.
	
	ofstream ofile;
	ofile.open(argv[3]);     //writing at output file eigenvalue 1 and eigenvector.
	ofile << "Eigenvalue#1: " << eigenvalue << endl;
	for (int i = 0; i < n; i++) {
		ofile << eigenvector.get_value(i, 0) << endl;
	}

	Matrix A2 = A.deflation_eigen2(eigenvector);   //constructed A2 is deflated using the member function.
	double eigenvalue_2;   //new eigenvalue 2 defined to use in power iteration.
	Matrix eigenvector2(n, 1);   //new eigenvector2 is constructed to use in power iteration.
	A2.power_iteration(tolerance, eigenvalue_2, eigenvector2);  //second eigenvalue is found by power iteration of deflated matrix.

	ofile << "Eigenvalue#2: " << eigenvalue_2 << endl;  //second eigenvalue is written on output file.
	ofile.close(); //output file  is closed.
	return 0;
}