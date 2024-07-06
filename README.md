# normalized power iteration with deflation
 Computing Eigenvalues and Eigenvectors using Normalized Power Iteration together with Deflation 
This program is coded in order to implement the normalized power iteration and deflation algorithms to find the dominant eigenvalue of the input matrix A, eigenvector and the second eigenvalue. This program takes 3 command line arguments, input matrix A.txt, tolerance value and output.txt file. Output contains dominant eigenvalue, eigenvector and second eigenvalue. Input (nxn) Matrix A is put in a normalized power iteration algorithm, dominant eigenvalue and its eigenvector found. Then using eigenvector and householder transformation, Matrix A is deflated. Normalized power iteration takes place again for the deflated (n-1 x n-1) matrix and second dominant eigenvalue found. During the code implementation, Object oriented programming (OOP) is used in order to avoid complex loops. Operations such as normalized power iteration, deflation, dot product, matrix multiplication and subtraction, normalizing a matrix and transpose are defined effectively with Matrix class member functions.
 
 
In order to compile the .cpp file in terminal, run the following:


g++ -o yigitokur2 yigitokur2.cpp


In order to run the compiled program with the specified command line arguments, run the following:


./deneme A.txt 1e-6 output.txt


