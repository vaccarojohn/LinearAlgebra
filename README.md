# Linear Algebra Library for C#
This package provides a library for standard linear algebra functions in C#. It provides data types and operators for matrices, vectors, linear spaces, and polynomials. It also contains methods to find the determinant, reduced row echelon form (RREF), inverse, eigenvalues/vectors, image, or kernel of a matrix. The library also provides a framework for complex number arithmetic within C# and can find complex eigenvalues and eigenvectors of a matrix.

### Example 1 (Matrix Properties)
```c#
using System;
using LinearAlgebra;

namespace Test {
    class Program {
        static void Main(string[] args) {
            Matrix myMatrix = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
            
            Console.WriteLine("Matrix Arithmetic:");
            Console.WriteLine(myMatrix + myMatrix);
            Console.WriteLine(myMatrix * myMatrix);
            
            Console.WriteLine("RREF:");
            Console.WriteLine(myMatrix.RREF());
            
            Console.WriteLine("Determinant:");
            Console.WriteLine(myMatrix.Determinant());
            
            Console.WriteLine("Trace:");
            Console.WriteLine(myMatrix.Trace());
            
            Console.WriteLine("Inverse:");
            Console.WriteLine(myMatrix.Inverse());
            
            Console.WriteLine("Transpose:");
            Console.WriteLine(myMatrix.Transpose);
            
            Console.WriteLine("Image:");
            Console.WriteLine(myMatrix.Image());
            
            Console.WriteLine("Kernel:");
            Console.WriteLine(myMatrix.Kernel());
            
            Console.WriteLine("Real Eigenvalues:");
            foreach (double eigenvalue in myMatrix.GetRealEigenvalues()) {
                Console.WriteLine(eigenvalue);
            }
            
            Console.WriteLine("All Eigenvalues:");
            foreach (ComplexNumber eigenvalue in myMatrix.GetComplexEigenvalues()) {
                Console.WriteLine(eigenvalue);
            }
            
            Console.WriteLine("Diagonalization: " + myMatrix.Diagonalize());
            
            //If your matrix has exactly two imaginary eigenvalues --- Console.WriteLine("Diagonalization: " + myMatrix.DiagonalizeToRotationMatrix());
        }
    }
}
```
