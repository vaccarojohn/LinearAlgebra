using System;
using System.Linq;
using System.Collections.Generic;

namespace LinearAlgebra
{
    /// <summary>
    /// A two-dimensional array representing an m x n matrix
    /// </summary>
    public class Matrix
    {
        /// <summary>
        /// The method of getting or changing a value at a specified position of the matrix
        /// </summary>
        /// <param name="row">The row number of the value to get or set (0 inclusive)</param>
        /// <param name="col">The column number of the value to get or set (0 inclusive)</param>
        public double this[int row, int col] { get { return rElements[row, col]; } set { rElements[row, col] = value; } }

        /// <summary>
        /// Returns the elements of the matrix as a two-dimensional double array
        /// </summary>
        public double[,] Elements { get { return rElements; } }

        /// <summary>
        /// The number of rows in the matrix
        /// </summary>
        public int Rows { get { return m; } }

        /// <summary>
        /// The number of columns in the matrix
        /// </summary>
        public int Cols { get { return n; } }

        /// <summary>
        /// Returns the rank of the matrix, or the dimension of its image
        /// </summary>
        public int Rank { get { return Image().Dimension; } }

        /// <summary>
        /// Returns the nullity of the matrix, or the dimension of its kernel
        /// </summary>
        public int Nullity { get { return Kernel().Dimension; } }

        /// <summary>
        /// Returns the transpose of the matrix
        /// </summary>
        public Matrix Transpose { get { return new Matrix(rElements, false); } }

        /// <summary>
        /// A boolean that is true when the matrix is square (i.e. the number of rows is equivalent to the number of columns)
        /// </summary>
        public bool IsSquare { get { return m == n; } }

        /// <summary>
        /// A boolean that is true when the matrix is an identity matrix (Note: The matrix must be square in order for this value to be true.)
        /// </summary>
        public bool IsIdentity
        {
            get
            {
                if (n != m)
                {
                    return false;
                }

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if ((i == j && rElements[i, j] != 1) || (i != j && rElements[i, j] != 0))
                        {
                            return false;
                        }
                    }
                }

                return true;
            }
        }

        /// <summary>
        /// A boolean that is true when the matrix can be inverted
        /// </summary>
        public bool IsInvertible { get { return (Determinant() != 0); } }

        /// <summary>
        /// A boolean that is true when the matrix is orthogonal (i.e. when the matrix is square and its transformation preserves magnitude of vectors)
        /// </summary>
        public bool IsOrthogonal
        {
            get
            {
                Space im = Image();
                return ((m == n) && im.Dimension == m && im.HasOrthonormalBasis);
            }
        }

        /// <summary>
        /// A boolean that is true when the matrix is symmetric (i.e. the matrix is equivalent to its transpose)
        /// </summary>
        public bool IsSymmetric
        {
            get
            {
                return (this == Transpose);
            }
        }

        /// <summary>
        /// A boolean that is true when the matrix is skew-symmetric (i.e. the matrix is equivalent to its transpose multiplied by -1)
        /// </summary>
        public bool IsSkewSymmetric
        {
            get
            {
                return (this == -Transpose);
            }
        }

        /// <summary>
        /// A boolean that is true when the matrix is diagonal (i.e. any non-zero values in the matrix are located in diagonal positions)
        /// </summary>
        public bool IsDiagonal
        {
            get
            {
                if (m != n)
                {
                    return false;
                }

                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j && rElements[i, j] != 0)
                        {
                            return false;
                        }
                    }
                }

                return true;
            }
        }

        /// <summary>
        /// A boolean that is true if the matrix is similar to a diagonal matrix
        /// </summary>
        public bool IsDiagonalizable
        {
            get
            {
                if (n != m)
                {
                    return false;
                }

                Eigenpair[] eigenpairs = GetRealEigenpairs();
                int multiplicity = 0;
                for (int i = 0; i < eigenpairs.Length; i++)
                {
                    multiplicity += eigenpairs[i].gemu;
                }

                if (multiplicity != m)
                {
                    return false;
                }

                return true;
            }
        }

        /// <summary>
        /// A boolean that is true if the matrix is similar to a rotation matrix
        /// </summary>
        public bool IsDiagonalizableToRotationMatrix
        {
            get
            {
                if (n != m)
                {
                    return false;
                }

                ComplexEigenpair[] eigenpairs = GetComplexEigenpairs();
                int multiplicityReal = 0;
                int multiplicityComplex = 0;
                for (int i = 0; i < eigenpairs.Length; i++)
                {
                    if (eigenpairs[i].eigenvalue.IsReal)
                    {
                        multiplicityReal += eigenpairs[i].eigenvectors.Length;
                    }
                    else
                    {
                        multiplicityComplex += eigenpairs[i].eigenvectors.Length;
                    }
                }

                if ((multiplicityReal + multiplicityComplex) != m || multiplicityComplex != 2)
                {
                    return false;
                }

                return true;
            }
        }

        private int m;
        private int n;
        private double[,] rElements;

        /// <summary>
        /// Initializes a zero square matrix in a specified number of dimensions
        /// </summary>
        /// <param name="dim">The dimension of the matrix</param>
        public static Matrix Zero(int dim)
        {
            return new Matrix(dim, dim);
        }

        /// <summary>
        /// Initializes a matrix representing the identity matrix in a specified number of dimensions
        /// </summary>
        /// <param name="dim">The dimension of the matrix</param>
        public static Matrix Identity(int dim)
        {
            double[] diagonals = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                diagonals[i] = 1;
            }

            return new Matrix(diagonals);
        }

        /// <summary>
        /// Initializes a matrix representing a scaling matrix in a specified number of dimensions
        /// </summary>
        /// <param name="dim">The dimension of the matrix</param>
        /// <param name="scaleFactor">The scale factor of the matrix</param>
        /// <returns></returns>
        public static Matrix Scale(int dim, int scaleFactor)
        {
            double[] diagonals = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                diagonals[i] = scaleFactor;
            }

            return new Matrix(diagonals);
        }

        /// <summary>
        /// Initializes a matrix representing orthogonal projection onto a subspace
        /// </summary>
        /// <param name="subspace">The subspace onto which vectors can be projected</param>
        public static Matrix Projection(Space subspace)
        {
            if (subspace.Dimension == 0)
            {
                throw new LinearAlgebraException("Cannot project onto the zero subspace.");
            }

            Space newSpace = subspace;
            newSpace.Orthonormalize();

            return (newSpace.ChangeOfBasisMatrix * newSpace.ChangeOfBasisMatrix.Transpose);
        }

        /// <summary>
        /// Initializes a matrix representing reflection across a subspace
        /// </summary>
        /// <param name="subspace">The subspace over which vectors can be reflected</param>
        public static Matrix Reflection(Space subspace)
        {
            if (subspace.Dimension == 0)
            {
                throw new LinearAlgebraException("Cannot reflect across the zero subspace.");
            }

            double[] diagonals = new double[subspace.ParentDimension];
            Vector[] newBasis = new Vector[subspace.ParentDimension];
            int i = 0;
            while (i < subspace.Dimension)
            {
                newBasis[i] = subspace.Basis[i];
                diagonals[i] = 1;
                i++;
            }

            Space perpSpace = subspace.PerpendicularSpace();
            while (i < subspace.ParentDimension)
            {
                newBasis[i] = perpSpace.Basis[i - subspace.Dimension];
                diagonals[i] = -1;
                i++;
            }

            Space newSpace = new Space(newBasis);
            return (newSpace.ChangeOfBasisMatrix * new Matrix(diagonals) * newSpace.ChangeOfBasisMatrix.Inverse());
        }

        /// <summary>
        /// Initializes a matrix representing counterclockwise rotation about a subspace
        /// </summary>
        /// <param name="subspace">The subspace about which vectors can be rotated (Note: The dimension of the subspace should be 2 less than the dimension of vectors in the subspace.)</param>
        /// <param name="angle">The angle of rotation in radians</param>
        public static Matrix Rotation(Space subspace, double angle)
        {
            if (subspace.ParentDimension != (subspace.Dimension + 2))
            {
                throw new LinearAlgebraException("The parent dimension of the subspace must be exactly " + (subspace.Dimension + 2) + " in order to get the rotation matrix.");
            }

            double[,] rElements = new double[subspace.ParentDimension, subspace.ParentDimension];
            Vector[] newBasis = new Vector[subspace.ParentDimension];
            int i = 0;
            while (i < subspace.Dimension)
            {
                newBasis[i] = subspace.Basis[i];
                for (int j = 0; j < subspace.ParentDimension; j++)
                {
                    if (i == j)
                    {
                        rElements[i, j] = 1;
                    }
                    else
                    {
                        rElements[i, j] = 0;
                    }
                }
                i++;
            }

            Space perpSpace = subspace.PerpendicularSpace();
            newBasis[i] = perpSpace.Basis[i - subspace.Dimension];
            newBasis[i + 1] = perpSpace.Basis[i + 1 - subspace.Dimension];
            rElements[i, i] = (Math.Cos(angle));
            rElements[i, i + 1] = (-Math.Sin(angle));
            rElements[i + 1, i] = (Math.Sin(angle));
            rElements[i + 1, i + 1] = (Math.Cos(angle));

            Space newSpace = new Space(newBasis);
            newSpace.Orthonormalize();
            return (newSpace.ChangeOfBasisMatrix * new Matrix(rElements) * newSpace.ChangeOfBasisMatrix.Inverse());
        }

        /// <summary>
        /// Initializes a matrix representing a shear
        /// </summary>
        /// <param name="dim">The dimension of the matrix</param>
        /// <param name="row">The row number of the shear value (0 inclusive)</param>
        /// <param name="col">The column number of the shear value (0 inclusive)</param>
        /// <param name="value">The value for the shear</param>
        public static Matrix Shear(int dim, int row, int col, double value)
        {
            Matrix shear = Identity(dim);
            shear[row, col] = value;
            return shear;
        }

        /// <summary>
        /// Rounds all values of a matrix to a specified number of digits
        /// </summary>
        /// <param name="a">The matrix to round</param>
        /// <param name="digits">The number of digits</param>
        public static Matrix Round(Matrix a, int digits)
        {
            double[,] rElements = a.rElements;
            for (int i = 0; i < a.m; i++)
            {
                for (int j = 0; j < a.n; j++)
                {
                    rElements[i, j] = Math.Round(rElements[i, j], digits);
                }
            }

            return new Matrix(rElements);
        }

        /// <summary>
        /// Flips the sign of all values in a matrix
        /// </summary>
        public static Matrix operator -(Matrix a)
        {
            return (-1 * a);
        }

        /// <summary>
        /// Returns a matrix representing the sum of two matrices
        /// </summary>
        public static Matrix operator +(Matrix a, Matrix b)
        {
            if (a.Rows == b.Rows && a.Cols == b.Cols)
            {
                double[,] rElements = new double[a.Rows, a.Cols];
                for (int i = 0; i < a.Rows; i++)
                {
                    for (int j = 0; j < a.Cols; j++)
                    {
                        rElements[i, j] = a[i, j] + b[i, j];
                    }
                }

                return new Matrix(rElements, true);
            }
            else
            {
                throw new LinearAlgebraException("Matrices must have the same number of rows and columns.");
            }
        }

        /// <summary>
        /// Returns a matrix representing the difference of two matrices
        /// </summary>
        public static Matrix operator -(Matrix a, Matrix b)
        {
            if (a.Rows == b.Rows && a.Cols == b.Cols)
            {
                double[,] rElements = new double[a.Rows, a.Cols];
                for (int i = 0; i < a.Rows; i++)
                {
                    for (int j = 0; j < a.Cols; j++)
                    {
                        rElements[i, j] = a[i, j] - b[i, j];
                    }
                }

                return new Matrix(rElements, true);
            }
            else
            {
                throw new LinearAlgebraException("Matrices must have the same number of rows and columns.");
            }
        }

        /// <summary>
        /// Returns a matrix representing the multiplication of a scalar and a matrix
        /// </summary>
        public static Matrix operator *(double a, Matrix b)
        {
            double[,] rElements = new double[b.Rows, b.Cols];
            for (int i = 0; i < b.Rows; i++)
            {
                for (int j = 0; j < b.Cols; j++)
                {
                    rElements[i, j] = a * b[i, j];
                }
            }

            return new Matrix(rElements, true);
        }

        /// <summary>
        /// Returns a matrix representing the multiplication of a scalar and a matrix
        /// </summary>
        public static Matrix operator *(Matrix a, double b)
        {
            double[,] rElements = new double[a.Rows, a.Cols];
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Cols; j++)
                {
                    rElements[i, j] = b * a[i, j];
                }
            }

            return new Matrix(rElements, true);
        }

        /// <summary>
        /// Returns a matrix representing the multiplication of a matrix and a vector (Note: The dimension of the vector must be equivalent to the number of columns in the matrix.)
        /// </summary>
        public static Vector operator *(Matrix a, Vector b)
        {
            if (a.Cols == b.Dimension)
            {
                double[] cElements = new double[a.Rows];
                for (int i = 0; i < a.Rows; i++)
                {
                    double el = 0;
                    for (int j = 0; j < a.Cols; j++)
                    {
                        el += a[i, j] * b[j];
                    }
                    cElements[i] = el;
                }

                return new Vector(cElements);
            }
            else
            {
                throw new LinearAlgebraException("The matrix must have a number of columns equivalent to the dimension of the vector.");
            }
        }

        /// <summary>
        /// Returns a matrix representing the multiplication of two matrices (Note: The first matrix must have a number of columns equivalent to the number of rows in the second matrix.)
        /// </summary>
        public static Matrix operator *(Matrix a, Matrix b)
        {
            if (a.Cols == b.Rows)
            {
                double[,] rElements = new double[a.Rows, b.Cols];
                for (int i = 0; i < a.Rows; i++)
                {
                    for (int j = 0; j < b.Cols; j++)
                    {
                        double el = 0;
                        for (int k = 0; k < a.Cols; k++)
                        {
                            el += a[i, k] * b[k, j];
                        }

                        rElements[i, j] = el;
                    }
                }

                return new Matrix(rElements, true);
            }
            else
            {
                throw new LinearAlgebraException("The second matrix must have a number of rows equivalent to the number of columns in the first matrix.");
            }
        }

        /// <summary>
        /// Returns a matrix representing the quotient of a matrix and a scalar
        /// </summary>
        public static Matrix operator /(Matrix a, double b)
        {
            if (b == 0)
            {
                throw new DivideByZeroException();
            }

            double[,] rElements = new double[a.Rows, a.Cols];
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Cols; j++)
                {
                    rElements[i, j] = (1 / b) * a[i, j];
                }
            }

            return new Matrix(rElements, true);
        }

        /// <summary>
        /// Determines whether two matrices are equal (Note: Matrices must have the same number of rows and columns and the same value at every position to be considered equal.)
        /// </summary>
        public static bool operator ==(Matrix a, Matrix b)
        {
            return a.Equals(b);
        }

        public static bool operator !=(Matrix a, Matrix b)
        {
            return !a.Equals(b);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Matrix b = (Matrix)obj;
            if (m != b.Rows || n != b.Cols)
            {
                return false;
            }

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (rElements[i, j] != b[i, j])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Initializes an empty matrix with zero values
        /// </summary>
        /// <param name="m">The number of rows in the matrix</param>
        /// <param name="n">The number of columns in the matrix</param>
        public Matrix(int m, int n)
        {
            this.m = m;
            this.n = n;
            rElements = new double[m, n];
        }

        /// <summary>
        /// Initializes a matrix with the elements of a two-dimensional array
        /// </summary>
        /// <param name="elements">The elements of the matrix</param>
        /// <param name="rows">Whether the elements in the array are formatted as rows or columns (e.g. { { first row }, { second row }, ... } vs. { { first col }, { second col }, ... } })</param>
        public Matrix(double[,] elements, bool rows = true)
        {
            if (rows)
            {
                m = elements.GetLength(0);
                n = elements.GetLength(1);
                rElements = elements;
            }
            else
            {
                m = elements.GetLength(1);
                n = elements.GetLength(0);
                rElements = new double[m, n];
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        rElements[i, j] = elements[j, i];
                    }
                }
            }
        }

        /// <summary>
        /// Initializes a matrix with values in an array on the diagonals
        /// </summary>
        /// <param name="diagonals">The diagonals of the matrix</param>
        public Matrix(double[] diagonals)
        {
            m = diagonals.Length;
            n = diagonals.Length;
            rElements = new double[m, n];
            for (int i = 0; i < diagonals.Length; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    rElements[i, j] = 0;
                }

                rElements[i, i] = diagonals[i];
                for (int j = (i + 1); j < n; j++)
                {
                    rElements[i, j] = 0;
                }
            }
        }

        /// <summary>
        /// Adds a row to the matrix
        /// </summary>
        /// <param name="row">The row number of the new row (0 inclusive)</param>
        /// <param name="elements">The elements of the new row</param>
        public void AddRow(int row, double[] elements)
        {
            if (elements.Length != n)
            {
                throw new LinearAlgebraException("The number of items in the row should be " + n + ".");
            }

            double[,] rElements = new double[m + 1, n];
            for (int i = 0; i <= m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == row)
                    {
                        rElements[i, j] = elements[j];
                    }
                    else if (i < row)
                    {
                        rElements[i, j] = this.rElements[i, j];
                    }
                    else
                    {
                        rElements[i, j] = this.rElements[i - 1, j];
                    }
                }
            }

            this.rElements = rElements;
            m += 1;
        }

        /// <summary>
        /// Adds a column to the matrix
        /// </summary>
        /// <param name="col">The column number of the new column (0 inclusive)</param>
        /// <param name="elements">The elements of the new column</param>
        public void AddCol(int col, double[] elements)
        {
            if (elements.Length != m)
            {
                throw new LinearAlgebraException("The number of items in the column should be " + m + ".");
            }

            double[,] rElements = new double[m, n + 1];
            for (int j = 0; j <= n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    if (j == col)
                    {
                        rElements[i, j] = elements[i];
                    }
                    else if (j < col)
                    {
                        rElements[i, j] = this.rElements[i, j];
                    }
                    else
                    {
                        rElements[i, j] = this.rElements[i, j - 1];
                    }
                }
            }

            this.rElements = rElements;
            n += 1;
        }

        /// <summary>
        /// Removes a row from the matrix
        /// </summary>
        /// <param name="row">The row number of the row to remove (0 inclusive)</param>
        public void RemoveRow(int row)
        {
            double[,] rElements = new double[m - 1, n];
            for (int i = 0; i < m; i++)
            {
                if (i == row)
                {
                    continue;
                }
                for (int j = 0; j < n; j++)
                {
                    if (i > row)
                    {
                        rElements[i - 1, j] = this.rElements[i, j];
                    }
                    else
                    {
                        rElements[i, j] = this.rElements[i, j];
                    }
                }
            }

            this.rElements = rElements;
            m -= 1;
        }

        /// <summary>
        /// Removes a column from the matrix
        /// </summary>
        /// <param name="col">The column number of the column to remove (0 inclusive)</param>
        public void RemoveCol(int col)
        {
            double[,] rElements = new double[m, n - 1];
            for (int j = 0; j < n; j++)
            {
                if (j == col)
                {
                    continue;
                }
                for (int i = 0; i < m; i++)
                {
                    if (j > col)
                    {
                        rElements[i, j - 1] = this.rElements[i, j];
                    }
                    else
                    {
                        rElements[i, j] = this.rElements[i, j];
                    }
                }
            }

            this.rElements = rElements;
            n -= 1;
        }

        /// <summary>
        /// Returns the trace of a square matrix, or the sum of its diagonals
        /// </summary>
        public double Trace()
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The trace only exists for square matrices.");
            }

            double tr = 0;
            for (int i = 0; i < m; i++)
            {
                tr += rElements[i, i];
            }

            return Math.Round(tr, 6);
        }

        /// <summary>
        /// Returns the determinant of a square matrix (computed using row-reduction)
        /// </summary>
        public double Determinant()
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The determinant only exists for square matrices.");
            }

            RREF(out double det);
            return Math.Round(det, 6);
        }

        /// <summary>
        /// Returns the inverse of the matrix
        /// </summary>
        public Matrix Inverse()
        {
            if (m != n)
            {
                throw new LinearAlgebraException("Only square matrices can have inverses.");
            }

            if (Determinant() == 0)
            {
                throw new LinearAlgebraException("The matrix is not invertible.");
            }

            int dim = m;
            Matrix myMatrix = new Matrix(rElements);
            for (int i = 0; i < dim; i++)
            {
                myMatrix.AddCol(myMatrix.Cols, Vector.StandardUnit(i, dim).ToArray());
            }

            myMatrix = myMatrix.RREF();
            for (int i = 1; i <= dim; i++)
            {
                myMatrix.RemoveCol(0);
            }

            return myMatrix;
        }

        /// <summary>
        /// Returns a new matrix representing the reduced row echelon form of this matrix
        /// </summary>
        public Matrix RREF()
        {
            double[,] rElements = (double[,])this.rElements.Clone();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (Math.Round(rElements[i, j], (-(int)Math.Sqrt(n * m) + 9).Clamp(0, 3)) == 0 && j >= i)
                    {
                        rElements[i, j] = 0;
                        for (int k = (i + 1); k < m; k++)
                        {
                            if (rElements[k, j] != 0)
                            {
                                for (int l = 0; l < n; l++)
                                {
                                    double placeholder = rElements[k, l];
                                    rElements[k, l] = rElements[i, l];
                                    rElements[i, l] = placeholder;
                                }

                                j = -1;
                                break;
                            }
                        }
                    }
                    else if (rElements[i, j] != 0)
                    {
                        bool zeroPossible = false;
                        for (int k = i - 1; k >= 0; k--)
                        {
                            if (rElements[k, j] == 1)
                            {
                                int l = 0;
                                while (rElements[k, l] == 0) l++;

                                if (l == j)
                                {
                                    double factor = rElements[i, j];
                                    for (int h = 0; h < n; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[k, h];
                                    }

                                    zeroPossible = true;
                                    break;
                                }
                            }
                        }
                        if (!zeroPossible)
                        {
                            double factor = rElements[i, j];
                            for (int k = 0; k < n; k++)
                            {
                                rElements[i, k] /= factor;
                            }
                            break;
                        }
                    }
                }
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    if (rElements[i, j] != 0)
                    {
                        double factor = rElements[i, j];
                        for (int l = Math.Min(j, m - 1); l > i; l--)
                        {
                            if (rElements[l, j] == 1)
                            {
                                int k = 0;
                                while (rElements[l, k] == 0) k++;

                                if (k == j)
                                {
                                    for (int h = 0; h < n; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[l, h];
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            Matrix myMatrix = new Matrix(rElements);
            return Round(myMatrix, 6);
        }

        /// <summary>
        /// Returns the image of the matrix
        /// </summary>
        public Space Image()
        {
            return new Space(ToVectorArray());
        }

        /// <summary>
        /// Returns the kernel of the matrix
        /// </summary>
        public Space Kernel()
        {
            Matrix rref = RREF();
            List<Vector> kernelVectors = new List<Vector>();
            bool[] columnsNormalized = new bool[rref.Cols];
            for (int j = 0; j < rref.Cols; j++)
            {
                bool allZero = true;
                for (int i = 0; i < rref.Rows; i++)
                {
                    if (rref[i, j] != 0)
                    {
                        allZero = false;
                        break;
                    }
                }

                if (allZero)
                {
                    columnsNormalized[j] = true;
                    kernelVectors.Add(Vector.StandardUnit(j, rref.Cols));
                }
            }
            for (int i = 0; i < rref.Rows; i++)
            {
                int leadingOneCol = -1;
                for (int j = 0; j < rref.Cols; j++)
                {
                    Vector kernelVec = Vector.Zero(rref.Cols);
                    if (rref[i, j] == 1 && leadingOneCol == -1)
                    {
                        leadingOneCol = j;
                    }
                    else if (Math.Round(rref[i, j], 5) != 0 && leadingOneCol != -1 && !columnsNormalized[j])
                    {
                        kernelVec[leadingOneCol] = Math.Round(rref[i, j], 6);
                        kernelVec[j] = -1;
                        for (int k = (i + 1); k < rref.Rows; k++)
                        {
                            if (rref[k, j] != 0)
                            {
                                int l = 0;
                                while (rref[k, l] != 1) l++;
                                kernelVec[l] = Math.Round(rref[k, j], 5);
                            }
                        }

                        columnsNormalized[j] = true;
                        kernelVectors.Add(kernelVec);
                    }
                }
            }

            if (kernelVectors.Count == 0)
            {
                return new Space(rref.Cols);
            }
            else
            {
                return new Space(kernelVectors.ToArray());
            }
        }

        /// <summary>
        /// Returns a matrix similar to this matrix with respect to a different basis
        /// </summary>
        /// <param name="basis">A space representing the change of basis matrix</param>
        public Matrix GetBasisMatrix(Space basis)
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The matrix must be square to transform to basis coordinates.");
            }
            else if (m != basis.Dimension)
            {
                throw new LinearAlgebraException("The basis must have exactly " + m + " linearly independent vectors.");
            }
            else if (basis.Dimension != basis.ParentDimension)
            {
                throw new LinearAlgebraException("The basis dimension and parent dimension must align.");
            }

            return (basis.ChangeOfBasisMatrix.Inverse() * this * basis.ChangeOfBasisMatrix);
        }

        /// <summary>
        /// Returns the characteristic polynomial of the matrix
        /// </summary>
        public Polynomial CharacteristicPolynomial()
        {
            return new CharacteristicMatrix(this).CharacteristicPolynomial();
        }

        /// <summary>
        /// Returns a list of the real eigenvalues of the matrix using Newton's Method (Note: The same eigenvalues can be listed multiple times depending on their respective multiplicities.)
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public double[] GetRealEigenvalues(params ComplexNumber[] knownEigenvalues)
        {
            if (m >= 10)
            {
                throw new LinearAlgebraException("This program only supports eigenvalue decomposition for matrices of 9x9 or less.");
            }

            return new CharacteristicMatrix(this).CharacteristicPolynomial().FindRealRoots(knownEigenvalues);
        }

        /// <summary>
        /// Returns a list of all eigenvalues of the matrix using Newton's Method (Note: The same eigenvalues can be listed multiple times depending on their respective multiplicities.)
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public ComplexNumber[] GetComplexEigenvalues(params ComplexNumber[] knownEigenvalues)
        {
            if (m >= 10)
            {
                throw new LinearAlgebraException("This program only supports eigenvalue decomposition for matrices of 9x9 or less.");
            }

            return new CharacteristicMatrix(this).CharacteristicPolynomial().FindComplexRoots(knownEigenvalues);
        }

        /// <summary>
        /// Returns a list of the real eigenvectors associated with a particular eigenvalue
        /// </summary>
        /// <param name="eigenvalue">The eigenvalue associated with the eigenvectors to find</param>
        public Vector[] GetRealEigenvectors(double eigenvalue)
        {
            if (!IsSquare)
            {
                throw new LinearAlgebraException("Only square matrices have eigenvectors.");
            }

            Space eigenspace = (this - eigenvalue * Identity(m)).Kernel();
            if (eigenspace.Dimension == 0)
            {
                throw new LinearAlgebraException(eigenvalue + " is not an eigenvalue of the matrix.");
            }

            return eigenspace.Basis;
        }

        /// <summary>
        /// Returns a list of all eigenvectors associated with a particular eigenvalue
        /// </summary>
        /// <param name="eigenvalue">The eigenvalue associated with the eigenvectors to find</param>
        public ComplexVector[] GetComplexEigenvectors(ComplexNumber eigenvalue)
        {
            if (!IsSquare)
            {
                throw new LinearAlgebraException("Only square matrices have eigenvectors.");
            }

            ComplexMatrix complexMatrix = new ComplexMatrix(this, eigenvalue);
            return complexMatrix.GetKernelVectors();
        }

        /// <summary>
        /// Returns the eigenspace of a particular real eigenvalue
        /// </summary>
        /// <param name="eigenvalue">The eigenvalue associated with the eigenspace to find</param>
        public Space GetEigenspace(double eigenvalue)
        {
            if (!IsSquare)
            {
                throw new LinearAlgebraException("Only square matrices have eigenspaces.");
            }

            Space eigenspace = (this - eigenvalue * Identity(m)).Kernel();
            if (eigenspace.Dimension == 0)
            {
                throw new LinearAlgebraException(eigenvalue + " is not an eigenvalue of the matrix.");
            }

            return eigenspace;
        }

        /// <summary>
        /// Returns a summary of the real eigenvalues and eigenvectors of the matrix, stored as eigenpairs
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public Eigenpair[] GetRealEigenpairs(params ComplexNumber[] knownEigenvalues)
        {
            List<double> eigenvalues = new List<double>(GetRealEigenvalues(knownEigenvalues));
            List<Eigenpair> eigenpairs = new List<Eigenpair>();
            foreach (IGrouping<double, double> eigenvalue in eigenvalues.GroupBy(i => i))
            {
                Eigenpair eigenpair = new Eigenpair(eigenvalue.Key, GetEigenspace(eigenvalue.Key), eigenvalue.Count());
                eigenpairs.Add(eigenpair);
            }

            return eigenpairs.ToArray();
        }

        /// <summary>
        /// Returns a summary of all eigenvalues and eigenvectors of the matrix, stored as eigenpairs
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public ComplexEigenpair[] GetComplexEigenpairs(params ComplexNumber[] knownEigenvalues)
        {
            List<ComplexNumber> eigenvalues = new List<ComplexNumber>(GetComplexEigenvalues(knownEigenvalues));
            List<ComplexEigenpair> eigenpairs = new List<ComplexEigenpair>();
            foreach (IGrouping<ComplexNumber, ComplexNumber> eigenvalue in eigenvalues.GroupBy(i => i))
            {
                ComplexEigenpair eigenpair = new ComplexEigenpair(eigenvalue.Key, GetComplexEigenvectors(eigenvalue.Key), eigenvalue.Count());
                eigenpairs.Add(eigenpair);
            }

            return eigenpairs.ToArray();
        }

        /// <summary>
        /// Returns a diagonal matrix similar to this matrix
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public Matrix Diagonalize(params ComplexNumber[] knownEigenvalues)
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            int eigenvectors = 0;
            List<double> eigenvalues = new List<double>();
            Eigenpair[] eigenpairs = GetRealEigenpairs(knownEigenvalues);
            for (int i = 0; i < eigenpairs.Length; i++)
            {
                for (int j = 0; j < eigenpairs[i].eigenvectors.Length; j++)
                {
                    eigenvalues.Add(eigenpairs[i].eigenvalue);
                    eigenvectors++;
                }
            }

            if (eigenvectors != m)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            return new Matrix(eigenvalues.ToArray());
        }

        /// <summary>
        /// Returns a diagonal matrix similar to this matrix
        /// </summary>
        /// <param name="ChangeOfBasisMatrix">The change of basis matrix for the diagonalization</param>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public Matrix Diagonalize(out Matrix ChangeOfBasisMatrix, params ComplexNumber[] knownEigenvalues)
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            List<double> eigenvalues = new List<double>();
            List<Vector> eigenvectors = new List<Vector>();
            Eigenpair[] eigenpairs = GetRealEigenpairs(knownEigenvalues);
            for (int i = 0; i < eigenpairs.Length; i++)
            {
                for (int j = 0; j < eigenpairs[i].eigenvectors.Length; j++)
                {
                    eigenvalues.Add(eigenpairs[i].eigenvalue);
                    eigenvectors.Add(eigenpairs[i].eigenvectors[j]);
                }
            }

            if (eigenvectors.Count != m)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            ChangeOfBasisMatrix = Vector.ToMatrix(eigenvectors.ToArray());
            return new Matrix(eigenvalues.ToArray());
        }

        /// <summary>
        /// Returns a rotation matrix similar to this matrix
        /// </summary>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public Matrix DiagonalizeToRotationMatrix(params ComplexNumber[] knownEigenvalues)
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            List<Vector> diagonalColumns = new List<Vector>();
            ComplexEigenpair[] eigenpairs = GetComplexEigenpairs(knownEigenvalues);
            bool complexEigenvalues = false;
            for (int i = 0; i < eigenpairs.Length; i++)
            {
                for (int j = 0; j < eigenpairs[i].eigenvectors.Length; j++)
                {
                    int pos = diagonalColumns.Count;
                    if (eigenpairs[i].eigenvalue.IsReal)
                    {
                        diagonalColumns.Add((double)eigenpairs[i].eigenvalue * Vector.StandardUnit(pos, m));
                    }
                    else if (!complexEigenvalues)
                    {
                        diagonalColumns.Add((eigenpairs[i].eigenvalue.a * Vector.StandardUnit(pos, m)) + (eigenpairs[i].eigenvalue.b * Vector.StandardUnit(pos + 1, m)));
                        diagonalColumns.Add((-eigenpairs[i].eigenvalue.b * Vector.StandardUnit(pos, m)) + (eigenpairs[i].eigenvalue.a * Vector.StandardUnit(pos + 1, m)));
                        complexEigenvalues = true;
                    }
                }
            }

            if (diagonalColumns.Count != m)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable to the rotation matrix. Try calling the Diagonalize() method instead.");
            }

            return Vector.ToMatrix(diagonalColumns.ToArray());
        }

        /// <summary>
        /// Returns a rotation matrix similar to this matrix
        /// </summary>
        /// <param name="ChangeOfBasisMatrix">The change of basis matrix for the diagonalization</param>
        /// <param name="knownEigenvalues">Any eigenvalues (real or complex) that do not need to be found with Newton's Method (Note: For imaginary eigenvalues, do NOT list their conjugates here.)</param>
        public Matrix DiagonalizeToRotationMatrix(out Matrix ChangeOfBasisMatrix, params ComplexNumber[] knownEigenvalues)
        {
            if (m != n)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable.");
            }

            List<Vector> diagonalColumns = new List<Vector>();
            List<Vector> changeOfBasisColumns = new List<Vector>();
            ComplexEigenpair[] eigenpairs = GetComplexEigenpairs(knownEigenvalues);
            bool complexEigenvalues = false;
            for (int i = 0; i < eigenpairs.Length; i++)
            {
                for (int j = 0; j < eigenpairs[i].eigenvectors.Length; j++)
                {
                    int pos = diagonalColumns.Count;
                    if (eigenpairs[i].eigenvalue.IsReal)
                    {
                        diagonalColumns.Add((double)eigenpairs[i].eigenvalue * Vector.StandardUnit(pos, m));
                        changeOfBasisColumns.Add(eigenpairs[i].eigenvectors[j].realVector);
                    }
                    else if (!complexEigenvalues)
                    {
                        diagonalColumns.Add((eigenpairs[i].eigenvalue.a * Vector.StandardUnit(pos, m)) + (eigenpairs[i].eigenvalue.b * Vector.StandardUnit(pos + 1, m)));
                        diagonalColumns.Add((-eigenpairs[i].eigenvalue.b * Vector.StandardUnit(pos, m)) + (eigenpairs[i].eigenvalue.a * Vector.StandardUnit(pos + 1, m)));
                        changeOfBasisColumns.Add(eigenpairs[i].eigenvectors[j].complexVector);
                        changeOfBasisColumns.Add(eigenpairs[i].eigenvectors[j].realVector);
                        complexEigenvalues = true;
                    }
                }
            }

            if (diagonalColumns.Count != m)
            {
                throw new LinearAlgebraException("The matrix is not diagonalizable to the rotation matrix. Try calling the Diagonalize() method instead.");
            }

            ChangeOfBasisMatrix = Vector.ToMatrix(changeOfBasisColumns.ToArray());
            return Vector.ToMatrix(diagonalColumns.ToArray());
        }

        /// <summary>
        /// Returns a list of the columns of the matrix as vectors
        /// </summary>
        public Vector[] ToVectorArray()
        {
            Vector[] vectors = new Vector[n];
            for (int j = 0; j < n; j++)
            {
                double[] elements = new double[m];
                for (int i = 0; i < m; i++)
                {
                    elements[i] = rElements[i, j];
                }

                vectors[j] = new Vector(elements);
            }
            return vectors;
        }

        /// <summary>
        /// Internal method for finding the determinant using row-reduction
        /// </summary>
        internal Matrix RREF(out double determinant)
        {
            if (n != m)
            {
                determinant = 0;
            }
            else
            {
                determinant = 1;
            }

            double[,] rElements = (double[,])this.rElements.Clone();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    rElements[i, j] = Math.Round(rElements[i, j], 4);
                    if (rElements[i, j] == 0 && j >= i)
                    {
                        for (int k = (i + 1); k < m; k++)
                        {
                            if (rElements[k, j] != 0)
                            {
                                determinant *= -1;
                                for (int l = 0; l < n; l++)
                                {
                                    double placeholder = rElements[k, l];
                                    rElements[k, l] = rElements[i, l];
                                    rElements[i, l] = placeholder;
                                }

                                j = -1;
                                break;
                            }
                        }
                    }
                    else if (rElements[i, j] != 0)
                    {
                        bool zeroPossible = false;
                        for (int k = i - 1; k >= 0; k--)
                        {
                            if (rElements[k, j] == 1)
                            {
                                int l = 0;
                                while (rElements[k, l] == 0) l++;

                                if (l == j)
                                {
                                    double factor = rElements[i, j];
                                    for (int h = 0; h < n; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[k, h];
                                    }

                                    zeroPossible = true;
                                    break;
                                }
                            }
                        }
                        if (!zeroPossible)
                        {
                            double factor = rElements[i, j];
                            determinant *= factor;
                            for (int k = 0; k < n; k++)
                            {
                                rElements[i, k] /= factor;
                            }
                            break;
                        }
                    }
                }
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = (i + 1); j < n; j++)
                {
                    if (rElements[i, j] != 0)
                    {
                        double factor = rElements[i, j];
                        for (int l = Math.Min(j, m - 1); l > i; l--)
                        {
                            if (rElements[l, j] == 1)
                            {
                                int k = 0;
                                while (rElements[l, k] == 0) k++;

                                if (k == j)
                                {
                                    for (int h = 0; h < n; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[l, h];
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            if (rElements[m - 1, n - 1] == 0)
            {
                determinant = 0;
            }

            return new Matrix(rElements, true);
        }

        public override string ToString()
        {
            string finalStr = "";
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    finalStr += (rElements[i, j] + " ");
                }

                finalStr += "\n";
            }

            return finalStr;
        }
    }

    /// <summary>
    /// An internal class used to find the characteristic polynomial of a matrix
    /// </summary>
    internal class CharacteristicMatrix
    {
        private Polynomial[,] rElements;
        private int dim;

        /// <summary>
        /// Initializes a characteristic matrix with given elements
        /// </summary>
        public CharacteristicMatrix(Polynomial[,] rElements)
        {
            if (rElements.GetLength(0) != rElements.GetLength(1))
            {
                throw new LinearAlgebraException("A matrix must be square in order to have a characteristic polynomial.");
            }

            this.rElements = rElements;
            this.dim = rElements.GetLength(0);
        }

        /// <summary>
        /// Initializes a characteristic matrix from a real matrix
        /// </summary>
        public CharacteristicMatrix(Matrix a)
        {
            if (!a.IsSquare)
            {
                throw new LinearAlgebraException("A matrix must be square in order to have a characteristic polynomial.");
            }

            dim = a.Rows;
            rElements = new Polynomial[a.Rows, a.Cols];
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Cols; j++)
                {
                    if (i == j)
                    {
                        rElements[i, j] = new Polynomial(new double[] { a[i, j], -1 });
                    }
                    else
                    {
                        rElements[i, j] = new Polynomial(new double[] { a[i, j] });
                    }
                }
            }
        }

        /// <summary>
        /// Removes a row and column from the matrix (used for finding matrix minors)
        /// </summary>
        /// <param name="row">The row to remove</param>
        /// <param name="col">The column to remove</param>
        public void RemoveRowAndCol(int row, int col)
        {
            Polynomial[,] rElements = new Polynomial[dim - 1, dim - 1];
            for (int i = 0; i < dim; i++)
            {
                if (i == row)
                {
                    continue;
                }
                for (int j = 0; j < dim; j++)
                {
                    if (j == col)
                    {
                        continue;
                    }
                    if (i > row)
                    {
                        if (j > col)
                        {
                            rElements[i - 1, j - 1] = this.rElements[i, j];
                        }
                        else
                        {
                            rElements[i - 1, j] = this.rElements[i, j];
                        }
                    }
                    else
                    {
                        if (j > col)
                        {
                            rElements[i, j - 1] = this.rElements[i, j];
                        }
                        else
                        {
                            rElements[i, j] = this.rElements[i, j];
                        }
                    }
                }
            }
            this.rElements = rElements;
            dim--;
        }

        /// <summary>
        /// Determines the characteristic polynomial of the matrix by using cofactor expansion
        /// </summary>
        public Polynomial CharacteristicPolynomial()
        {
            if (dim == 2)
            {
                return (rElements[0, 0] * rElements[1, 1]) - (rElements[1, 0] * rElements[0, 1]);
            }

            Polynomial det = new Polynomial(dim);
            for (int i = 0; i < dim; i++)
            {
                CharacteristicMatrix subMatrix = new CharacteristicMatrix(rElements);
                subMatrix.RemoveRowAndCol(0, i);
                if (i % 2 == 0)
                {
                    det += rElements[0, i] * subMatrix.CharacteristicPolynomial();
                }
                else
                {
                    det -= rElements[0, i] * subMatrix.CharacteristicPolynomial();
                }
            }

            return det;
        }
    }

    /// <summary>
    /// An internal class used to find complex eigenvectors that correspond to a given eigenvalue
    /// </summary>
    internal class ComplexMatrix
    {
        private ComplexNumber[,] rElements;
        private int dim;

        /// <summary>
        /// Initializes a complex matrix with given elements
        /// </summary>
        public ComplexMatrix(ComplexNumber[,] rElements)
        {
            if (rElements.GetLength(0) != rElements.GetLength(1))
            {
                throw new LinearAlgebraException("A matrix must be square in order to be a complex matrix.");
            }

            this.rElements = rElements;
            this.dim = rElements.GetLength(0);
        }

        /// <summary>
        /// Initializes a complex matrix from a real matrix and a complex eigenvalue
        /// </summary>
        public ComplexMatrix(Matrix a, ComplexNumber eigenvalue)
        {
            if (!a.IsSquare)
            {
                throw new LinearAlgebraException("A matrix must be square in order to have a compex matrix.");
            }

            this.dim = a.Rows;
            this.rElements = new ComplexNumber[a.Rows, a.Cols];
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Cols; j++)
                {
                    if (i == j)
                    {
                        rElements[i, j] = (a[i, j] - eigenvalue);
                    }
                    else
                    {
                        rElements[i, j] = a[i, j];
                    }
                }
            }
        }

        /// <summary>
        /// Finds the reduced row echelon form of this matrix
        /// </summary>
        public ComplexMatrix RREF()
        {
            ComplexNumber[,] rElements = (ComplexNumber[,])this.rElements.Clone();
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    if (ComplexNumber.Round(rElements[i, j], (-dim + 9).Clamp(0, 6)) == 0 && j >= i)
                    {
                        rElements[i, j] = 0;
                        for (int k = (i + 1); k < dim; k++)
                        {
                            if (rElements[k, j] != 0)
                            {
                                for (int l = 0; l < dim; l++)
                                {
                                    ComplexNumber placeholder = rElements[k, l];
                                    rElements[k, l] = rElements[i, l];
                                    rElements[i, l] = placeholder;
                                }

                                j = -1;
                                break;
                            }
                        }
                    }
                    else if (rElements[i, j] != 0)
                    {
                        bool zeroPossible = false;
                        for (int k = i - 1; k >= 0; k--)
                        {
                            if (rElements[k, j] == 1)
                            {
                                int l = 0;
                                while (rElements[k, l] == 0) l++;

                                if (l == j)
                                {
                                    ComplexNumber factor = rElements[i, j];
                                    for (int h = 0; h < dim; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[k, h];
                                    }

                                    zeroPossible = true;
                                    break;
                                }
                            }
                        }
                        if (!zeroPossible)
                        {
                            ComplexNumber factor = rElements[i, j];
                            for (int k = 0; k < dim; k++)
                            {
                                rElements[i, k] /= factor;
                            }
                            break;
                        }
                    }
                }
            }
            for (int i = 0; i < dim; i++)
            {
                for (int j = (i + 1); j < dim; j++)
                {
                    if (rElements[i, j] != 0)
                    {
                        ComplexNumber factor = rElements[i, j];
                        for (int l = Math.Min(j, dim - 1); l > i; l--)
                        {
                            if (rElements[l, j] == 1)
                            {
                                int k = 0;
                                while (rElements[l, k] == 0) k++;

                                if (k == j)
                                {
                                    for (int h = 0; h < dim; h++)
                                    {
                                        rElements[i, h] -= factor * rElements[l, h];
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            ComplexMatrix result = new ComplexMatrix(rElements);
            result.Round(6);
            return result;
        }

        /// <summary>
        /// Determines the kernel vectors, or eigenvectors, of the matrix
        /// </summary>
        public ComplexVector[] GetKernelVectors()
        {
            ComplexMatrix rref = RREF();
            List<ComplexVector> kernelVectors = new List<ComplexVector>();
            bool[] columnsNormalized = new bool[rref.dim];
            for (int j = 0; j < rref.dim; j++)
            {
                bool allZero = true;
                for (int i = 0; i < rref.dim; i++)
                {
                    if (rref.rElements[i, j] != 0)
                    {
                        allZero = false;
                        break;
                    }
                }
                if (allZero)
                {
                    columnsNormalized[j] = true;
                    kernelVectors.Add(new ComplexVector(Vector.StandardUnit(j, rref.dim), Vector.Zero(rref.dim)));
                }
            }

            for (int i = 0; i < rref.dim; i++)
            {
                int leadingOneCol = -1;
                for (int j = 0; j < rref.dim; j++)
                {
                    ComplexVector kernelVec = ComplexVector.Zero(rref.dim);
                    if (rref.rElements[i, j] == 1 && leadingOneCol == -1)
                    {
                        leadingOneCol = j;
                    }
                    else if (ComplexNumber.Round(rref.rElements[i, j], 5) != 0 && leadingOneCol != -1 && !columnsNormalized[j])
                    {
                        ComplexNumber el = ComplexNumber.Round(rref.rElements[i, j], 6);
                        kernelVec.realVector[leadingOneCol] = el.a;
                        kernelVec.complexVector[leadingOneCol] = el.b;
                        kernelVec.realVector[j] = -1;
                        kernelVec.complexVector[j] = 0;
                        for (int k = (i + 1); k < rref.dim; k++)
                        {
                            if (rref.rElements[k, j] != 0)
                            {
                                int l = 0;
                                while (rref.rElements[k, l] != 1) l++;

                                ComplexNumber el2 = ComplexNumber.Round(rref.rElements[k, j], 5);
                                kernelVec.realVector[l] = el2.a;
                                kernelVec.complexVector[l] = el2.b;
                            }
                        }
                        columnsNormalized[j] = true;
                        kernelVectors.Add(kernelVec);
                    }
                }
            }

            return kernelVectors.ToArray();
        }

        public void Round(int digits)
        {
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    rElements[i, j] = ComplexNumber.Round(rElements[i, j], digits);
                }
            }
        }

        public override string ToString()
        {
            string finalStr = "";
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    finalStr += (rElements[i, j] + " ");
                }
                finalStr += "\n";
            }

            return finalStr;
        }
    }
}