using System;

namespace LinearAlgebra
{
    /// <summary>
    /// An ordered collection of doubles representing a vector in n-dimensional space
    /// </summary>
    public class Vector
    {

        /// <summary>
        /// The method of getting or changing a value at a specified position in the vector
        /// </summary>
        /// <param name="n">The position of the value to get or set in the range (0 inclusive)</param>
        public double this[int n] { get { return elements[n]; } set { elements[n] = value; } }

        /// <summary>
        /// Returns the number of values in the vector
        /// </summary>
        public int Dimension { get { return elements.Length; } }

        /// <summary>
        /// Returns the magnitude of the vector
        /// </summary>
        public double Magnitude { get { return Math.Sqrt(SquareMagnitude); } }

        /// <summary>
        /// Returns the square magnitude of the vector (alternatively v • v)
        /// </summary>
        public double SquareMagnitude
        {
            get
            {
                double magSq = 0;
                for (int i = 0; i < elements.Length; i++)
                {
                    magSq += elements[i] * elements[i];
                }

                return Math.Round(magSq, 6);
            }
        }

        private double[] elements;

        /// <summary>
        /// Returns the zero vector in a specified number of dimensions
        /// </summary>
        /// <param name="dim">The dimension of the vector</param>
        public static Vector Zero(int dim)
        {
            double[] list = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                list[i] = 0;
            }

            return new Vector(list);
        }

        /// <summary>
        /// Returns a standard unit vector in a specified number of dimensions
        /// </summary>
        /// <param name="pos">The position of the 1 in the vector</param>
        /// <param name="dim">The dimension of the vector</param>
        public static Vector StandardUnit(int pos, int dim)
        {
            if (dim < 1)
            {
                throw new LinearAlgebraException("The dimension of the vector must be at least 1");
            }

            if (pos >= dim || pos < 0)
            {
                throw new LinearAlgebraException("The position of the 1 should be in the range [0, " + (dim - 1) + "].");
            }

            double[] list = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                if (i == pos)
                {
                    list[i] = 1;
                }
                else
                {
                    list[i] = 0;
                }
            }

            return new Vector(list);
        }

        /// <summary>
        /// Converts a list of vectors into a matrix (Note: All vectors in the list must have the same dimension.)
        /// </summary>
        public static Matrix ToMatrix(params Vector[] vectors)
        {
            if (vectors.Length != 0)
            {
                int vectorDimension = vectors[0].Dimension;
                double[,] elements = new double[vectors.Length, vectorDimension];
                for (int i = 0; i < vectors.Length; i++)
                {
                    if (vectors[i].Dimension != vectorDimension)
                    {
                        throw new LinearAlgebraException("All vectors in the array must have the same dimension in order to be converted into a matrix.");
                    }

                    for (int j = 0; j < vectors[i].Dimension; j++)
                    {
                        elements[i, j] = vectors[i][j];
                    }
                }

                return new Matrix(elements, false);
            }
            else
            {
                throw new LinearAlgebraException("There must be at least one vector in the vector list to create a matrix.");
            }
        }

        /// <summary>
        /// Returns a vector representing the cross product in three dimensions between two vectors
        /// </summary>
        public static Vector CrossProduct(Vector a, Vector b)
        {
            if (a.Dimension == 3 && b.Dimension == 3)
            {
                return new Vector(new double[3] { Math.Round((a[1] * b[2]) - (a[2] * b[1]), 6), Math.Round((a[2] * b[0]) - (a[0] * b[2]), 6), Math.Round((a[0] * b[1]) - (a[1] * b[0]), 6) });
            }
            else
            {
                throw new LinearAlgebraException("Both vectors must have a dimension of 3 to use the cross product.");
            }
        }

        /// <summary>
        /// Flips the sign of all values in the vector
        /// </summary>
        public static Vector operator -(Vector a)
        {
            return (-1 * a);
        }

        /// <summary>
        /// Returns a vector representing the sum of two vectors
        /// </summary>
        public static Vector operator +(Vector a, Vector b)
        {
            if (a.Dimension != b.Dimension)
            {
                throw new LinearAlgebraException("Vectors must be of the same dimension.");
            }

            double[] elements = new double[a.Dimension];
            for (int i = 0; i < a.Dimension; i++)
            {
                elements[i] = a[i] + b[i];
            }

            return new Vector(elements);
        }

        /// <summary>
        /// Returns a vector representing the difference between two vectors
        /// </summary>
        public static Vector operator -(Vector a, Vector b)
        {
            if (a.Dimension != b.Dimension)
            {
                throw new LinearAlgebraException("Vectors must be of the same dimension.");
            }

            double[] elements = new double[a.Dimension];
            for (int i = 0; i < a.Dimension; i++)
            {
                elements[i] = a[i] - b[i];
            }

            return new Vector(elements);
        }

        /// <summary>
        /// Returns a vector representing the multiplication of a vector with a scalar
        /// </summary>
        public static Vector operator *(Vector a, double b)
        {
            double[] elements = new double[a.Dimension];
            for (int i = 0; i < a.Dimension; i++)
            {
                elements[i] = b * a[i];
            }

            return new Vector(elements);
        }

        /// <summary>
        /// Returns a vector representing the multiplication of a vector with a scalar
        /// </summary>
        public static Vector operator *(double a, Vector b)
        {
            double[] elements = new double[b.Dimension];
            for (int i = 0; i < b.Dimension; i++)
            {
                elements[i] = a * b[i];
            }

            return new Vector(elements);
        }

        /// <summary>
        /// Returns the dot product of two vectors
        /// </summary>
        public static double operator *(Vector a, Vector b)
        {
            if (a.Dimension != b.Dimension)
            {
                throw new LinearAlgebraException("Vectors must be of the same dimension.");
            }

            double result = 0;
            for (int i = 0; i < a.Dimension; i++)
            {
                result += a[i] * b[i];
            }

            return result;
        }

        /// <summary>
        /// Returns a vector representing the multiplication of a vector with the inverse of a scalar
        /// </summary>
        public static Vector operator /(Vector a, double b)
        {
            if (b == 0)
            {
                throw new DivideByZeroException();
            }

            double[] elements = new double[a.Dimension];
            for (int i = 0; i < a.Dimension; i++)
            {
                elements[i] = (1 / b) * a[i];
            }

            return new Vector(elements);
        }

        /// <summary>
        /// Determines whether two vectors are equal (Note: Vectors must have the same dimension and values to be considered equal.)
        /// </summary>
        public static bool operator ==(Vector a, Vector b)
        {
            return a.Equals(b);
        }

        public static bool operator !=(Vector a, Vector b)
        {
            return !a.Equals(b);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Vector b = (Vector)obj;

            if (elements.Length != b.Dimension)
            {
                return false;
            }

            for (int i = 0; i < elements.Length; i++)
            {
                if (elements[i] != b[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Initializes a vector with the specified values
        /// </summary>
        /// <param name="elements">Ordered list of values for the vector</param>
        public Vector(params double[] elements)
        {
            this.elements = elements;
        }

        /// <summary>
        /// Returns the angle in radians between this vector and another vector
        /// </summary>
        public double AngleWith(Vector other)
        {
            return Math.Round(Math.Acos((this * other) / (Magnitude * other.Magnitude)), 6);
        }

        /// <summary>
        /// Projects another vector onto this vector
        /// </summary>
        /// <param name="other">The vector to project onto this vector</param>
        public Vector ProjectFrom(Vector other)
        {
            return ((other * this) / SquareMagnitude) * this;
        }

        /// <summary>
        /// Reflects another vector across this vector
        /// </summary>
        /// <param name="other">The vector to reflect across this vector</param>
        public Vector ReflectFrom(Vector other)
        {
            return (-1 * other) + 2 * ProjectFrom(other);
        }

        /// <summary>
        /// Changes the dimension of a vector, truncating values or adding 0s when necessary
        /// </summary>
        /// <param name="newDim">The new dimension of the vector</param>
        public Vector Resize(int newDim)
        {
            double[] newElements = new double[newDim];
            for (int i = 0; i < newDim; i++)
            {
                if (i > elements.Length)
                {
                    newElements[i] = 0;
                }
                else
                {
                    newElements[i] = elements[i];
                }
            }

            return new Vector(newElements);
        }

        /// <summary>
        /// Returns the elements of the vector as a double[] array
        /// </summary>
        public double[] ToArray()
        {
            return elements;
        }

        /// <summary>
        /// Converts the vector into matrix with one column
        /// </summary>
        public Matrix ToMatrix()
        {
            double[,] elements = new double[1, this.elements.Length];
            for (int i = 0; i < this.elements.Length; i++)
            {
                elements[0, i] = this.elements[i];
            }

            return new Matrix(elements, false);
        }

        /// <summary>
        /// Represents the vector as a string (ex: "<![CDATA[<1, 1, 1>]]>")
        /// </summary>
        public override string ToString()
        {
            string finalStr = "<";
            for (int i = 0; i < elements.Length; i++)
            {
                if (i == (elements.Length - 1))
                {
                    finalStr += (elements[i] + ">");
                }
                else
                {
                    finalStr += (elements[i] + ", ");
                }
            }

            return finalStr;
        }
    }

    /// <summary>
    /// A vector containing a real and imaginary part
    /// </summary>
    public struct ComplexVector
    {
        /// <summary>
        /// The real part of the complex vector 
        /// </summary>
        public readonly Vector realVector;

        /// <summary>
        /// The imaginary part of the complex vector
        /// </summary>
        public readonly Vector complexVector;

        /// <summary>
        /// The zero vector in a specified number of dimensions, represented as a complex vector
        /// </summary>
        /// <param name="dim">The dimension of the vector</param>
        public static ComplexVector Zero(int dim)
        {
            return new ComplexVector(Vector.Zero(dim), Vector.Zero(dim));
        }

        /// <summary>
        /// Flips the sign of all values of the complex vector
        /// </summary>
        public static ComplexVector operator -(ComplexVector a)
        {
            return new ComplexVector(-a.realVector, -a.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the sum of two complex vectors
        /// </summary>
        public static ComplexVector operator +(ComplexVector a, ComplexVector b)
        {
            return new ComplexVector(a.realVector + b.realVector, a.complexVector + b.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the sum of a complex vector and a real vector
        /// </summary>
        public static ComplexVector operator +(ComplexVector a, Vector b)
        {
            return new ComplexVector(a.realVector + b, a.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the sum of a complex vector and a real vector
        /// </summary>
        public static ComplexVector operator +(Vector a, ComplexVector b)
        {
            return new ComplexVector(a + b.realVector, b.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the difference between two complex vectors
        /// </summary>
        public static ComplexVector operator -(ComplexVector a, ComplexVector b)
        {
            return new ComplexVector(a.realVector - b.realVector, a.complexVector - b.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the difference between a complex vector and a real vector
        /// </summary>
        public static ComplexVector operator -(ComplexVector a, Vector b)
        {
            return new ComplexVector(a.realVector - b, a.complexVector);
        }

        /// <summary>
        /// Returns a complex vector representing the difference between a real vector and a complex vector
        /// </summary>
        public static ComplexVector operator -(Vector a, ComplexVector b)
        {
            return new ComplexVector(a - b.realVector, -b.complexVector);
        }

        /// <summary>
        /// Initializes a complex vector given a real and imaginary part
        /// </summary>
        /// <param name="realVector">The real part of the complex vector</param>
        /// <param name="complexVector">The imaginary part of the complex vector</param>
        public ComplexVector(Vector realVector, Vector complexVector)
        {
            this.realVector = realVector;
            this.complexVector = complexVector;
        }

        public override string ToString()
        {
            return realVector.ToString() + " + i" + complexVector.ToString();
        }
    }
}