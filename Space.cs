using System;

namespace LinearAlgebra
{
    /// <summary>
    /// A linear space representing the spanning set of a particular basis
    /// </summary>
    public class Space
    {
        /// <summary>
        /// Returns a list of vectors in the stored basis of this linear space
        /// </summary>
        public Vector[] Basis { get { return basis; } }

        /// <summary>
        /// Returns the number of vectors in a basis of this linear space
        /// </summary>
        public int Dimension
        {
            get
            {
                try
                {
                    return basis.Length;
                }
                catch (NullReferenceException)
                {
                    return 0;
                }
            }
        }

        /// <summary>
        /// Returns the dimension of each basis vector
        /// </summary>
        public int ParentDimension { get { return parentDimension; } }

        /// <summary>
        /// Returns the change of basis matrix converting from basis coordinates into standard coordinates
        /// </summary>
        public Matrix ChangeOfBasisMatrix { get { return Vector.ToMatrix(basis); } }

        /// <summary>
        /// A boolean that is true if the stored basis is orthonormal
        /// </summary>
        public bool HasOrthonormalBasis
        {
            get
            {
                if (Dimension == 0)
                {
                    return false;
                }

                for (int i = 0; i < Dimension; i++)
                {
                    if (Math.Round(basis[i].SquareMagnitude, 6) != 1)
                    {
                        return false;
                    }

                    for (int j = 0; j < i; j++)
                    {
                        if (Math.Round(basis[i] * basis[j], 6) != 0)
                        {
                            return false;
                        }
                    }
                }

                return true;
            }
        }

        private int parentDimension;
        private Vector[] basis;

        /// <summary>
        /// Initializes a space that represents the origin in a specified number of dimensions
        /// </summary>
        /// <param name="n">The dimension of the space</param>
        public static Space Origin(int n)
        {
            return new Space(n);
        }

        /// <summary>
        /// Initializes a space that includes all vectors of a particular dimension
        /// </summary>
        /// <param name="n">The dimension of the vectors</param>
        public static Space R(int n)
        {
            Vector[] spanningSet = new Vector[n];

            for (int i = 0; i < n; i++)
            {
                spanningSet[i] = Vector.StandardUnit(i, n);
            }

            return new Space(spanningSet);
        }

        /// <summary>
        /// Determinates whether two spaces are equivalent (i.e. a vector is in space A if and only if it is also in space B)
        /// </summary>
        public static bool operator ==(Space a, Space b)
        {
            return a.Equals(b);
        }

        public static bool operator !=(Space a, Space b)
        {
            return !a.Equals(b);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Space b = (Space)obj;
            if ((Dimension != b.Dimension) || (ParentDimension != b.ParentDimension))
            {
                return false;
            }

            if (Dimension == 0)
            {
                return true;
            }

            Vector[] myVectors = new Vector[Dimension + b.Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                myVectors[i] = basis[i];
            }

            for (int i = 0; i < b.Dimension; i++)
            {
                myVectors[Dimension + i] = b.Basis[i];
            }

            return (Vector.ToMatrix(myVectors).Kernel().Dimension == b.Dimension);
        }

        /// <summary>
        /// Initializes an empty space (containing only the origin) that is a subset of all vectors of a particular dimension
        /// </summary>
        /// <param name="parentDimension">The dimension of the vectors in the space</param>
        public Space(int parentDimension)
        {
            if (parentDimension < 1)
            {
                throw new LinearAlgebraException("The parent dimension must be 1 or greater.");
            }

            this.parentDimension = parentDimension;
        }

        /// <summary>
        /// Initializes a space derived from a spanning set of vectors
        /// </summary>
        /// <param name="spanningSet">The list of vectors spanning the space (Note: it is not required for all vectors in the spanning set to be linearly independent)</param>
        public Space(params Vector[] spanningSet)
        {
            if (spanningSet.Length == 0)
            {
                throw new LinearAlgebraException("There must be at least one vector in the spanning set.");
            }
            else if (spanningSet.Length == 1)
            {
                parentDimension = spanningSet[0].Dimension;
                basis = spanningSet;
            }
            else
            {
                Matrix rref = Vector.ToMatrix(spanningSet).RREF();
                Vector[] basis = new Vector[spanningSet.Length];
                int basisVectorsIndex = 0;
                for (int i = 0; i < rref.Rows; i++)
                {
                    int k = 0;
                    while (k < rref.Cols && rref[i, k] == 0) k++;
                    if (k < rref.Cols)
                    {
                        basis[basisVectorsIndex] = spanningSet[k];
                        basisVectorsIndex++;
                    }
                }

                parentDimension = basis[0].Dimension;
                if (basisVectorsIndex != 0)
                {
                    this.basis = new Vector[basisVectorsIndex];
                    for (int i = 0; i < basisVectorsIndex; i++)
                    {
                        this.basis[i] = basis[i];
                    }
                }
            }
        }

        /// <summary>
        /// Returns true if the space contains a specified vector
        /// </summary>
        /// <param name="vector">The specified vector</param>
        public bool Contains(Vector vector)
        {
            if (vector.Dimension != parentDimension)
            {
                return false;
            }
            else if (vector.SquareMagnitude == 0)
            {
                return true;
            }
            else if (Dimension == 0)
            {
                return false;
            }

            Vector[] myVectors = new Vector[Dimension + 1];
            for (int i = 0; i < Dimension; i++)
            {
                myVectors[i] = basis[i];
            }

            myVectors[Dimension] = vector;
            Matrix rref = Vector.ToMatrix(myVectors).RREF();
            for (int i = 0; i < rref.Rows; i++)
            {
                Vector uv = Vector.StandardUnit(i, rref.Cols);
                if (rref * uv != uv)
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Converts a vector from standard coordinates into basis coordinates with respect to the stored basis
        /// </summary>
        /// <param name="vector">The vector to convert into basis coordinates</param>
        public Vector ToBasisCoordinates(Vector vector)
        {
            if (Dimension != ParentDimension)
            {
                throw new LinearAlgebraException("The basis does not span the entirety of R^" + ParentDimension + ", so vectors cannot be put in basis coordinates with respect to this basis.");
            }

            if (vector.Dimension != Dimension)
            {
                throw new LinearAlgebraException("The dimension of the vector must equal the dimension of the space.");
            }

            return ChangeOfBasisMatrix.Inverse() * vector;
        }

        /// <summary>
        /// Converts a vector from basis coordinates into standard coordinates with respect to the stored basis
        /// </summary>
        /// <param name="vector">The vector to convert into standard coordinates</param>
        public Vector ToStandardCoordinates(Vector vector)
        {
            if (vector.Dimension != Dimension)
            {
                throw new LinearAlgebraException("The dimension of the vector must equal the dimension of the space.");
            }

            return ChangeOfBasisMatrix * vector;
        }

        /// <summary>
        /// Returns the space perpendicular to this space
        /// </summary>
        /// <returns></returns>
        public Space PerpendicularSpace()
        {
            if (Dimension == 0)
            {
                Vector[] spanningSet = new Vector[parentDimension];
                for (int i = 0; i < ParentDimension; i++)
                {
                    spanningSet[i] = Vector.StandardUnit(i, parentDimension);
                }

                return new Space(spanningSet);
            }

            return Vector.ToMatrix(basis).Transpose.Kernel();
        }

        /// <summary>
        /// Projects a vector onto this space
        /// </summary>
        /// <param name="other">The vector to project onto this space</param>
        /// <param name="orthonormalizeThisSpace">Whether this space's basis should be permanently orthonormalized</param>
        public Vector ProjectFrom(Vector other, bool orthonormalizeThisSpace = false)
        {
            Space mySpace;
            if (orthonormalizeThisSpace)
            {
                mySpace = this;
            }
            else
            {
                mySpace = (Space)MemberwiseClone();
            }

            mySpace.Orthonormalize();
            Vector projection = Vector.Zero(mySpace.ParentDimension);
            for (int i = 0; i < mySpace.Dimension; i++)
            {
                projection += (other * mySpace.Basis[i]) * mySpace.Basis[i];
            }

            return projection;
        }

        /// <summary>
        /// Reflects a vector across this space
        /// </summary>
        /// <param name="other">The vector to reflect across this space</param>
        /// <param name="orthonormalizeThisSpace">Whether this space's basis should be permanently orthonormalized</param>
        public Vector ReflectFrom(Vector other, bool orthonormalizeThisSpace = false)
        {
            return (-1 * other) + (2 * ProjectFrom(other, orthonormalizeThisSpace));
        }

        /// <summary>
        /// Orthonormalizes this space's basis using the Gram-Schmidt process
        /// </summary>
        public void Orthonormalize()
        {
            for (int i = 0; i < Dimension; i++)
            {
                Vector parallelVector = Vector.Zero(parentDimension);
                for (int j = 0; j < i; j++)
                {
                    parallelVector += ((basis[i] * basis[j]) * basis[j]);
                }

                basis[i] -= parallelVector;
                basis[i] /= basis[i].Magnitude;
            }
        }

        public override string ToString()
        {
            if (Dimension != 0)
            {
                string finalStr = "Linear space with basis: ";
                for (int i = 0; i < Dimension; i++)
                {
                    finalStr += basis[i].ToString();
                    if (i != (Dimension - 1))
                    {
                        finalStr += ", ";
                    }
                }

                finalStr += ".";
                return finalStr;
            }
            else
            {

                return "Zero space with parent dimension: " + parentDimension;
            }
        }
    }
}