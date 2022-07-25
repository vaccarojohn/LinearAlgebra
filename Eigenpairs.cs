namespace LinearAlgebra
{
    /// <summary>
    /// A struct providing a representation of the eigenvector(s), eigenspace, and algebraic/geometric multiplicities associated with a particular real eigenvalue
    /// </summary>
    public struct Eigenpair
    {
        /// <summary>
        /// The eigenvalue of the pair
        /// </summary>
        public readonly double eigenvalue;

        /// <summary>
        /// A list of linearly independent eigenvectors corresponding to the given eigenvalue
        /// </summary>
        public readonly Vector[] eigenvectors;

        /// <summary>
        /// The eigenspace associated with the eigenvalue
        /// </summary>
        public readonly Space eigenspace;

        /// <summary>
        /// The algebraic multiplicity of the eigenvalue, or the multiplicity of the eigenvalue in a matrix's characteristic polynomial
        /// </summary>
        public readonly int almu;

        /// <summary>
        /// The geometric multiplicity of the eigenvalue, or dimension of its eigenspace
        /// </summary>
        public readonly int gemu;

        /// <summary>
        /// Initializes an eigenpair with a particular eigenvalue, eigenspace, and algebraic multiplicity
        /// </summary>
        /// <param name="eigenvalue">The eigenvalue of the pair</param>
        /// <param name="eigenspace">The eigenspace of the eigenvalue</param>
        /// <param name="almu">The algebraic multiplicity of the eigenvalue</param>
        public Eigenpair(double eigenvalue, Space eigenspace, int almu)
        {
            if (almu < eigenspace.Dimension)
            {
                throw new LinearAlgebraException("The algebraic multiplicity of the eigenvalue must be at least its geometric multiplicity (the dimension of its eigenspace).");
            }

            this.eigenvalue = eigenvalue;
            this.eigenspace = eigenspace;
            this.eigenvectors = eigenspace.Basis;
            this.almu = almu;
            this.gemu = eigenspace.Dimension;
        }

        public override string ToString()
        {
            string myStr = "Eigenpair: λ [lambda] = " + eigenvalue + ", and ";
            for (int i = 0; i < eigenvectors.Length; i++)
            {
                string vectorStr = "v" + i + " = " + eigenvectors[i];
                if (i < (eigenvectors.Length - 1))
                {
                    vectorStr += ", ";
                }

                myStr += vectorStr;
            }
            return myStr;
        }
    }

    /// <summary>
    /// A struct providing a representation of the eigenvector(s) and algebraic/geometric multiplicities associated with a particular complex eigenvalue
    /// </summary>
    public struct ComplexEigenpair
    {
        /// <summary>
        /// The eigenvalue of the pair
        /// </summary>
        public readonly ComplexNumber eigenvalue;

        /// <summary>
        /// A list of linearly independent eigenvectors corresponding to the given eigenvalue
        /// </summary>
        public readonly ComplexVector[] eigenvectors;

        /// <summary>
        /// The algebraic multiplicity of the eigenvalue, or the multiplicity of the eigenvalue in a matrix's characteristic polynomial
        /// </summary>
        public readonly int almu;

        /// <summary>
        /// The geometric multiplicity of the eigenvalue, or the number of linearly independent eigenvectors associated with it
        /// </summary>
        public readonly int gemu;

        /// <summary>
        /// Initializes a complex eigenpair with a particular eigenvalue, list of eigenvectors, and algebraic multiplicity
        /// </summary>
        /// <param name="eigenvalue">The eigenvalue of the pair</param>
        /// <param name="eigenvectors">The eigenvectors associated with the eigenvalue</param>
        /// <param name="almu">The algebraic multiplicity of the eigenvalue</param>
        public ComplexEigenpair(ComplexNumber eigenvalue, ComplexVector[] eigenvectors, int almu)
        {
            this.eigenvalue = eigenvalue;
            this.eigenvectors = eigenvectors;
            this.almu = almu;
            this.gemu = eigenvectors.Length;
        }

        public override string ToString()
        {
            string myStr = "Eigenpair: λ [lambda] = " + eigenvalue + ", and ";
            for (int i = 0; i < eigenvectors.Length; i++)
            {
                string vectorStr = "v" + i + " = " + eigenvectors[i];

                if (i < (eigenvectors.Length - 1))
                {
                    vectorStr += ", ";
                }

                myStr += vectorStr;
            }
            return myStr;
        }
    }
}