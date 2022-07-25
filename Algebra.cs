using System;
using System.Collections.Generic;

namespace LinearAlgebra
{
    /// <summary>
    /// A list of coefficients representing a polynomial function
    /// </summary>
    public struct Polynomial
    {
        /// <summary>
        /// The method of getting or changing a coefficient in the polynomial
        /// </summary>
        /// <param name="n">The position of the coefficient (0 inclusive)</param>
        public double this[int n] { get { return values[n]; } set { values[n] = value; } }

        /// <summary>
        /// The degree of the polynomial
        /// </summary>
        public int Degree { get { return degree; } }

        private double[] values;
        private int degree;

        /// <summary>
        /// Flips the sign on all polynomial coefficients
        /// </summary>
        public static Polynomial operator -(Polynomial a)
        {
            return (-1 * a);
        }

        /// <summary>
        /// Returns a polynomial representing the sum of two polynomials
        /// </summary>
        public static Polynomial operator +(Polynomial a, Polynomial b)
        {
            double[] values = new double[Math.Max(a.degree, b.degree) + 1];
            for (int i = 0; i <= Math.Max(a.degree, b.degree); i++)
            {
                if (i > a.degree)
                {
                    values[i] = b.values[i];
                }
                else if (i > b.degree)
                {
                    values[i] = a.values[i];
                }
                else
                {
                    values[i] = a.values[i] + b.values[i];
                }
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the difference between two polynomials
        /// </summary>
        public static Polynomial operator -(Polynomial a, Polynomial b)
        {
            double[] values = new double[Math.Max(a.degree, b.degree) + 1];
            for (int i = 0; i <= Math.Max(a.degree, b.degree); i++)
            {
                if (i > a.degree)
                {
                    values[i] = -1 * b.values[i];
                }
                else if (i > b.degree)
                {
                    values[i] = a.values[i];
                }
                else
                {
                    values[i] = a.values[i] - b.values[i];
                }
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the multiplication of a polynomial by a constant
        /// </summary>
        public static Polynomial operator *(Polynomial a, double b)
        {
            double[] values = new double[a.degree + 1];
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = b * a[i];
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the multiplication of a polynomial by a constant
        /// </summary>
        public static Polynomial operator *(double a, Polynomial b)
        {
            double[] values = new double[b.degree + 1];
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = a * b[i];
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the multiplication of two polynomials
        /// </summary>
        public static Polynomial operator *(Polynomial a, Polynomial b)
        {
            double[] values = new double[a.degree + b.degree + 1];
            for (int i = 0; i <= (a.degree + b.degree); i++)
            {
                values[i] = 0;
            }

            for (int i = 0; i <= a.degree; i++)
            {
                for (int j = 0; j <= b.degree; j++)
                {
                    values[i + j] += a.values[i] * b.values[j];
                }
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the multiplication of a polynomial with the inverse of a constant
        /// </summary>
        public static Polynomial operator /(Polynomial a, double b)
        {
            double[] values = a.values;
            for (int i = 0; i < values.Length; i++)
            {
                values[i] /= b;
            }

            return new Polynomial(values);
        }

        /// <summary>
        /// Returns a polynomial representing the quotient of two polynomials (Note: The operator will throw an exceptionn if the first polynomial is not divisible by the second one.)
        /// </summary>
        public static Polynomial operator /(Polynomial a, Polynomial b)
        {
            if (a.degree < b.degree)
            {
                throw new LinearAlgebraException("The polynomial " + a + " is not divisible by " + b + ".");
            }

            Polynomial dividend = (Polynomial)a.MemberwiseClone();
            Polynomial result = new Polynomial(a.degree - b.degree);
            for (int i = (a.degree - b.degree); i >= 0; i--)
            {
                double diff = dividend[i + b.degree] / b[b.degree];
                result[i] = diff;
                Polynomial diffPoly = new Polynomial(i);
                diffPoly[i] = diff;
                dividend -= diffPoly * b;
            }

            for (int i = 0; i <= dividend.degree; i++)
            {
                if (dividend[i] != 0)
                {
                    throw new LinearAlgebraException("The polynomial " + a + " is not divisible by " + b + ".");
                }
            }

            return result;
        }

        /// <summary>
        /// Determines whether two polynomials are equal (Note: Polynomials must have the same degree and coefficients to be considered equal.)
        /// </summary>
        public static bool operator ==(Polynomial a, Polynomial b)
        {
            return a.Equals(b);
        }

        public static bool operator !=(Polynomial a, Polynomial b)
        {
            return !a.Equals(b);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Polynomial b = (Polynomial)obj;
            if (degree != b.degree)
            {
                return false;
            }

            for (int i = 0; i < values.Length; i++)
            {
                if (values[i] != b.values[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Initializes a polynomial to a specified degree while filling coefficients with 0s
        /// </summary>
        /// <param name="degree">The degree of the new polynomial</param>
        public Polynomial(int degree)
        {
            if (degree < 0)
            {
                throw new LinearAlgebraException("The degree of the polynomial must be at least 0.");
            }

            this.degree = degree;
            this.values = new double[degree + 1];
            for (int i = 0; i <= degree; i++)
            {
                values[i] = 0;
            }
        }

        /// <summary>
        /// Initializes a polynomial with specified coefficients
        /// </summary>
        /// <param name="values">Coefficients listed in increasing degree order (Ex: 3x + 4x^2 - 7x^4 would be initialized as new Polynomial(0, 3, 4, 0, -7))</param>
        public Polynomial(params double[] values)
        {
            if (values.Length == 0)
            {
                this.degree = 0;
                this.values = new double[1] { 0 };
            }
            else
            {
                this.degree = values.Length - 1;
                this.values = values;
            }
        }

        /// <summary>
        /// Evaluates a polynomial for a given value of x
        /// </summary>
        /// <param name="x">The x-value of the polynomial</param>
        public double Evaluate(double x)
        {
            double val = 0;
            for (int i = 0; i <= degree; i++)
            {
                val += (values[i] * Math.Pow(x, i));
            }

            return val;
        }

        /// <summary>
        /// Evaluates a polynomial for a given value of x
        /// </summary>
        /// <param name="x">The x-value of the polynomial</param>
        public ComplexNumber Evaluate(ComplexNumber x)
        {
            ComplexNumber val = 0;
            for (int i = 0; i <= degree; i++)
            {
                val += (values[i] * ComplexNumber.Pow(x, i));
            }

            return val;
        }

        /// <summary>
        /// Attempts to divide a polynomial by a root
        /// </summary>
        /// <param name="r">The root of the polynomial</param>
        /// <param name="ignoreRemainder">When this value is false, an exception will be thrown if the polynomial is not divisible by the root.</param>
        public void DivideByRoot(double r, bool ignoreRemainder = false)
        {
            if (degree == 0)
            {
                throw new LinearAlgebraException("Cannot divide by a root when the degree of the polynomial is 0.");
            }

            double[] newElements = new double[degree];
            newElements[degree - 1] = values[degree];
            for (int i = (degree - 1); i > 0; i--)
            {
                newElements[i - 1] = newElements[i] * r + values[i];
            }

            if (Math.Round(newElements[0] * r + values[0], 7) != 0 && !ignoreRemainder)
            {
                throw new LinearAlgebraException(r + " is not a root of " + ToString() + ".");
            }

            this.values = newElements;
            this.degree--;
        }

        /// <summary>
        /// Attempts to divide a polynomial by a root
        /// </summary>
        /// <param name="r">The root of the polynomial</param>
        /// <param name="ignoreRemainder">When this value is false, an exception will be thrown if the polynomial is not divisible by the root.</param>
        public void DivideByRoot(ComplexNumber r, bool ignoreRemainder = false)
        {
            if ((degree == 0 || degree == 1) && !r.IsReal)
            {
                throw new LinearAlgebraException("Cannot divide by a complex root when the degree of the polynomial is 0 or 1.");
            }

            if (r.IsReal)
            {
                DivideByRoot(r.a, ignoreRemainder);
            }
            else
            {
                ComplexNumber[] newElements1 = new ComplexNumber[degree];
                newElements1[degree - 1] = new ComplexNumber(values[degree], 0, false);
                for (int i = (degree - 1); i > 0; i--)
                {
                    newElements1[i - 1] = newElements1[i] * r + values[i];
                }

                double[] newElements2 = new double[degree - 1];
                newElements2[degree - 2] = newElements1[degree - 1].a;
                for (int i = (degree - 2); i > 0; i--)
                {
                    newElements2[i - 1] = (newElements2[i] * r.Conjugate + newElements1[i]).a;
                }

                if ((double)ComplexNumber.Round(newElements2[0] * r.Conjugate + newElements1[0], 7) != 0 && !ignoreRemainder)
                {
                    throw new LinearAlgebraException(r + " is not a root of " + ToString() + ".");
                }

                this.values = newElements2;
                this.degree -= 2;
            }
        }

        /// <summary>
        /// Returns the derivative of the polynomial
        /// </summary>
        public Polynomial Derivative()
        {
            if (degree == 0)
            {
                return new Polynomial();
            }

            double[] newValues = new double[degree];
            for (int i = 1; i < values.Length; i++)
            {
                newValues[i - 1] = values[i] * i;
            }

            return new Polynomial(newValues);
        }

        /// <summary>
        /// Returns the anti-derivative of the polynomial (Note: The final constant term is set to 0 by default.)
        /// </summary>
        public Polynomial Antiderivative()
        {
            double[] newValues = new double[degree + 2];
            newValues[0] = 0;
            for (int i = 0; i < values.Length; i++)
            {
                newValues[i + 1] = values[i] / (i + 1);
            }

            return new Polynomial(newValues);
        }

        /// <summary>
        /// Returns a list of real roots of the polynomial using Newton's Method (Note: The same roots can be listed multiple times depending on their respective multiplicities.)
        /// </summary>
        /// <param name="knownRoots">Any roots (real or complex) that do not need to be found with Newton's Method (Note: For imaginary roots, do NOT list their conjugates here.)</param>
        public double[] FindRealRoots(params ComplexNumber[] knownRoots)
        {
            List<double> roots = new List<double>();
            Polynomial myPolynomial = (Polynomial)MemberwiseClone();
            foreach (ComplexNumber root in knownRoots)
            {
                if (root.IsReal)
                {
                    roots.Add((double)root);
                }
                myPolynomial.DivideByRoot(root);
            }

            while (myPolynomial.values[0] == 0)
            {
                roots.Add(0);
                myPolynomial.DivideByRoot(0);
            }

            while (myPolynomial.degree > 1)
            {
                bool stuck = true;
                ComplexNumber approximationComplex = ComplexNumber.i;
                for (int i = 0; i < 100; i++)
                {
                    ComplexNumber resComplex = myPolynomial.Evaluate(approximationComplex);
                    ComplexNumber derComplex = myPolynomial.Derivative().Evaluate(approximationComplex);
                    if (ComplexNumber.Round(derComplex, 6) == 0)
                    {
                        approximationComplex += ComplexNumber.i;
                        continue;
                    }

                    if (ComplexNumber.Round(resComplex, 6) == 0)
                    {
                        myPolynomial.DivideByRoot(approximationComplex, true);
                        if (approximationComplex.IsReal)
                        {
                            roots.Add((double)ComplexNumber.Round(approximationComplex, 6));
                        }
                        stuck = false;
                        break;
                    }

                    approximationComplex -= (resComplex / derComplex);
                }

                if (stuck)
                {
                    approximationComplex += 1;
                }
            }

            if (myPolynomial.degree == 1)
            {
                roots.Add(Math.Round(-(myPolynomial.values[0] / myPolynomial.values[1]), 6));
            }

            return roots.ToArray();
        }

        /// <summary>
        /// Returns a list of all roots of the polynomial using Newton's Method (Note: The same roots can be listed multiple times depending on their respective multiplicities.)
        /// </summary>
        /// <param name="knownRoots">Any roots (real or complex) that do not need to be found with Newton's Method (Note: For imaginary roots, do NOT list their conjugates here.)</param>
        public ComplexNumber[] FindComplexRoots(params ComplexNumber[] knownRoots)
        {
            List<ComplexNumber> roots = new List<ComplexNumber>();
            Polynomial myPolynomial = (Polynomial)MemberwiseClone();
            foreach (ComplexNumber root in knownRoots)
            {
                roots.Add(root);
                if (!root.IsReal)
                {
                    roots.Add(root.Conjugate);
                }
                myPolynomial.DivideByRoot(root);
            }

            while (myPolynomial.values[0] == 0)
            {
                roots.Add(0);
                myPolynomial.DivideByRoot(0);
            }

            while (myPolynomial.degree > 1)
            {
                bool stuck = true;
                ComplexNumber approximationComplex = ComplexNumber.i;
                for (int i = 0; i < 100; i++)
                {
                    ComplexNumber resComplex = myPolynomial.Evaluate(approximationComplex);
                    ComplexNumber derComplex = myPolynomial.Derivative().Evaluate(approximationComplex);
                    if (ComplexNumber.Round(derComplex, 6) == 0)
                    {
                        approximationComplex += ComplexNumber.i;
                        continue;
                    }

                    if (ComplexNumber.Round(resComplex, 6) == 0)
                    {
                        myPolynomial.DivideByRoot(approximationComplex, true);
                        roots.Add(ComplexNumber.Round(approximationComplex, 6));
                        if (!approximationComplex.IsReal)
                        {
                            roots.Add(ComplexNumber.Round(approximationComplex, 6).Conjugate);
                        }
                        stuck = false;
                        break;
                    }

                    approximationComplex -= (resComplex / derComplex);
                }

                if (stuck)
                {
                    approximationComplex += 1;
                }
            }

            if (myPolynomial.degree == 1)
            {
                roots.Add(Math.Round(-(myPolynomial.values[0] / myPolynomial.values[1]), 6));
            }

            return roots.ToArray();
        }

        /// <summary>
        /// Converts the polynomial to a vector of coefficients
        /// </summary>
        public Vector ToVector()
        {
            return new Vector(values);
        }

        public override string ToString()
        {
            string str = values[0].ToString();
            for (int i = 1; i <= degree; i++)
            {
                if (values[i] >= 0)
                {
                    if (i == 1)
                    {
                        str += (" + " + (values[i] == 1 ? "" : values[i].ToString()) + "x");
                    }
                    else
                    {
                        str += (" + " + (values[i] == 1 ? "" : values[i].ToString()) + "x^" + i);
                    }
                }
                else if (values[i] < 0)
                {
                    if (i == 1)
                    {
                        str += (" - " + (values[i] == -1 ? "" : (-values[i]).ToString()) + "x");
                    }
                    else
                    {
                        str += (" - " + (values[i] == -1 ? "" : (-values[i]).ToString()) + "x^" + i);
                    }
                }
            }

            return str;
        }
    }

    /// <summary>
    /// A complex number consisting of a real and imaginary part
    /// </summary>
    public struct ComplexNumber
    {
        /// <summary>
        /// The real part of the complex number
        /// </summary>
        public readonly double a;

        /// <summary>
        /// The imaginary part of the complex number
        /// </summary>
        public readonly double b;

        /// <summary>
        /// The magnitude of the complex number
        /// </summary>
        public readonly double r;

        /// <summary>
        /// The angle with the x-axis in the complex plane from -π [-pi] to π [pi], inclusive (expressed in radians)
        /// </summary>
        public readonly double theta;

        /// <summary>
        /// Initializes a complex number with two parts
        /// </summary>
        /// <param name="a">The real part of the number if initializing in rectangular form or the magnitude of the number if initializing in polar form</param>
        /// <param name="b">The imaginary part of the number if initializing in rectangular form or the angle (in radians) of the number if initializing in polar form</param>
        /// <param name="polar">Whether to initializes the number in polar or rectangular form form</param>
        public ComplexNumber(double a, double b, bool polar = false)
        {
            if (polar)
            {
                this.a = a * Math.Cos(b);
                this.b = a * Math.Sin(b);
                this.r = a;
                if (b < -Math.PI)
                {
                    b += (2 * Math.PI * Math.Ceiling(-b / (Math.PI)));
                }
                else if (b > Math.PI)
                {
                    b -= (2 * Math.PI * Math.Floor(b / (Math.PI)));
                }

                this.theta = b;
            }
            else
            {
                this.a = a;
                this.b = b;
                this.r = Math.Sqrt(a * a + b * b);
                if (r == 0)
                {
                    this.theta = 0;
                }
                else if (b >= 0)
                {
                    this.theta = Math.Acos(a / r);
                }
                else
                {
                    this.theta = (-Math.Acos(a / r));
                }
            }
        }

        /// <summary>
        /// Represents the number 0 as a complex number
        /// </summary>
        public static ComplexNumber Zero = new ComplexNumber(0, 0);

        /// <summary>
        /// Represents the number i as a complex number
        /// </summary>
        public static ComplexNumber i = new ComplexNumber(0, 1);

        /// <summary>
        /// Returns the square root of a real number expressed as a complex number
        /// </summary>
        public static ComplexNumber Sqrt(double a)
        {
            if (a < 0)
            {
                return new ComplexNumber(0, Math.Sqrt(-a));
            }
            else
            {
                return new ComplexNumber(Math.Sqrt(a), 0);
            }
        }

        /// <summary>
        /// Returns the square root of a complex number (where the angle is reduced in half)
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static ComplexNumber Sqrt(ComplexNumber a)
        {
            return new ComplexNumber(Math.Sqrt(a.r), a.theta / 2, true);
        }

        /// <summary>
        /// Raises a complex number to a specified real power
        /// </summary>
        /// <param name="a">The complex base</param>
        /// <param name="b">The real exponent</param>
        /// <returns></returns>
        public static ComplexNumber Pow(ComplexNumber a, double b)
        {
            return new ComplexNumber(Math.Pow(a.r, b), b * a.theta, true);
        }

        /// <summary>
        /// Rounds a complex number to a specified number of digits
        /// </summary>
        /// <param name="value">The complex number to be rounded</param>
        /// <param name="digits">The number of digits</param>
        public static ComplexNumber Round(ComplexNumber value, int digits)
        {
            return new ComplexNumber(Math.Round(value.a, digits), Math.Round(value.b, digits));
        }

        /// <summary>
        /// Flips the sign on both the real and imaginary parts of the complex number
        /// </summary>
        public static ComplexNumber operator -(ComplexNumber a)
        {
            return new ComplexNumber(-a.a, -a.b);
        }

        /// <summary>
        /// Returns a complex number representing the sum of two other complex numbers
        /// </summary>
        public static ComplexNumber operator +(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.a + b.a, a.b + b.b);
        }

        /// <summary>
        /// Returns a complex number representing the sum of a real and a complex number
        /// </summary>
        public static ComplexNumber operator +(double a, ComplexNumber b)
        {
            return new ComplexNumber(a + b.a, b.b);
        }

        /// <summary>
        /// Returns a complex number representing the sum of a real and a complex number
        /// </summary>
        public static ComplexNumber operator +(ComplexNumber a, double b)
        {
            return new ComplexNumber(a.a + b, a.b);
        }

        /// <summary>
        /// Returns a complex number representing the difference between two complex numbers
        /// </summary>
        public static ComplexNumber operator -(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.a - b.a, a.b - b.b);
        }

        /// <summary>
        /// Returns a complex number representing the difference between a real and a complex number
        /// </summary>
        public static ComplexNumber operator -(double a, ComplexNumber b)
        {
            return new ComplexNumber(a - b.a, -b.b);
        }

        /// <summary>
        /// Returns a complex number representing the difference between a complex and a real number
        /// </summary>
        public static ComplexNumber operator -(ComplexNumber a, double b)
        {
            return new ComplexNumber(a.a - b, a.b);
        }

        /// <summary>
        /// Returns a complex number representing the product of a complex number and a real number
        /// </summary>
        public static ComplexNumber operator *(ComplexNumber a, double b)
        {
            return new ComplexNumber(b * a.a, b * a.b);
        }

        /// <summary>
        /// Returns a complex number representing the product of a complex number and a real number
        /// </summary>
        public static ComplexNumber operator *(double a, ComplexNumber b)
        {
            return new ComplexNumber(a * b.a, a * b.b);
        }

        /// <summary>
        /// Returns a complex number representing the product of two complex numbers
        /// </summary>
        public static ComplexNumber operator *(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber((a.a * b.a) - (a.b * b.b), (a.a * b.b + a.b * b.a));
        }

        /// <summary>
        /// Returns a complex number representing the quotient of a complex number and a real number
        /// </summary>
        public static ComplexNumber operator /(ComplexNumber a, double b)
        {
            return new ComplexNumber(a.a / b, a.b / b);
        }

        /// <summary>
        /// Returns a complex number representing the quotient of a real number and a complex number
        /// </summary>
        public static ComplexNumber operator /(double a, ComplexNumber b)
        {
            return ((a * b.Conjugate) / (b.a * b.a + b.b * b.b));
        }

        /// <summary>
        /// Returns a complex number representing the quotient of two complex numbers
        /// </summary>
        public static ComplexNumber operator /(ComplexNumber a, ComplexNumber b)
        {
            return ((a * b.Conjugate) / (b.a * b.a + b.b * b.b));
        }

        /// <summary>
        /// Determines whether two complex numbers are equal
        /// </summary>
        public static bool operator ==(ComplexNumber a, ComplexNumber b)
        {
            return a.Equals(b);
        }

        public static bool operator !=(ComplexNumber a, ComplexNumber b)
        {
            return !a.Equals(b);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            ComplexNumber b = (ComplexNumber)obj;
            return (this.a == b.a && this.b == b.b);
        }

        /// <summary>
        /// Converts a real number to a complex number with an imaginary part of 0
        /// </summary>
        public static implicit operator ComplexNumber(double a) => new ComplexNumber(a, 0);

        /// <summary>
        /// Converts a complex number to a real number by dropping the imaginary part 
        /// </summary>
        public static explicit operator double(ComplexNumber a) => a.a;

        /// <summary>
        /// Returns the conjugate of the complex number
        /// </summary>
        public ComplexNumber Conjugate
        {
            get
            {
                return new ComplexNumber(a, -b);
            }
        }

        /// <summary>
        /// Determines whether the complex number is also a real number
        /// </summary>
        public bool IsReal { get { return Math.Round(b, 6) == 0; } }

        /// <summary>
        /// Raises the complex number to a real power
        /// </summary>
        /// <param name="a">The power, or exponent</param>
        public ComplexNumber Pow(double a)
        {
            return new ComplexNumber(Math.Pow(r, a), a * theta, true);
        }

        public override string ToString()
        {
            if (b < 0)
            {
                return a + " - " + (-b) + "i";
            }
            else if (b > 0)
            {
                return a + " + " + b + "i";
            }
            else
            {
                return a.ToString();
            }
        }
    }

    /// <summary>
    /// The main type of arithmetic exception thrown by linear algebra methods
    /// </summary>
    public class LinearAlgebraException : Exception
    {
        public LinearAlgebraException() { }
        public LinearAlgebraException(string message) : base(message) { }
        public LinearAlgebraException(string message, Exception inner) : base(message, inner) { }
    }

    /// <summary>
    /// An extension class providing a clamp method for any IComparable data type
    /// </summary>
    public static class LinearAlgebraExtensions
    {
        /// <summary>
        /// Constrains an object between a minimum and maximum value based on its comparison method
        /// </summary>
        public static T Clamp<T>(this T val, T min, T max) where T : IComparable<T>
        {
            if (val.CompareTo(min) < 0)
            {
                return min;
            }
            else if (val.CompareTo(max) > 0)
            {
                return max;
            }
            else
            {
                return val;
            }
        }
    }
}