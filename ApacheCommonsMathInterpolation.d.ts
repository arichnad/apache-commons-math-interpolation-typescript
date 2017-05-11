declare namespace ApacheCommonsMathInterpolation {
    /**
     * An interface representing a univariate real function.
     */
    interface UnivariateFunction {
        /**
         * Compute the value of the function.
         *
         * @param x Point at which the function value should be computed.
         * @return the value of the function.
         */
        value(x: number): number;
    }
    /**
     * Interface representing a univariate real interpolating function.
     */
    interface UnivariateInterpolator {
        /**
         * Compute an interpolating function for the dataset.
         *
         * @param xval Arguments for the interpolation points.
         * @param yval Values for the interpolation points.
         * @return a function which interpolates the dataset.
         */
        interpolate(xvals: number[], yvals: number[]): UnivariateFunction;
    }
    /**
     * Constant function.
     */
    class Constant implements UnivariateFunction {
        private readonly c;
        /**
         * @param c Constant.
         */
        constructor(c: number);
        value(x: number): number;
    }
    /**
     * Computes a cubic spline interpolation for the data set using the Akima
     * algorithm, as originally formulated by Hiroshi Akima in his 1970 paper
     * "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures."
     * J. ACM 17, 4 (October 1970), 589-602. DOI=10.1145/321607.321609
     * http://doi.acm.org/10.1145/321607.321609
     * <p>
     * This implementation is based on the Akima implementation in the CubicSpline
     * class in the Math.NET Numerics library. The method referenced is
     * CubicSpline.InterpolateAkimaSorted
     * </p>
     * <p>
     * The {@link #interpolate(double[], double[]) interpolate} method returns a
     * {@link PolynomialSplineFunction} consisting of n cubic polynomials, defined
     * over the subintervals determined by the x values, {@code x[0] < x[i] ... < x[n]}.
     * The Akima algorithm requires that {@code n >= 5}.
     * </p>
     */
    class AkimaSplineInterpolator implements UnivariateInterpolator {
        /** The minimum number of points that are needed to compute the function. */
        private static readonly MINIMUM_NUMBER_POINTS;
        /**
         * Computes an interpolating function for the data set.
         *
         * @param xvals the arguments for the interpolation points
         * @param yvals the values for the interpolation points
         * @return a function which interpolates the data set
         */
        interpolate(xvals: number[], yvals: number[]): PolynomialSplineFunction;
        /**
         * Three point differentiation helper, modeled off of the same method in the
         * Math.NET CubicSpline class. This is used by both the Apache Math and the
         * Math.NET Akima Cubic Spline algorithms
         *
         * @param xvals x values to calculate the numerical derivative with
         * @param yvals y values to calculate the numerical derivative with
         * @param indexOfDifferentiation index of the elemnt we are calculating the derivative around
         * @param indexOfFirstSample index of the first element to sample for the three point method
         * @param indexOfSecondsample index of the second element to sample for the three point method
         * @param indexOfThirdSample index of the third element to sample for the three point method
         * @return the derivative
         */
        private differentiateThreePoint(xvals, yvals, indexOfDifferentiation, indexOfFirstSample, indexOfSecondsample, indexOfThirdSample);
        /**
         * Creates a Hermite cubic spline interpolation from the set of (x,y) value
         * pairs and their derivatives. This is modeled off of the
         * InterpolateHermiteSorted method in the Math.NET CubicSpline class.
         *
         * @param xvals x values for interpolation
         * @param yvals y values for interpolation
         * @param firstDerivatives first derivative values of the function
         * @return polynomial that fits the function
         */
        private interpolateHermiteSorted(xvals, yvals, firstDerivatives);
    }
    /**
     * Computes a natural (also known as "free", "unclamped") cubic spline interpolation for the data set.
     * <p>
     * The {@link #interpolate(double[], double[])} method returns a {@link PolynomialSplineFunction}
     * consisting of n cubic polynomials, defined over the subintervals determined by the x values,
     * {@code x[0] < x[i] ... < x[n].}  The x values are referred to as "knot points."</p>
     * <p>
     * The value of the PolynomialSplineFunction at a point x that is greater than or equal to the smallest
     * knot point and strictly less than the largest knot point is computed by finding the subinterval to which
     * x belongs and computing the value of the corresponding polynomial at <code>x - x[i] </code> where
     * <code>i</code> is the index of the subinterval.  See {@link PolynomialSplineFunction} for more details.
     * </p>
     * <p>
     * The interpolating polynomials satisfy: <ol>
     * <li>The value of the PolynomialSplineFunction at each of the input x values equals the
     *  corresponding y value.</li>
     * <li>Adjacent polynomials are equal through two derivatives at the knot points (i.e., adjacent polynomials
     *  "match up" at the knot points, as do their first and second derivatives).</li>
     * </ol></p>
     * <p>
     * The cubic spline interpolation algorithm implemented is as described in R.L. Burden, J.D. Faires,
     * <u>Numerical Analysis</u>, 4th Ed., 1989, PWS-Kent, ISBN 0-53491-585-X, pp 126-131.
     * </p>
     *
     */
    class SplineInterpolator implements UnivariateInterpolator {
        /**
         * Computes an interpolating function for the data set.
         * @param x the arguments for the interpolation points
         * @param y the values for the interpolation points
         * @return a function which interpolates the data set
         */
        interpolate(x: number[], y: number[]): PolynomialSplineFunction;
    }
    /**
     * Implements a linear function for interpolation of real univariate functions.
     */
    class LinearInterpolator implements UnivariateInterpolator {
        /**
         * Computes a linear interpolating function for the data set.
         *
         * @param x the arguments for the interpolation points
         * @param y the values for the interpolation points
         * @return a function which interpolates the data set
         */
        interpolate(x: number[], y: number[]): PolynomialSplineFunction;
    }
    /**
     * Represents a polynomial spline function.
     * <p>
     * A <strong>polynomial spline function</strong> consists of a set of
     * <i>interpolating polynomials</i> and an ascending array of domain
     * <i>knot points</i>, determining the intervals over which the spline function
     * is defined by the constituent polynomials.  The polynomials are assumed to
     * have been computed to match the values of another function at the knot
     * points.  The value consistency constraints are not currently enforced by
     * <code>PolynomialSplineFunction</code> itself, but are assumed to hold among
     * the polynomials and knot points passed to the constructor.</p>
     * <p>
     * N.B.:  The polynomials in the <code>polynomials</code> property must be
     * centered on the knot points to compute the spline function values.
     * See below.</p>
     * <p>
     * The domain of the polynomial spline function is
     * <code>[smallest knot, largest knot]</code>.  Attempts to evaluate the
     * function at values outside of this range generate IllegalArgumentExceptions.
     * </p>
     * <p>
     * The value of the polynomial spline function for an argument <code>x</code>
     * is computed as follows:
     * <ol>
     * <li>The knot array is searched to find the segment to which <code>x</code>
     * belongs.  If <code>x</code> is less than the smallest knot point or greater
     * than the largest one, an <code>IllegalArgumentException</code>
     * is thrown.</li>
     * <li> Let <code>j</code> be the index of the largest knot point that is less
     * than or equal to <code>x</code>.  The value returned is
     * {@code polynomials[j](x - knot[j])}</li></ol>
     */
    class PolynomialSplineFunction {
        /**
         * Spline segment interval delimiters (knots).
         * Size is n + 1 for n segments.
         */
        private knots;
        /**
         * The polynomial functions that make up the spline.  The first element
         * determines the value of the spline over the first subinterval, the
         * second over the second, etc.   Spline function values are determined by
         * evaluating these functions at {@code (x - knot[i])} where i is the
         * knot segment to which x belongs.
         */
        private polynomials;
        /**
         * Number of spline segments. It is equal to the number of polynomials and
         * to the number of partition points - 1.
         */
        private n;
        /**
         * Construct a polynomial spline function with the given segment delimiters
         * and interpolating polynomials.
         * The constructor copies both arrays and assigns the copies to the knots
         * and polynomials properties, respectively.
         *
         * @param knots Spline segment interval delimiters.
         * @param polynomials Polynomial functions that make up the spline.
         */
        constructor(knots: number[], polynomials: PolynomialFunction[]);
        /**
         * Compute the value for the function.
         * See {@link PolynomialSplineFunction} for details on the algorithm for
         * computing the value of the function.
         *
         * @param v Point for which the function value should be computed.
         * @return the value.
         */
        value(v: number): number;
    }
    /**
     * Immutable representation of a real polynomial function with real coefficients.
     * <p>
     * <a href="http://mathworld.wolfram.com/HornersMethod.html">Horner's Method</a>
     * is used to evaluate the function.</p>
     */
    class PolynomialFunction {
        /**
         * The coefficients of the polynomial, ordered by degree -- i.e.,
         * coefficients[0] is the constant term and coefficients[n] is the
         * coefficient of x^n where n is the degree of the polynomial.
         */
        private coefficients;
        /**
         * Construct a polynomial with the given coefficients.  The first element
         * of the coefficients array is the constant term.  Higher degree
         * coefficients follow in sequence.  The degree of the resulting polynomial
         * is the index of the last non-null element of the array, or 0 if all elements
         * are null.
         * <p>
         * The constructor makes a copy of the input array and assigns the copy to
         * the coefficients property.</p>
         *
         * @param c Polynomial coefficients.
         */
        constructor(c: number[]);
        /**
         * Compute the value of the function for the given argument.
         * <p>
         *  The value returned is </p><p>
         *  {@code coefficients[n] * x^n + ... + coefficients[1] * x  + coefficients[0]}
         * </p>
         *
         * @param x Argument for which the function value should be computed.
         * @return the value of the polynomial at the given point.
         *
         * @see org.hipparchus.analysis.UnivariateFunction#value(double)
         */
        value(x: number): number;
        /**
         * Uses Horner's Method to evaluate the polynomial with the given coefficients at
         * the argument.
         *
         * @param coefficients Coefficients of the polynomial to evaluate.
         * @param argument Input value.
         * @return the value of the polynomial.
         */
        protected static evaluate(coefficients: number[], argument: number): number;
    }
}
