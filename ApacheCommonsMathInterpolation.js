//http://www.source-code.biz/snippets/typescript/akima/ApacheCommonsMathInterpolation.ts
// This is a TypeScript version of some interpolation Java classes of the
// Apache Commons Math library. It currently contains the following main
// classes and related classes and interfaces:
//
//  - AkimaSplineInterpolator
//  - SplineInterpolator
//  - LinearInterpolator
//
// The source code has been manually translated from Java to Typescript.
// Unneeded class methods have been omitted.
//
// Copyright 2016 Christian d'Heureuse, Inventec Informatik AG, Zurich, Switzerland
// www.source-code.biz, www.inventec.ch/chdh
//
// -----------------------------------------------------------------------------
//
// Copyright 2001-2016 The Apache Software Foundation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
var ApacheCommonsMathInterpolation;
(function (ApacheCommonsMathInterpolation) {
    var EPSILON = 2.2204460492503130808472633361816E-16; // (with ES6 we could use Number.EPSILON)
    /**
     * Constant function.
     */
    var Constant = (function () {
        /**
         * @param c Constant.
         */
        function Constant(c) {
            this.c = c;
        }
        Constant.prototype.value = function (x) {
            return this.c;
        };
        return Constant;
    }());
    ApacheCommonsMathInterpolation.Constant = Constant;
    //------------------------------------------------------------------------------
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
    var AkimaSplineInterpolator = (function () {
        function AkimaSplineInterpolator() {
        }
        /**
         * Computes an interpolating function for the data set.
         *
         * @param xvals the arguments for the interpolation points
         * @param yvals the values for the interpolation points
         * @return a function which interpolates the data set
         */
        AkimaSplineInterpolator.prototype.interpolate = function (xvals, yvals) {
            if (!xvals || !yvals) {
                throw new Error("Null argument.");
            }
            if (xvals.length != yvals.length) {
                throw new Error("Dimension mismatch for xvals and yvals.");
            }
            if (xvals.length < AkimaSplineInterpolator.MINIMUM_NUMBER_POINTS) {
                throw new Error("Number of points is too small.");
            }
            MathArrays.checkOrder(xvals);
            var numberOfDiffAndWeightElements = xvals.length - 1;
            var differences = Array(numberOfDiffAndWeightElements);
            var weights = Array(numberOfDiffAndWeightElements);
            for (var i = 0; i < differences.length; i++) {
                differences[i] = (yvals[i + 1] - yvals[i]) / (xvals[i + 1] - xvals[i]);
            }
            for (var i = 1; i < weights.length; i++) {
                weights[i] = Math.abs(differences[i] - differences[i - 1]);
            }
            // Prepare Hermite interpolation scheme.
            var firstDerivatives = Array(xvals.length);
            for (var i = 2; i < firstDerivatives.length - 2; i++) {
                var wP = weights[i + 1];
                var wM = weights[i - 1];
                if (Math.abs(wP) < EPSILON && Math.abs(wM) < EPSILON) {
                    var xv = xvals[i];
                    var xvP = xvals[i + 1];
                    var xvM = xvals[i - 1];
                    firstDerivatives[i] = (((xvP - xv) * differences[i - 1]) + ((xv - xvM) * differences[i])) / (xvP - xvM);
                }
                else {
                    firstDerivatives[i] = ((wP * differences[i - 1]) + (wM * differences[i])) / (wP + wM);
                }
            }
            firstDerivatives[0] = this.differentiateThreePoint(xvals, yvals, 0, 0, 1, 2);
            firstDerivatives[1] = this.differentiateThreePoint(xvals, yvals, 1, 0, 1, 2);
            firstDerivatives[xvals.length - 2] = this.differentiateThreePoint(xvals, yvals, xvals.length - 2, xvals.length - 3, xvals.length - 2, xvals.length - 1);
            firstDerivatives[xvals.length - 1] = this.differentiateThreePoint(xvals, yvals, xvals.length - 1, xvals.length - 3, xvals.length - 2, xvals.length - 1);
            return this.interpolateHermiteSorted(xvals, yvals, firstDerivatives);
        };
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
        AkimaSplineInterpolator.prototype.differentiateThreePoint = function (xvals, yvals, indexOfDifferentiation, indexOfFirstSample, indexOfSecondsample, indexOfThirdSample) {
            var x0 = yvals[indexOfFirstSample];
            var x1 = yvals[indexOfSecondsample];
            var x2 = yvals[indexOfThirdSample];
            var t = xvals[indexOfDifferentiation] - xvals[indexOfFirstSample];
            var t1 = xvals[indexOfSecondsample] - xvals[indexOfFirstSample];
            var t2 = xvals[indexOfThirdSample] - xvals[indexOfFirstSample];
            var a = (x2 - x0 - (t2 / t1 * (x1 - x0))) / (t2 * t2 - t1 * t2);
            var b = (x1 - x0 - a * t1 * t1) / t1;
            return (2 * a * t) + b;
        };
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
        AkimaSplineInterpolator.prototype.interpolateHermiteSorted = function (xvals, yvals, firstDerivatives) {
            if (xvals.length != yvals.length) {
                throw new Error("Dimension mismatch");
            }
            if (xvals.length != firstDerivatives.length) {
                throw new Error("Dimension mismatch");
            }
            var minimumLength = 2;
            if (xvals.length < minimumLength) {
                throw new Error("Not enough points.");
            }
            var size = xvals.length - 1;
            var polynomials = Array(size);
            var coefficients = Array(4);
            for (var i = 0; i < polynomials.length; i++) {
                var w = xvals[i + 1] - xvals[i];
                var w2 = w * w;
                var yv = yvals[i];
                var yvP = yvals[i + 1];
                var fd = firstDerivatives[i];
                var fdP = firstDerivatives[i + 1];
                coefficients[0] = yv;
                coefficients[1] = firstDerivatives[i];
                coefficients[2] = (3 * (yvP - yv) / w - 2 * fd - fdP) / w;
                coefficients[3] = (2 * (yv - yvP) / w + fd + fdP) / w2;
                polynomials[i] = new PolynomialFunction(coefficients);
            }
            return new PolynomialSplineFunction(xvals, polynomials);
        };
        return AkimaSplineInterpolator;
    }()); // end class AkimaSplineInterpolator
    /** The minimum number of points that are needed to compute the function. */
    AkimaSplineInterpolator.MINIMUM_NUMBER_POINTS = 5;
    ApacheCommonsMathInterpolation.AkimaSplineInterpolator = AkimaSplineInterpolator;
    //------------------------------------------------------------------------------
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
    var SplineInterpolator = (function () {
        function SplineInterpolator() {
        }
        /**
         * Computes an interpolating function for the data set.
         * @param x the arguments for the interpolation points
         * @param y the values for the interpolation points
         * @return a function which interpolates the data set
         */
        SplineInterpolator.prototype.interpolate = function (x, y) {
            if (x.length != y.length) {
                throw new Error("Dimension msmatch.");
            }
            if (x.length < 3) {
                throw new Error("Number of points is too small.");
            }
            // Number of intervals.  The number of data points is n + 1.
            var n = x.length - 1;
            MathArrays.checkOrder(x);
            // Differences between knot points
            var h = Array(n);
            for (var i = 0; i < n; i++) {
                h[i] = x[i + 1] - x[i];
            }
            var mu = Array(n);
            var z = Array(n + 1);
            mu[0] = 0;
            z[0] = 0;
            var g = 0;
            for (var i = 1; i < n; i++) {
                g = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / g;
                z[i] = (3 * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1]) + y[i - 1] * h[i]) /
                    (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / g;
            }
            // cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
            var b = Array(n);
            var c = Array(n + 1);
            var d = Array(n);
            z[n] = 0;
            c[n] = 0;
            for (var j = n - 1; j >= 0; j--) {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
                d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
            }
            var polynomials = Array(n);
            var coefficients = Array(4);
            for (var i = 0; i < n; i++) {
                coefficients[0] = y[i];
                coefficients[1] = b[i];
                coefficients[2] = c[i];
                coefficients[3] = d[i];
                polynomials[i] = new PolynomialFunction(coefficients);
            }
            return new PolynomialSplineFunction(x, polynomials);
        };
        return SplineInterpolator;
    }()); // end class SplineInterpolator
    ApacheCommonsMathInterpolation.SplineInterpolator = SplineInterpolator;
    //------------------------------------------------------------------------------
    /**
     * Implements a linear function for interpolation of real univariate functions.
     */
    var LinearInterpolator = (function () {
        function LinearInterpolator() {
        }
        /**
         * Computes a linear interpolating function for the data set.
         *
         * @param x the arguments for the interpolation points
         * @param y the values for the interpolation points
         * @return a function which interpolates the data set
         */
        LinearInterpolator.prototype.interpolate = function (x, y) {
            if (x.length != y.length) {
                throw new Error("Dimension msmatch.");
            }
            if (x.length < 2) {
                throw new Error("Number of points is too small.");
            }
            // Number of intervals.  The number of data points is n + 1.
            var n = x.length - 1;
            MathArrays.checkOrder(x);
            // Slope of the lines between the datapoints.
            var m = Array(n);
            for (var i = 0; i < n; i++) {
                m[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            }
            var polynomials = Array(n);
            var coefficients = Array(2);
            for (var i = 0; i < n; i++) {
                coefficients[0] = y[i];
                coefficients[1] = m[i];
                polynomials[i] = new PolynomialFunction(coefficients);
            }
            return new PolynomialSplineFunction(x, polynomials);
        };
        return LinearInterpolator;
    }()); // end class LinearInterpolator
    ApacheCommonsMathInterpolation.LinearInterpolator = LinearInterpolator;
    //------------------------------------------------------------------------------
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
    var PolynomialSplineFunction = (function () {
        /**
         * Construct a polynomial spline function with the given segment delimiters
         * and interpolating polynomials.
         * The constructor copies both arrays and assigns the copies to the knots
         * and polynomials properties, respectively.
         *
         * @param knots Spline segment interval delimiters.
         * @param polynomials Polynomial functions that make up the spline.
         */
        function PolynomialSplineFunction(knots, polynomials) {
            if (!knots || !polynomials) {
                throw new Error("Null argument.");
            }
            if (knots.length < 2) {
                throw new Error("Illegal argument knots.");
            }
            if (knots.length - 1 != polynomials.length) {
                throw new Error("Dimension mismatch.");
            }
            MathArrays.checkOrder(knots);
            this.n = knots.length - 1;
            this.knots = knots.slice();
            this.polynomials = polynomials.slice();
        }
        /**
         * Compute the value for the function.
         * See {@link PolynomialSplineFunction} for details on the algorithm for
         * computing the value of the function.
         *
         * @param v Point for which the function value should be computed.
         * @return the value.
         */
        PolynomialSplineFunction.prototype.value = function (v) {
            // Patched to allow an extended function domain.
            //      if (v < this.knots[0] || v > this.knots[this.n]) {
            //          throw new Error("Argument value out of range.");
            //      }
            var i = Arrays.binarySearch(this.knots, v);
            if (i < 0) {
                i = -i - 2;
            }
            // This will handle the case where v is the last knot value
            // There are only n-1 polynomials, so if v is the last knot
            // then we will use the last polynomial to calculate the value.
            // Patched to allow an extended function domain.
            //      if ( i >= this.polynomials.length ) {
            //          i--;
            //      }
            i = Math.max(0, Math.min(i, this.polynomials.length - 1));
            return this.polynomials[i].value(v - this.knots[i]);
        };
        return PolynomialSplineFunction;
    }()); // end class PolynomialSplineFunction
    ApacheCommonsMathInterpolation.PolynomialSplineFunction = PolynomialSplineFunction;
    //------------------------------------------------------------------------------
    /**
     * Immutable representation of a real polynomial function with real coefficients.
     * <p>
     * <a href="http://mathworld.wolfram.com/HornersMethod.html">Horner's Method</a>
     * is used to evaluate the function.</p>
     */
    var PolynomialFunction = (function () {
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
        function PolynomialFunction(c) {
            var n = c.length;
            if (n == 0) {
                throw new Error("Empty polynomials coefficients array");
            }
            while ((n > 1) && (c[n - 1] == 0)) {
                --n;
            }
            this.coefficients = c.slice();
        }
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
        PolynomialFunction.prototype.value = function (x) {
            return PolynomialFunction.evaluate(this.coefficients, x);
        };
        /**
         * Uses Horner's Method to evaluate the polynomial with the given coefficients at
         * the argument.
         *
         * @param coefficients Coefficients of the polynomial to evaluate.
         * @param argument Input value.
         * @return the value of the polynomial.
         */
        PolynomialFunction.evaluate = function (coefficients, argument) {
            var n = coefficients.length;
            if (n == 0) {
                throw new Error("Empty polynomials coefficients array.");
            }
            var result = coefficients[n - 1];
            for (var j = n - 2; j >= 0; j--) {
                result = argument * result + coefficients[j];
            }
            return result;
        };
        return PolynomialFunction;
    }()); // end class PolynomialFunction
    ApacheCommonsMathInterpolation.PolynomialFunction = PolynomialFunction;
    //------------------------------------------------------------------------------
    /**
     * Arrays utilities.
     */
    var MathArrays = (function () {
        function MathArrays() {
        }
        /**
         * Check that the given array is sorted in strictly increasing order.
         *
         * @param val Values.
         */
        MathArrays.checkOrder = function (val) {
            var previous = val[0];
            var max = val.length;
            for (var index = 1; index < max; index++) {
                if (val[index] <= previous) {
                    throw new Error("Non-monotonic sequence exception.");
                }
                previous = val[index];
            }
        };
        return MathArrays;
    }()); // end class MathArrays
    //------------------------------------------------------------------------------
    // Corresponds to java.util.Arrays.
    var Arrays = (function () {
        function Arrays() {
        }
        Arrays.binarySearch = function (a, key) {
            var low = 0;
            var high = a.length - 1;
            while (low <= high) {
                var mid = (low + high) >>> 1;
                var midVal = a[mid];
                if (midVal < key)
                    low = mid + 1;
                else if (midVal > key)
                    high = mid - 1;
                else if (midVal == key)
                    return mid;
                else {
                    throw new Error("Invalid number encountered in binary search.");
                }
            }
            return -(low + 1); // key not found.
        };
        return Arrays;
    }()); // end class Arrays
})(ApacheCommonsMathInterpolation || (ApacheCommonsMathInterpolation = {})); // end namespace ApacheCommonsMath
