// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

/*
 * Metric.java
 *
 * Created on February 1, 2005, 12:20 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */

/**
 * A class for representing totally arbitrary metric. Metric is given
 * in the form of a (symmetric) matrix. The eigenvectors of the matrix
 * are computed since these are used for computing (geometric) products
 * in this metric.
 *
 * @author  fontijne
 */
(function() {
    var fun = require('./util/pattern_match').fun;
    var Metric = require('./_metric').Metric;
    var BasisBlade = require('./_basis_blade').BasisBlade;
    Metric.prototype = {
        toEigenBasis: fun(
            [BasisBlade, function(a) { return this.transform(a, this.invEigMatrix); }],
            [Array, function(a) {
                var result = [];
                for (var i = 0; i < a.length; i++) {
                    result = result.concat(this.toEigenBasis(a[i]));
                }
                return BasisBlade.simplify(result);
            }]
        ),
        toMetricBasis: fun(
            [BasisBlade, function(a) { return this.transform(a, this.eig.matrix); }],
            [Array, function(a) {
                var result = [];
                for (var i = 0; i < a.length; i++) {
                    result = result.concat(this.toMetricBasis(a[i]));
                }
                return BasisBlade.simplify(result);
            }]
        ),
        transform: function(a, M) {
            var A = [];
            A.push(BasisBlade.fromScalar(a.scale)); // start with just scalar;

            // for each 1 bit: convert to list of blades
            var i = 0;
            var b = a.bitmap;
            while (b !== 0) {
                if ((b & 1) !== 0) {
                    // take column 'i' out of the matrix, wedge it to 'A'
                    var tmp = [];
                    for (var j = 0; j < M.length; j++) {
                        if (M[j][i] !== 0) {
                            var m = M[j][i];
                            for (var k = 0; k < A.length; k++) {
                                tmp.push(BasisBlade.op(A[k], new BasisBlade((1 << j), m)));
                            }
                        }
                    }
                    A = tmp;
                }

                b >>= 1;
                i++;
            }
            return A;
        }
    };

    exports.Metric = Metric;
})();


