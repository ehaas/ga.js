(function() {
    "use strict";
    var util = require('./util/bits');
    var InnerProductTypes = require('./inner_product_types');
    var BasisBlade = require('./_basis_blade').BasisBlade;
    var Metric = require('./metric').Metric;
    var fun = require('./util/pattern_match').fun;
    var MathU = require('./util/math_util').MathU;
    function minusOnePow(i) {
        return ((i & 1) === 0) ? 1 : -1;
    }

    BasisBlade.fromScalar = function(s) {
        return new BasisBlade(0, s);
    }
    BasisBlade.fromBitmap = function(b) {
        return new BasisBlade(b, 1.0);
    }

    BasisBlade.prototype = {
        equals: fun(
            [BasisBlade, function(o) { this.bitmap === o.bitmap && this.scale === o.scale; }],
            [Number, function(n) { return this.bitmap === 0 && this.scale === n; }],
            [fun.Any, function() { return false; }]
        ),
        grade: function() {
            return util.Bits.bitCount(this.bitmap);
        },
        reverse: function() {
            // multiplier = (-1)^(x(x-1)/2)
            return new BasisBlade(this.bitmap, minusOnePow((this.grade() * (this.grade() - 1)) / 2) * this.scale);
        },
        gradeInversion: function() {
            // multiplier is (-1)^x
            return new BasisBlade(this.bitmap, minusOnePow(this.grade()) * this.scale);
        },
        cliffordConjugate: function() {
            // multiplier is ((-1)^(x(x+1)/2)
            return new BasisBlade(this.bitmap, minusOnePow((this.grade() * (this.grade() + 1)) / 2) * scale);
        },
        round: function(multipleOf, epsilon) {
            /**
             * Rounds the scalar part of this to the nearest multiple X of multipleOf,
             * if |X - what| <= epsilon. This is useful when eigenbasis is used to perform products in arbitrary
             * metric, which leads to small roundof errors. You don't want to keep these roundof errors if you
             * are computing a multiplication table.
             *
             * @returns a new basis blade if a change is required.
             */
            var a = util.DoubleU.round(this.scale, multipleOf, epsilon);
            if (a !== this.scale) {
                return new BasisBlade(this.bitmap, a);
            }
            return this;
        },
        clone: function() {
            return new BasisBlade(this.bitmap, this.scale);
        },
        toString: function(bvNames) {
            var result = '', i = 1, b = this.bitmap;
            while (b !== 0) {
                if ((b & 1) !== 0) {
                    if (result.length > 0) {
                        result += "^";
                    }
                    if ((bvNames === undefined) || (i > bvNames.length) || (bvNames[i - 1] === undefined)) {
                        result += ("e" + i);
                    } else {
                        result += bvNames[i - 1];
                    }
                }
                b >>= 1;
                i++;
            }

            return (result.length === 0) ? this.scale.toString() : this.scale + "*" + result;

        }
    };

    BasisBlade.canonicalReorderingSign = function(a, b) {
        a = a >>> 1;
        var sum = 0;
        while (a !== 0) {
            sum = sum + util.Bits.bitCount(a & b);
            a = a >>> 1;
        }
        return ((sum & 1) === 0) ? 1.0 : -1.0;
    }

    BasisBlade.op = BasisBlade.outerProduct = function(a, b) {
        return BasisBlade.gp_op(a, b, true);
    }

    BasisBlade.gp = BasisBlade.geometricProduct = fun(
        [BasisBlade, BasisBlade, function(a, b) {
            return BasisBlade.gp_op(a, b, false);
        }],
        [BasisBlade, BasisBlade, Array, function(a, b, m) {
            // compute the geometric product in Euclidean metric:
            var i, bitmap, result = BasisBlade.geometricProduct(a, b);

            // compute the meet (bitmap of annihilated vectors):
            bitmap = a.bitmap & b.bitmap;

            // change the scale according to the metric:
            i = 0;
            while (bitmap !== 0) {
                if ((bitmap & 1) !== 0) result.scale *= m[i];
                i++;
                bitmap = bitmap >> 1;
            }

            return result;
        }],
        [BasisBlade, BasisBlade, Metric, function(a, b, M) {
            // convert argument to eigenbasis
            var A = M.toEigenBasis(a);
            var B = M.toEigenBasis(b);

            var result = [];
            for (var i = 0; i < A.length; i++) {
                for (var j = 0; j < B.length; j++) {
                    result.push(BasisBlade.gp(A[i], B[j], M.eigenMetric));
                }
            }
            return M.toMetricBasis(BasisBlade.simplify(result));
        }]
    );

    BasisBlade.gp_op = function(a, b, outer) {
        if (outer && ((a.bitmap & b.bitmap) !== 0)) {
            return new BasisBlade(0, 0.0);
        }
        var bitmap = a.bitmap ^ b.bitmap;
        var sign = BasisBlade.canonicalReorderingSign(a.bitmap, b.bitmap);
        return new BasisBlade(bitmap, sign * a.scale * b.scale);
    }

    BasisBlade.ip = BasisBlade.innerProduct = fun(
        [BasisBlade, BasisBlade, Number, function(a, b, type) {
            return innerProductFilter(a.grade(), b.grade(), BasisBlade.geometricProduct(a, b), type);
        }],
        [BasisBlade, BasisBlade, Array, Number, function(a, b, m, type) {
            return innerProductFilter(a.grade(), b.grade(), BasisBlade.geometricProduct(a, b, m), type);
        }],
        [BasisBlade, BasisBlade, Metric, Number, function(a, b, M, type) {
            return innerProductFilter(a.grade(), b.grade(), BasisBlade.geometricProduct(a, b, M), type);
        }]
    );

    function bladesComparator(b1, b2) {
        if (b1.bitmap < b2.bitmap) return -1;
        else if (b1.bitmap > b2.bitmap) return 1;
        else return b1.scale - b2.scale;
    }

    BasisBlade.simplify = function(A) {
        if (A.length === 0) return A;

        A.sort(bladesComparator);
        var result = [];
        var current = A[0].clone();
        for (var i = 1; i < A.length; i++) {
            var b = A[i];
            if (b.bitmap === current.bitmap) {
                current.scale += b.scale;
            } else {
                if (current.scale !== 0) {
                    current.scale = MathU.normalize(current.scale);
                    result.push(current);
                }
                current = b.clone();
            }
        }
        if (current.scale !== 0) {
            current.scale = MathU.normalize(current.scale);
            result.push(current);
        }

        return result;
    };

    var innerProductFilter = fun(
        [Number, Number, Array, Number, function(ga, gb, R, type) {
            var result = [];
            for (var i = 0; i < R.length; i++) {
                var B = innerProductFilter(ga, gb, R[i], type);
                if (B.scale !== 0.0) {
                    result.push(B);
                }
            }
            return result;
        }],
        [Number, Number, BasisBlade, Number, function(ga, gb, r, type) {
            switch(type) {
                case InnerProductTypes.LEFT_CONTRACTION:
                    if ((ga > gb) || (r.grade() !== (gb-ga)))
                        return new BasisBlade();
                    else return r;
                case InnerProductTypes.RIGHT_CONTRACTION:
                    if ((ga < gb) || (r.grade() !== (ga-gb)))
                        return new BasisBlade();
                    else return r;
                case InnerProductTypes.HESTENES_INNER_PRODUCT:
                    if ((ga === 0) || (gb === 0)) return new BasisBlade();
                    // drop through to MODIFIED_HESTENES_INNER_PRODUCT
                case InnerProductTypes.MODIFIED_HESTENES_INNER_PRODUCT:
                    if (Math.abs(ga - gb) === r.grade()) return r;
                    else return new BasisBlade();
                default:
                    return null;
            }
        }]
    );

    exports.BasisBlade = BasisBlade;
})();
