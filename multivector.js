(function() {
    var Bits = require('./util/bits').Bits;
    var Metric = require('./metric').Metric;
    var BasisBlade = require('./basis_blade').BasisBlade;
    var InnerProductTypes = require('./inner_product_types');
    var MathU = require('./util/math_util').MathU;
    var numeric = require('./numeric');
    function slice(a) { return Array.prototype.slice.call(a, 0); };
    var fun = require('./util/pattern_match').fun;

    function cloneArrayElements(a) {
        for (var i = 0; i < a.length; i++) {
            a[i] = a[i].clone();
        }
        return a;
    }
    function getter(attr) {
        return function(o) { return o[attr]; };
    }

    var Multivector = fun(
        [Array, function(a) { this.blades = slice(a); }],
        [BasisBlade, function(b) { this.blades = [b.clone()] }],
        [Number, function(a) { return new Multivector(BasisBlade.fromScalar(a)); }],
        [function() { this.blades = []; }]
    );

    Multivector.prototype = {
        equals: fun(
            [Multivector, function(x) {
                var zero = this.subtract(x);
                return zero.blades.length === 0;
            }],
            [Number, function(a) {
                if (a === 0 && this.blades.length === 0) return true;

                var couldMatch = false, fail = false;
                this.blades.forEach(function(blade) {
                    if (blade.bitmap === 0 && blade.scale === a) {
                        couldMatch = true;
                    }
                    if ((blade.bitmap !== 0 && blade.scale !== 0.0)) {
                        fail = true;
                    }
                });
                return couldMatch && !fail;

            }],
            [fun.Any, function() { return false; }]
        ),
        gp: fun(
            [Multivector, function(x) {
                var result = [];
                for (var i = 0; i < this.blades.length; i++) {
                    for (var j = 0; j < x.blades.length; j++) {
                        result.push(BasisBlade.gp(this.blades[i], x.blades[j]));
                    }
                }
                return new Multivector(simplify(result));
            }],
            [Number, function(a) {
                if (a === 0.0) return new Multivector();
                var result = [];
                for (var i = 0; i < this.blades.length; i++) {
                    var b = this.blades[i];
                    result.push(new BasisBlade(b.bitmap, b.scale * a));
                }
                return new Multivector(result);
            }],
            [Multivector, Metric, function(x, M) {
                var result = [];
                for (var i = 0; i < this.blades.length; i++) {
                    for (var j = 0; j < x.blades.length; j++) {
                        result = result.concat(BasisBlade.gp(this.blades[i], x.blades[j], M));
                    }
                }
                return new Multivector(simplify(result));
            }]
        ),
        op: function(x) {
            var result = []
            for (var i = 0; i < this.blades.length ; i++) {
                var B1 = this.blades[i];
                for (var j = 0; j < x.blades.length; j++) {
                    var B2 = x.blades[j];
                    result.push(BasisBlade.op(B1, B2));
                }
            }
            return new Multivector(simplify(result));
        },
        simplify: function() {
            simplify(this.blades);
            return this;
        },
        add: fun(
            [Multivector, function(x) {
                var result = this.blades.concat(x.blades);
                return new Multivector(simplify(cloneArrayElements(result)));
            }],
            [Number, function(a) {
                var result = slice(this.blades);
                result.push(BasisBlade.fromScalar(a));
                return new Multivector(simplify(cloneArrayElements(result)));
            }]
        ),
        subtract: fun(
            [Multivector, function(x) {
                var result = this.blades.concat(x.gp(-1.0).blades);
                return new Multivector(simplify(cloneArrayElements(result)));
            }],
            [Number, function(a) { return this.add(-a); }]
        ),
        exp: fun(
            [function() { return this._exp(null, 12); }],
            [Metric, function(M) { return this._exp(M, 12); }]
        ),
        _exp: function(M, order) {
            // check out this^2 for special cases
            var A2 = this._gp(this, M).compress();
            if (A2.isNull(1e-8)) {
                // special case A^2 = 0
                return this.add(1);
            } else if (A2.isScalar()) {
                var a2 = A2.scalarPart();
                // special case A^2 = +-alpha^2
                if (a2 < 0) {
                    var alpha = Math.sqrt(-a2);
                    return this.gp(Math.sin(alpha) / alpha).add(Math.cos(alpha));
                } else {
                    //hey: todo what if a2 == 0?
                    var alpha = Math.sqrt(a2);
                    return this.gp(MathU.sinh(alpha) / alpha).add(MathU.cosh(alpha));
                }
            }
            return this.expSeries(M, order);
        },
        expSeries: function(M, order) {
            var scale = 1;
            var max = this.norm_e();
            if (max > 1.0) scale <<= 1;
            while (max > 1.0) {
                max = max / 2;
                scale <<= 1;
            }
            var scaled = this.gp(1.0 / scale);
            var result = new Multivector(1.0);
            var tmp = new Multivector(1.0);

            for (var i = 1; i < order; i++) {
                tmp = tmp._gp(scaled.gp(1.0 / i), M);
                result = result.add(tmp);
            }
            while (scale > 1) {
                result = result._gp(result, M);
                scale >>>= 1;
            }
            return result;
        },
        _gp: function(x, M) {
            if (!M) return this.gp(x);
            return this.gp(x, M);
        },
        _scp: function(x, M) {
            if (!M) return this.scp(x);
            return this.scp(x, M);
        },
        _versorInverse: function(M) {
            if (!M) return this.versorInverse();
            return this.versorInverse(M);
        },
        clone: function() {
            var M = new Multivector(slice(this.blades));
            M.bladesSorted = this.bladesSorted;
            return M;
        },
        toString: function(bvNames) {
            if (this.blades.length === 0) return "0";
            var result = '';
            this.blades.forEach(function(b, i) {
                var s = b.toString(bvNames);
                if (i === 0) {
                    result += s;
                } else if (s[0] === '-') {
                    result += " - ";
                    result += s.substr(1);
                } else {
                    result += " + ";
                    result += s;
                }
            });
            return result;
        },
        grade: function() {
            var g = -1;
            for (var i = 0; i < this.blades.length; i++) {
                var b = this.blades[i];
                if (g < 0) g = b.grade();
                else if (g !== b.grade()) return -1;
            }
            return (g < 0) ? 0 : g;
        },
        scalarPart: function() {
            return this.blades.reduce(function(sum, b) { return sum + (b.bitmap === 0 ? b.scale : 0); }, 0);
        },
        reverse: function() {
            return new Multivector(this.blades.map(function(b) { return b.reverse(); }));
        },
        gradeInversion: function() {
            return new Multivector(this.blades.map(function(b) { return b.gradeInversion(); }));
        },
        cliffordConjugate: function() {
            return new Multivector(this.blades.map(function(b) { return b.cliffordConjugate(); }));
        },
        ip: fun(
            [Multivector, Number, function(x, type) {
                var result = [];
                for (var i = 0; i < this.blades.length; i++) {
                    var b1 = this.blades[i];
                    for (var j = 0; j < x.blades.length; j++) {
                        var b2 = x.blades[j];
                        result.push(BasisBlade.ip(b1, b2, type));
                    }
                }
                return new Multivector(simplify(result));
            }],
            [Multivector, Metric, Number, function(x, M, type) {
                var result = [];
                for (var i = 0; i < this.blades.length; i++) {
                    var b1 = this.blades[i];
                    for (var j = 0; j < x.blades.length; j++) {
                        var b2 = x.blades[j];
                        result = result.concat(BasisBlade.ip(b1, b2, M, type));
                    }
                }
                return new Multivector(simplify(result));
            }]
        ),
        scalarProduct: fun(
            [Multivector, function(x) {
                return this.ip(x, InnerProductTypes.LEFT_CONTRACTION).scalarPart();
            }],
            [Multivector, Metric, function(x, M) {
                return this.ip(x, M, InnerProductTypes.LEFT_CONTRACTION).scalarPart();
            }]
        ),
        scp: function() {
            return this.scalarProduct.apply(this, slice(arguments));
        },
        norm_e: function() {
            var s = this.scp(this.reverse());
            if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'

            return Math.sqrt(s);
        },
        norm_e2: function() {
            var s = this.scp(this.reverse());
            if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'

            return s;
        },
        unit_e: function() {
            return this.unit_r();
        },
        unit_r: fun(
            [function() {
                var s = this.scp(this.reverse());
                if (s === 0.0) throw new Error("null multivector");
                return this.gp(1 / Math.sqrt(Math.abs(s)));
            }],
            [Metric, function(M) {
                var s = this.scp(this.reverse(), M);
                if (s === 0.0) throw new Error("null multivector");
                return this.gp(1 / Math.sqrt(Math.abs(s)));
            }]
        ),
        isNull: function(epsilon) {
            if (epsilon) {
                var s = this.norm_e2();
                return (s < epsilon * epsilon);
            }
            this.simplify();
            return this.blades.length === 0;
        },
        isScalar: function() {
            if (this.isNull()) return true;
            if (this.blades.length === 1) {
                return (this.blades[0].bitmap === 0);
            }
            return false;
        },
        dual: fun(
            [Number, function(dim) {
                var I = new Multivector(new BasisBlade((1 << dim) - 1, 1.0));
                return this.ip(I.versorInverse(), InnerProductTypes.LEFT_CONTRACTION);
            }],
            [Metric, function(M) {
                var I = new Multivector(new BasisBlade((1 << M.eigenMetric.length)-1, 1.0));
                return this.ip(I.versorInverse(), M, InnerProductTypes.LEFT_CONTRACTION);
            }]
        ),
        versorInverse: fun(
            //assumes `this` is a versor! no check is performed
            [function() {
                var R = this.reverse();
                var s = this.scp(R);
                if (s === 0.0) throw new Error("non-invertible multivector");
                return R.gp(1.0 / s);
            }],
            [Metric, function(M) {
                var R = this.reverse();
                var s = this.scp(R, M);
                if (s === 0.0) throw new Error("non-invertible multivector");
                return R.gp(1.0 / s);
            }]
        ),
        versorProduct: fun(
            //assumes V is a versor! no check is performed
            [Multivector, function(V) {
                return V.gp(this).gp(V.versorInverse());
            }],
            [Multivector, Metric, function(V, M) {
                return V.gp(this, M).gp(V.versorInverse(M), M);
            }]
        ),
        generalInverse: function(metric) {
            var dim = this.spaceDim();
            var i, M = [];
            for (i = 0; i < (1<<dim); i++) {
                M.push([]);
                for (var j = 0; j < (1<<dim); j++) {
                    M[i][j] = 0;
                }
            }
            var B = [];
            for (i = 0; i < (1 << dim); i++) {
                B[i] = BasisBlade.fromBitmap(i);
            }
            for (i = 0; i < this.blades.length; i++) {
                var b = this.blades[i];
                for (var j = 0; j < (1 << dim); j++) {
                    if (metric === undefined) {
                        addToMatrix(M, b, B[j], BasisBlade.gp(b, B[j]));
                    } else if (metric instanceof Metric) {
                        addToMatrix(M, b, B[j], BasisBlade.gp(b, B[j], metric));
                    } else if (Array.isArray(metric)) {
                        addToMatrix(M, b, B[j], BasisBlade.gp(b, B[j], metric));
                    }
                }
            }

            // try to invert matrix (can fail, then we throw an exception)
            var IM = null;
            try {
                IM = numeric.inv(M);
            } catch (e) {
                throw new Error("Multivector is not invertible");
            }

            // reconstruct multivector from first column of matrix
            var result = [];
            for (var j = 0; j < (1 << dim); j++) {
                var v = IM[j][0];
                if (v !== 0) {
                    B[j].scale = MathU.normalize(v);
                    result.push(B[j]);
                }
            }
            return new Multivector(result);

        },
        largestCoordinate: function() {
            this.simplify();
            return Math.max(0, Math.max.apply(Math, this.blades.map(getter('scale'))));
        },
        compress: function(epsilon) {
            epsilon = epsilon || 1e-13;
            this.simplify();

            // find maximum magnitude:
            var maxMag = 0.0;
            for (var i = 0; i < this.blades.length; i++) {
                maxMag = Math.max(Math.abs(this.blades[i].scale), maxMag);
            }
            if (maxMag === 0.0) {
                this.blades = [];
            } else {
                // premultiply maxMag
                maxMag = epsilon; // used to read *=
                var newBlades = [];
                this.blades.forEach(function(b) {
                    if (Math.abs(b.scale) >= maxMag) {
                        newBlades.push(b)
                    }
                });
                this.blades = newBlades;
            }
            return this;
        },
        topGradeIndex: function() {
            var maxG = 0;
            for (var i = 0; i < this.blades.length; i++) {
                maxG = Math.max(this.blades[i].grade(), maxG);
            }
            return maxG;
        },
        extractGrade: function(G) {
            if (Array.isArray(G)) {
                var maxG = Math.max(0, Math.max.apply(Math, G));

                // create boolean array of what grade to keep
                var keep = [];
                G.forEach(function(g) {
                    keep[g] = true;
                });

                // extract the grade, store in result:
                var result = [];
                this.blades.forEach(function(b) {
                    var g = b.grade();
                    if (g > maxG) return;
                    if (keep[g]) result.push(b.clone());
                });
                return new Multivector(result);
            } else {
                return this.extractGrade([G]);
            }
        },
        gradeUsage: function() {
            var gu = 0;
            this.blades.forEach(function(b) {
                gu |= 1 << b.grade();
            });
            return gu;
        },
        largestGradePart: function() {
            this.simplify();

            var maxGP = null;
            var maxNorm = -1.0;
            var gu = this.gradeUsage();
            for (var i = 0; i <= this.topGradeIndex(); i++) {
                if ((gu & (1 << i)) === 0) continue;
                var GP = this.extractGrade(i);
                var n = GP.norm_e();
                if (n > maxNorm) {
                    maxGP = GP;
                    maxNorm = n;
                }
            }
            return (maxGP == null) ? new Multivector() : maxGP;

        },
        spaceDim: function() {
            var maxD = 0;
            this.blades.forEach(function(b) {
                maxD = Math.max(Bits.highestOneBit(b.bitmap), maxD);
            });
            return maxD + 1;
        }
    };



    function bladesComparator(b1, b2) {
        if (b1.bitmap < b2.bitmap) return -1;
        else if (b1.bitmap > b2.bitmap) return 1;
        else return 0;
    }

    function simplify(list) {
        var prevBlade = null, i = 0;

        list.sort(bladesComparator);
        while (true) {
            if (i >= list.length) break;
            var curBlade = list[i];
            if (curBlade.scale === 0.0) {
                list.splice(i, 1);
                prevBlade = null;
            } else if (prevBlade !== null && (prevBlade.bitmap === curBlade.bitmap)) {
                prevBlade.scale += curBlade.scale;
                list.splice(i, 1);
            } else {
                prevBlade = curBlade;
                i++;
            }
        }
        for (i = 0; i < list.length; i++) {
            list[i].scale = MathU.normalize(list[i].scale);
        }

        //remove null blades
        return list.filter(function(b) { return b.scale !== 0; });
    }

    var addToMatrix = fun(
        [Array, BasisBlade, BasisBlade, BasisBlade, function(M, alpha, beta, gamma) {
            var v = M[gamma.bitmap][beta.bitmap];
            M[gamma.bitmap][beta.bitmap] = v + gamma.scale;
        }],
        [Array, BasisBlade, BasisBlade, Array, function(M, alpha, beta, gamma) {
            for (var i = 0; i < gamma.length; i++) {
                addToMatrix(M, alpha, beta, gamma[i]);
            }
        }]
    );


    Multivector.getBasisVector = function(idx) {
        return new Multivector(BasisBlade.fromBitmap(1 << idx));
    };
    exports.Multivector = Multivector;
})();