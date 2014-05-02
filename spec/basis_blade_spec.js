var BasisBlade = require('../basis_blade').BasisBlade;
var Multivector = require('../multivector').Multivector;
var InnerProductTypes = require('../inner_product_types');
var MeetJoin = require('../meet_join').MeetJoin;
var Metric = require('../metric').Metric;

var e1 = Multivector.getBasisVector(0);
var e2 = Multivector.getBasisVector(1);
var e3 = Multivector.getBasisVector(2);
var e4 = Multivector.getBasisVector(3);

beforeEach(function() {
  this.addMatchers({
    toBeEquivalent: function(expected) {
      return this.actual.equals(expected);
    }
  });
});

describe('Clone', function() {
    it('should clone', function() {
        expect(e1.clone()).toBeEquivalent(e1);
        expect(e1.clone()).not.toBe(e1);
    });
});

describe('Equality', function() {
    it('should have basis vectors equal to themselves', function() {
        expect(e1).toBeEquivalent(e1);
        expect(e1).toBeEquivalent(Multivector.getBasisVector(0));
    });
    it('should have basis vectors not equal to other basis vectors', function() {
        expect(e1).not.toBeEquivalent(e2);
    });
});

describe('Addition', function() {
    it('should work when adding basis vectors', function() {
        expect(e1.add(e2)).toBeEquivalent(e2.add(e1));
        expect(e1.add(e1)).toBeEquivalent(e1.gp(2));
    });
    it('should work when adding scalars', function() {
        expect(e1.add(42).scalarPart()).toBe(42);
    });
});
describe('Subtraction', function() {
    it('should work when subtracting basis vectors', function() {
        expect(e1.subtract(e1)).toBeEquivalent(0);
    });
    it('should work when subtracting scalars', function() {
        expect(e1.subtract(42).scalarPart()).toBe(-42);
    });
});

describe('Outer product', function() {
    describe('axioms: ', function() {
        it('should be antisymmetric', function() {
            expect(e1.op(e2)).toBeEquivalent(e2.op(e1).gp(-1));

            var b1 = e1.add(e2);
            expect(b1.op(e1)).toBeEquivalent(e1.op(b1).gp(-1));
        });
        it('should scale', function() {
            var scalar = 2;
            var left = e1.op(e2.gp(scalar)); // e1 ^ (βe2)
            var right = e1.op(e2).gp(scalar); // β(e1 ^ e2)
            expect(left).toBeEquivalent(right);
        });
        it('should distribute', function() {
            var left = e1.op(e2.add(e3)); //e1 ^ (e2 + e3)
            var right = e1.op(e2).add(e1.op(e3)); // e1^e2 + e1^e3
            expect(left).toBeEquivalent(right);
        });
    });

    it('should be zero for a basis blade with itself', function() {
        expect(e1.op(e1)).toBeEquivalent(0);
    });
    it('should raise the grade properly', function() {
        expect(e1.op(e2).grade()).toBe(e1.grade() + e2.grade());
        expect(e1.op(e2).op(e3).grade()).toBe(e1.grade() + e2.grade() + e3.grade());
    });
    it('should stringify correctly', function() {
        expect(e1.op(e2).toString()).toBe('1*e1^e2');
        expect(e1.op(e2).op(e3).toString()).toBe('1*e1^e2^e3');
    });

    it('should be associative for basis blades', function() {
        expect(e1.op(e2).op(e3)).toBeEquivalent(e1.op(e2.op(e3)));
    });
    it('should be reversible', function() {
        expect(e1.op(e2).reverse()).toBeEquivalent(e2.op(e1));
        expect(e1.op(e2).op(e3).reverse()).toBeEquivalent(e3.op(e2).op(e1));
    });
});

describe('Scalar product', function() {
    describe('axioms: ', function() {
        it('should be symmetric', function() {
            var v1 = e1.add(e2);
            var v2 = e1.add(e3);
            expect(v1.scp(v2)).toBe(v2.scp(v1));

            expect(v1.reverse().scp(v2.reverse())).toBe(v1.scp(v2));
            expect(v1.reverse().scp(v2.reverse())).toBe(v2.scp(v1));
        });

        it('should compute the angle between subspaces', function() {
            expect(e1.scp(e2) / (e1.norm_e() * e2.norm_e())).toBe(0);

            var v1 = e1.add(e2);
            expect(v1.scp(e2.reverse()) / (v1.norm_e() * e2.norm_e())).toBeCloseTo(Math.cos(Math.PI / 4), 15);

            var v2 = e2.add(e3);
            expect(v1.scp(v2.reverse()) / (v1.norm_e() * v2.norm_e())).toBeCloseTo(0.5, 15);
        });
    });

    describe('in a Euclidean metric', function() {
        it('basis vector times itself', function() {
            expect(e1.scp(e1)).toBe(1);
        });
        it('basis vector times a different basis vector', function() {
            expect(e1.scp(e2)).toBe(0);
        });
    });
});

describe('Geometric product', function() {
    it('should work with scalars', function() {
        expect(e1.gp(2)).toBeEquivalent(e1.add(e1));
        expect(e1.gp(0)).toBeEquivalent(0);
    });
    it('2-blades should square to -1 in euclidean metric', function() {
        var I = e1.op(e2).op(e3);
        expect(I.gp(I)).toBeEquivalent(-1);
    });
});

describe('Contraction', function() {
    describe('axioms: ', function() {
        it('should match inner product', function() {
            var scale = 2;
            var scalar = new Multivector(new BasisBlade(0, scale));
            var v1 = e1.add(e2), v2 = e2.add(e3);
            expect(scalar.ip(e1, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(e1.gp(scale));

            expect(scalar.ip(v1, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(v1.gp(scale));

            expect(scalar.ip(scalar, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(scale * scale);

            expect(v1.ip(v2, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(v1.scp(v2));
        });

        it('should contract to zero against a scalar', function() {
            var scale = 2;
            var scalar = new Multivector(new BasisBlade(0, scale));
            var v1 = e1.add(e2);

            expect(e1.ip(scalar, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(0);
            expect(v1.ip(scalar, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(0);
        });

        it('should be separable', function() {
            var p1 = e1;
            var p2 = e2;
            var p3 = e3;
            expect(p1.op(p2).ip(p3, InnerProductTypes.LEFT_CONTRACTION)).toBeEquivalent(p1.ip(p2.ip(p3, InnerProductTypes.LEFT_CONTRACTION), InnerProductTypes.LEFT_CONTRACTION));
        });

        it('should be distributive', function() {
            var p1 = e1.op(e2);
            var p2 = e1.op(e2.gp(3));
            var p3 = e1.op(e2.add(e3));

            var res = (p1.add(p2).ip(p3, InnerProductTypes.LEFT_CONTRACTION));
            expect(res).toBeEquivalent(p1.ip(p3, InnerProductTypes.LEFT_CONTRACTION).add(p2.ip(p3, InnerProductTypes.LEFT_CONTRACTION)));
        });

        it('should work', function() {
            var p1 = e1.op(e2.gp(2));
            var p2 = e1.op(e2).op(e3);
            var res = p1.ip(p2, InnerProductTypes.LEFT_CONTRACTION);
            expect(res).toBeEquivalent(e3.gp(-2));
        });

        it('implicit definition', function() {
            // (X ^ A) % B == X % (A << B)
            var A = e1.op(e2.gp(2));
            var X = e3.op(e4).add(e2.op(e3));
            var B = e1.op(e2).op(e3).op(e4).gp(2);

            //left contraction
            expect(X.op(A).scp(B)).toBe(X.scp(A.ip(B, InnerProductTypes.LEFT_CONTRACTION)));

            //right contraction
            expect(B.scp(A.op(X))).toBe(B.ip(A, InnerProductTypes.RIGHT_CONTRACTION).scp(X));

        });
    });
});

describe('Dual', function() {
    describe('in 2d', function() {
        it('should match the definition', function() {
            var dual = e1.dual(2);
            expect(dual).toBeEquivalent(e2.gp(-1));
            expect(dual.dual(2)).toBeEquivalent(e1.gp(-1));
        });
    });


    describe('in 3d', function() {
        it('should match the definition', function() {
            var dual = e1.dual(3);
            expect(dual).toBeEquivalent(e2.op(e3).gp(-1));
            expect(dual.dual(3)).toBeEquivalent(e1.gp(-1));
        });
    });

    describe('axioms', function() {
        it('should work in 3d', function() {
            expect(e1.op(e2).dual(3)).toBeEquivalent(e1.ip(e2.dual(3), InnerProductTypes.LEFT_CONTRACTION));
        });
    });
});


describe('Meet and Join', function() {
    it('should compute', function() {
        var mj = new MeetJoin(e1.op(e2), e3.gp(2).add(e2.gp(1).add(e1.gp(3))));
    });
});
