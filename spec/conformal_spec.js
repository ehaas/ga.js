var Conformal = require('../conformal').Conformal;
var Multivector = require('../multivector').Multivector;
var InnerProductTypes = require('../inner_product_types');
var MeetJoin = require('../meet_join').MeetJoin;
var Metric = require('../metric').Metric;

describe('Conformal model', function() {
    var metric = Conformal.metric;
    var no = Conformal.bv.no;
    var e1 = Conformal.bv.e1;
    var e2 = Conformal.bv.e2;
    var e3 = Conformal.bv.e3;
    var ni = Conformal.bv.ni;
    var bvNames = Conformal.bvNames;
    it('should calculate scalar product correctly', function() {
        expect(no.scp(no, metric)).toBe(0);
        expect(no.scp(e2, metric)).toBe(0);
        expect(no.scp(ni, metric)).toBe(-1);
    });

    it('should make conformal points', function() {
        var p = Conformal.point(1, 2, 3);
        expect(p.gp(p, metric)).toBeEquivalent(0);
        expect(ni.scp(p, metric)).toBe(-1);
        var scaled_point = p.gp(3);
        expect(ni.scp(scaled_point, metric)).toBe(-3);
    });
    it('should make conformal planes', function() {
        var plane = Conformal.dualPlane(1, 2, 3);
        expect(plane.gp(plane, metric)).not.toBe(0);
        expect(ni.gp(-1).scp(plane, metric)).toBe(0);
    });
    it('should make dual spheres', function() {
        var radius = 4;
        var p = Conformal.point(1, 2, 3);
        var sphere = Conformal.dualSphere(p, radius);

        expect(sphere.gp(sphere, metric)).toBeEquivalent(radius * radius);
        expect(ni.gp(-1).scp(sphere, metric)).toBe(1);

        var iSphere = Conformal.dualImaginarySphere(p, radius);
        expect(iSphere.gp(iSphere, metric)).toBeEquivalent(-radius * radius);
        expect(ni.gp(-1).scp(iSphere, metric)).toBe(1);
    });
    it('should invert versors properly', function() {
        var T = Conformal.translationVersor(2, 4, 6);
        var T_1 = T.versorInverse(metric);
        expect(T.generalInverse(metric)).toBeEquivalent(T_1);
        expect(T_1.generalInverse(metric)).toBeEquivalent(T);
    });
    it('should handle pseudoscalar properly', function() {
        var I = no.op(e1).op(e2).op(e3).op(ni);
        var I_1 = no.op(e1.op(e2).op(e3).reverse()).op(ni);

        expect(I.generalInverse(metric)).toBeEquivalent(I_1);
        expect(I.gp(I_1, metric)).toBeEquivalent(1);
    });

    describe('translations', function() {
        it('should translate the origin', function() {
            var p = Conformal.point(2, 4, 6);

            var T = Conformal.translationVersor(2, 4, 6);
            var T_1 = T.versorInverse(metric);
            expect(T.gp(T_1, metric)).toBeEquivalent(1);
            expect(T.gp(no, metric).gp(T_1, metric)).toBeEquivalent(p);
            expect(no.versorProduct(T, metric)).toBeEquivalent(p);
        });
        it('should correctly compute the exponential', function() {
            var T = Conformal.translationVersor(1, 2, 3);

            var t = e1.add(e2.gp(2)).add(e3.gp(3));
            var exponent = t.gp(ni, metric).gp(-1/2);
            expect(exponent.exp(metric)).toBeEquivalent(T);
        });
    });

    describe('rotations', function() {
        it('should do simple rotations correctly', function() {
            var p = Conformal.point(1, 2, 3);
            var p_r = Conformal.point(-1, -2, 3);
            var R = e1.gp(e2, metric);

            var rotatedPoint = R.gp(p, metric).gp(R.versorInverse(metric), metric);
            expect(p_r).toBeEquivalent(rotatedPoint);
        });
    });

    describe('versor composition', function() {
        it('should compose properly', function() {
            var T = Conformal.translationVersor(1, 2, 3);

            var n1 = e1.add(e2.gp(3));
            var n2 = e1.add(e2).add(e3.gp(2));
            var R = n1.gp(n2, metric);
            var rotate_then_move = R.versorProduct(T, metric);

            var n1_moved = n1.versorProduct(T, metric);
            var n2_moved = n2.versorProduct(T, metric);
            var move_then_rotate = n1_moved.gp(n2_moved, metric);

            expect(move_then_rotate).toBeEquivalent(rotate_then_move);
        });
    });

    describe('Classify', function() {
        it('should classify points', function() {
            var p = Conformal.point(1, 2, 3);
            var classification = Conformal.classify(p);
            expect(classification.type).toEqual('point');
        });
        it('should classify dual planes', function() {
            var p = Conformal.dualPlane(1, 2, 3);
            var classification = Conformal.classify(p);
            expect(classification.type).toEqual('dualPlane');
        });
        it('should classify dual spheres', function() {
            var p = Conformal.point(1, 2, 3);
            var sphere = Conformal.dualSphere(p, 2);
            var classification = Conformal.classify(sphere);
            expect(classification.type).toEqual('dualRealSphere');
        });
        it('should classify dual imaginary spheres', function() {
            var p = Conformal.point(1, 2, 3);
            var imaginarySphere = Conformal.dualImaginarySphere(p, 2);
            var classification = Conformal.classify(imaginarySphere);
            expect(classification.type).toEqual('dualImaginarySphere');
        });
    });


});
