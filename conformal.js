(function() {
    var Multivector = require('./multivector').Multivector;
    var Metric = require('./metric').Metric;
    var fun = require('./util/pattern_match').fun;

    var Conformal = {};

    var metric = new Metric([
        [0, 0, 0, 0, -1],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [-1, 0, 0, 0, 0]
    ]);
    Conformal.metric = metric;

    var no = Multivector.getBasisVector(0);
    var e1 = Multivector.getBasisVector(1);
    var e2 = Multivector.getBasisVector(2);
    var e3 = Multivector.getBasisVector(3);
    var ni = Multivector.getBasisVector(4);

    Conformal.point = fun(
        [Number, Number, Number, function(x, y, z) {
            var p = e1.gp(x).add(e2.gp(y)).add(e3.gp(z));
            return no.add(p).add(p.gp(p).gp(0.5).gp(ni));
        }],
        [Multivector, function(v) {
            // var p = e1.gp(x).add(e2.gp(y)).add(e3.gp(z));
            // return no.add(p).add(p.gp(p).gp(0.5).gp(ni));
        }]
    );

    Conformal.bv = {
        no: no,
        e1: e1,
        e2: e2,
        e3: e3,
        ni: ni
    };
    Conformal.bvNames = ['no', 'e1', 'e2', 'e3', 'ni'];


    Conformal.dualPlane = function(x, y, z) {
        var n = e1.gp(x).add(e2.gp(y)).add(e3.gp(z));
        return n.add(ni);
    };
    Conformal.dualSphere = function(p, r) {
        return p.subtract(ni.gp(r * r / 2));
    }
    Conformal.dualImaginarySphere = function(p, r) {
        return p.add(ni.gp(r*r/2));
    }

    Conformal.translationVersor = function(x, y, z) {
        var n = e1.gp(x).add(e2.gp(y)).add(e3.gp(z));
        return n.gp(ni).gp(-1/2).add(1);
    };

    Conformal.classify = function(x) {
        var square = x.gp(x, metric);
        var ninf = ni.gp(-1);
        var probe = ninf.scp(x, metric);
        if (square.equals(0) && probe !== 0) {
            return {type: 'point'};
        } else if (!square.equals(0) && probe === 0) {
            return {type: 'dualPlane'};
        } else if (square.scalarPart() > 0 && probe !== 0) {
            return {type: 'dualRealSphere'};
        } else if (square.scalarPart() < 0 && probe !== 0) {
            return {type: 'dualImaginarySphere'};
        }
    };
    exports.Conformal = Conformal;

})();
