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
 * MeetJoin.java
 *
 * Created on October 12, 2005, 2:10 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */


/**
 * Computes the meet and join.
 * Usage:
 * new MeetJoin(a, b).getMeet();
 * Or:
 * MeetJoin MJ = new MeetJoin(a, b);
 * Multivector M = MJ.getMeet();
 * Multivector J = MJ.getJoin();
 *
 * @author  fontijne
 */
/**
 * <p>Both meet and join are computed simultaneously by this
 * algorithm by Ian Bell, which turns out to be somewhat
 * more efficient than computing the join and deriving the
 * meet from that.
 *
 * @return the meet and join of this with 'b'
 * The first array element returned is the meet (intersection),
 * the second array element is the join (union).
 */
(function() {
    var BasisBlade = require('./basis_blade').BasisBlade;
    var Multivector = require('./multivector').Multivector;
    var InnerProductTypes = require('./inner_product_types');

    function MeetJoin(a, b) {
        var la = a.largestCoordinate();
        var lb = b.largestCoordinate();

        var smallEpsilon = 10e-9;
        var largeEpsilon = 10e-4;

        // step one: check for near-zero input
        if ((la < smallEpsilon) || (lb < smallEpsilon)) {
            this.meet = new Multivector();
            this.join = new Multivector();
        }

        // check grade of input
        var ga = a.grade();
        var gb = b.grade();
        if (ga < 0) { // if we are not handed homogeneous multivectors, take the grade parts with the largest norm
            a = a.largestGradePart();
            ga = a.grade();
        }
        if (gb < 0) {
            b = b.largestGradePart();
            gb = b.grade();
        }
        // normalize (approximately) and swap (optionally)
        var ca, cb;
        if (ga <= gb) {
            // just normalize:
            ca = a.gp(1.0 / la);
            cb = b.gp(1.0 / lb);
        } else {
            // also swap:
            ca = b.gp(1.0 / lb);
            cb = a.gp(1.0 / la);

            var tempg = ga;
            ga = gb;
            gb = tempg;
        }
        // compute delta product & 'normalize'
        var d, _d = MeetJoin.deltaProduct(ca, cb);
        var gd = _d.grade();
        var ld = _d.largestCoordinate();
        d = _d.gp(1.0 / ld);

        // if delta product is scalar, we're done:
        if (gd === 0) {
            // meet = 1
            this.meet = ca;
            // join = computed from meet
            this.join = ca.op(ca.versorInverse().ip(cb, InnerProductTypes.LEFT_CONTRACTION));
            return;
        }

        // if grade of delta product is equal to ga + gb, we're done, too
        if (gd === ga + gb) {
            // a and b entirely disjoint
            // meet = 1
            this.meet = new Multivector(1.0);
            // join = computed from meet
            this.join = ca.op(cb);
            return;
        }
        // dimension of space we are working in:
        var dim = Math.max(ca.spaceDim(), cb.spaceDim());
        var I = new Multivector(new BasisBlade((1 << dim) - 1, 1.0));

        // init join to pseudoscalar
        var j = I;
        var Ej = dim - ((ga + gb + gd) >> 1); // compute excessity of join

        // check join excessity
        if (Ej === 0) {
            // join computed
            this.join = j;
            // meet = computed from join
            this.meet = cb.ip(j.versorInverse(), InnerProductTypes.LEFT_CONTRACTION).ip(ca, InnerProductTypes.LEFT_CONTRACTION);
            return;
        }

        // init meet
        var m = new Multivector(1.0);
        var Em = ((ga + gb - gd) >> 1); // compute excessity of meet

        // init s, the dual of the delta product:
        var s = d.ip(I.versorInverse(), InnerProductTypes.LEFT_CONTRACTION);

        // precompute inverse of ca
        var cai = ca.versorInverse();

        // todo: maybe we can improve: search only the largest basis blade of the not-delta product?
        for (var i = 0; i < dim; i++) {
            // compute next factor 'c'
            var c;

            // project 'tmpc' onto 's' (the dual of the delta product)
            // project using MHIP because 's' may just be a scalar some times?
            var tmpc = new Multivector(new BasisBlade(1 << i, 1.0)).ip(s, InnerProductTypes.MHIP);
            c = tmpc.ip(s, InnerProductTypes.MHIP); // no need to invert 's' here

            // todo: then this naughty step could be avoided:

            // check if 'c' is an OK candidate:
            if (c.largestCoordinate() < largeEpsilon)
                continue;

            // compute projection, rejection of 'c' wrt to 'ca'
            var cp, cr; // c projected, c rejected
            tmpc = c.ip(ca, InnerProductTypes.LEFT_CONTRACTION);
            cp = tmpc.ip(cai, InnerProductTypes.LEFT_CONTRACTION); // use correct inverse because otherwise cr != c - cp
            cr = c.subtract(cp);

            // if 'c' has enough of it in 'ca', then add to meet
            if (cp.largestCoordinate() > largeEpsilon) {
                m = m.op(cp);
                Em--; // reduce excessity of meet
                if (Em === 0) { // has the meet been computed???
                    this.meet = m;
                    // join = computed from meet
                    this.join = ca.op(ca.versorInverse().ip(cb, InnerProductTypes.LEFT_CONTRACTION));
                    return;
                }
            }

            if (cr.largestCoordinate() > largeEpsilon) {
                j = cr.ip(j, InnerProductTypes.LEFT_CONTRACTION);
                Ej--; // reduce excessity of join
                if (Ej === 0) { // has the join been computed???
                    this.join = j;
                    // meet = computed from join
                    this.meet = cb.ip(j.versorInverse(), InnerProductTypes.LEFT_CONTRACTION).ip(ca, InnerProductTypes.LEFT_CONTRACTION);
                    return;
                }
            }
        }
        throw new Error("meet & join algorithm failed!");

    }
    MeetJoin.deltaProduct = function(a, b) {
        var D = a.gp(b).compress();
        return D.extractGrade(D.topGradeIndex());
    };

    exports.MeetJoin = MeetJoin;
})();
