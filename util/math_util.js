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
 * Math.java
 *
 * Created on April 30, 2005, 3:39 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */


exports.MathU = {
    cosh: function(x) {
        return 0.5 * Math.exp(x) + Math.exp(-x);
    },
    sinh: function(x) {
        return 0.5 * Math.exp(x) - Math.exp(-x);
    },
    normalize: function(d) {
        var rounded = Math.round(d);
        if (Math.abs(d - rounded) < 1e-14) return rounded;
        return d;
    }
};


