(function() {
    var calcEigr = require('./eigen/real_eig').calcEigr;
    var fun = require('./util/pattern_match').fun;
    function cloneMatrix(a) {
        var output = [];
        a.forEach(function(row) {
            output.push(row.slice(0));
        });
        return output;
    }
    function transpose(m) {
        var output = [];
        for (var i = 0; i < m.length; i++) {
            output.push([]);
        }
        for (var i = 0; i < m.length; i++) {
            for (var j = 0; j < m[i].length; j++) {
                output[i][j] = m[j][i];
            }
        }
        return output;
    }
    function isDiagonal(m) {
        for (var i = 0; i < m.length; i++) {
            for (var j = 0; j < m[i].length; j++) {
                if (i !== j && m[i][j] !== 0.0) return false;
            }
        }
        return true;
    }
    function matrixEqual(m1, m2) {
        if (m1.length !== m2.length) return false;
        for (var i = 0; i < m1.length; i++) {
            if (m1[i].length !== m2[i].length) return false;
            for (var j = 0; j < m1[i].length; j++) {
                if (m1[i][j] !== m2[i][j]) return false;
            }
        }
        return true;
    }
    function isSymmetric(m) {
        return matrixEqual(m, transpose(m));
    }

    function Metric(m) {
        this.matrix = cloneMatrix(m);
        if (!isSymmetric(this.matrix)) {
            throw new Error("The metric matrix must be symmetric");
        }

        // compute eigen value decomposition
        this.eig = calcEigr(cloneMatrix(this.matrix));

        this.invEigMatrix = transpose(this.eig.matrix);

        this.eigenMetric = this.eig.wr.slice(0);

        this.isDiagonal = isDiagonal(this.matrix);
         if (!this.isDiagonal) {
             this.isEuclidean = this.isAntiEuclidean = false;
         } else {
             this.isEuclidean = this.isAntiEuclidean = true;
             for (var i = 0; i < this.matrix.length; i++) {
                if (this.matrix[i][i] !== 1.0) {
                    this.isEuclidean = false;
                }
                if (this.matrix[i][i] !== -1.0) {
                    this.isAntiEuclidean = false;
                }
            }
         }
    }
    exports.Metric = Metric;
})();
