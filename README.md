Geometric algebra implementation in JavaScript.  Based on the Java reference implementation from [Geometric Algebra for Computer Science](http://www.amazon.com/Geometric-Algebra-Computer-Science-Revised/dp/0123749425). See also http://www.geometricalgebra.net/

Uses eigenvalue calculator code from http://www.akiti.ca/EigR12Solver.html

Emphasis is on correctness and legibility, not speed or efficiency.

To run tests (requires jasmine-node):
    jasmine-node spec

To compile the parser (requires [PEG.js](http://pegjs.majda.cz/)):
    pegjs --cache parser.pegjs

TODO
===============================
* compiler/interpreter for symbolic expressions
* render to WebGL
* multivector logarithm
* full suite of conformal model constructs
