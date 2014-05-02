var fun = require('../util/pattern_match').fun;

describe('Pattern Matcher', function() {
    it('should work for regular functions', function() {
        var f = fun(
            [Array, function(a) { return a[0]; }],
            [Number, function(n) { return n * 2; }],
            [Number, Number, function(a, b) { return a + b; }],
            [Number, String, function(n, s) { return "The string is: " + s;}],
            [Boolean, function(b) { return !b; }],
            [String, function(s) { return s === 'ok'; }],
            [fun.Any, function(x) { return 42; }]
        );

        expect(f([42])).toBe(42);
        expect(f([1,2,3])).toBe(1);
        expect(f(5, 12)).toBe(17);
        expect(f(42, "hello, world")).toBe("The string is: hello, world");
        expect(f(42)).toBe(84);
        expect(f(true)).toBe(false);
        expect(f('ok')).toBe(true);
        expect(f({})).toBe(42);
    });

    it('should work for "instance" functions', function() {
        function Foo(bar) {
            this.bar = bar;
        }
        Foo.prototype = {
            getBar: function() {
                return this.bar;
            },
            doStuff: fun(
                [Number, function(x) { return this.getBar() * x; }],
                [Array, function(a) { return a.length === this.getBar(); }],
                [Foo, function(foo) { return "foo!"; }]
            )
        };
        var bar = 2;
        function Baz() {
        }
        Baz.prototype = new Foo(bar);
        var foo = new Foo(bar);
        var baz = new Baz(bar);
        expect(foo.getBar()).toBe(bar);
        expect(foo.doStuff(42)).toBe(bar * 42);
        expect(foo.doStuff(new Array(bar))).toBe(true);
        expect(foo.doStuff([1])).toBe(false);
        expect(baz.doStuff(2)).toBe(4);
        expect(baz.doStuff(baz)).toBe("foo!");
    });
});

