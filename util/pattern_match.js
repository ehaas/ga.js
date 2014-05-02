(function() {
    var slice = Array.prototype.slice, toString = Object.prototype.toString;
    function isNumber(n) {
        return typeof n === 'number' && isFinite(n);
    }
    function isFunction(f) {
        return toString.call(f) === '[object Function]';
    }
    function isString(s) {
        return toString.call(s) === '[object String]';
    }
    function isMatch(type, arg) {
        switch(type) {
            case Number: return isNumber(arg);
            case Array: return Array.isArray(arg);
            case Object: return toString.call(arg) === '[object Object]';
            case Boolean: return arg === true || arg === false;
            case String: return isString(arg);
            case Function: return isFunction(arg);
            case Matcher.Any: return true;
            default:
                if (isFunction(type)) {
                    return type.prototype && type.prototype.isPrototypeOf && type.prototype.isPrototypeOf(arg);
                } else if(isString(type)) {
                    return type === arg;
                }
        }
    }
    function Matcher() {
        var defs = slice.call(arguments, 0);
        return function() {
            var args = slice.call(arguments, 0), i;
            for (i = 0; i < defs.length; i++) {
                var def = defs[i], j, allMatched = true;
                for (j = 0; j < def.length - 1; j++) {
                    if (!isMatch(def[j], args[j])) {
                        allMatched = false;
                        break;
                    }
                }
                if (allMatched && def.length - 1 === args.length) {
                    var fn = def[def.length - 1];
                    return fn.apply(this, args);
                }
            }
            throw new Error("No matching def!");
        }
    }
    Matcher.Any = function() { };


    exports.fun = Matcher;
})();
