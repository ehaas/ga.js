{
	var scope = {};
}
start
  = additive

additive
  = left:geometric ws "+" ws right:additive { return ['+', left, right]; }
  / geometric

geometric "Geometric product"
  = left:outer ws "*" ws right:geometric { return ['*', left, right]; }
  / outer

outer "Outer product"
  = left:scalar ws "^" ws right:outer { return ['^', left, right]; }
  / scalar

scalar
  = left:lcontract ws "%" ws right:scalar { return ['%', left, right]; }
  / lcontract

lcontract "Left contraction"
  = left:rcontract ws "<<" ws right:lcontract { return ['<<', left, right]; }
  / rcontract

rcontract "Right contraction"
  = left:meet ws ">>" ws right:rcontract { return ['>>', left, right]; }
  / meet

meet
  = left:join ws "&" ws right:meet { return ['&', left, right]; }
  / join

join
  = left:primary ws "|" ws right:join { return ['|', left, right]; }
  / primary

primary
  = integer
  / ref
  / "(" additive:additive ")" { return additive; }

ref "Reference"
  = id:(refstart refchar*) { return ['ref', Array.prototype.concat.apply([], id).join('')]; }

refstart
  = letter:[a-z]i { return letter; }

refchar
  = char:[a-z0-9_]i { return char; }

integer "integer"
  = digits:[0-9]+ { return parseInt(digits.join(''), 10); }

ws "Whitespace"
  = [\t ]*