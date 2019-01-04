"""
Misplaced and miscellaneous functions
"""
module Utils

using SpecialFunctions: gamma
using StatsFuns: gammacdf
using HypergeometricFunctions: _₂F₁

"""
Wrapper around hypergeometric function 2F1
"""
function hyp2f1(a, b, c, z)
    return _₂F₁(a, b, c, z)
end

"""
Wrapper around regularized lower incomplete gamma function from Rmath
"""
function gamma_inc(a, x)
    return gammacdf(a, 1, x)
end


"""
Redefined incomplete beta function in terms of hypergeometric function 2F1.
"""
function B_inc(a, b, z)
    F(z) = hyp2f1(a, 1.0 - b, a + 1.0, z)
    return @. a ^ -1.0  * z ^ a * F(z)
end

end
