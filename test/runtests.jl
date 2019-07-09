using Rita
using Test
using QuickCheck
using LinearAlgebra: norm

import QuickCheck.generator
generator(::Type{Rita.Dir}, size) = Rita.Dir(rand(0:2*pi), rand(-pi:pi))

# Some simple 2d snellius in order to test some simple cases
function snellius_2d(in_angle::Float64, n1::Float64, n2::Float64)::Float64
    asin(sin(in_angle)*n1/n2)
end

@test Rita.refract_dir(Rita.Dir([0,-1, 0]), Rita.Dir([0, 0, -1]), 1., 1.) ≈ Rita.Dir([0,-1, 0])
@test Rita.refract_dir(Rita.Dir([0, 0, 1]), Rita.Dir([0, 0, -1]), 1., 1.) ≈ Rita.Dir([0, 0, 1])
# no total reflection 
condproperty((dir::Rita.Dir, n::Rita.Dir, n1::Float64, n2::Float64) -> norm(Rita.refract_dir(dir, n, n1, n2)) ≈ 1.,
             100, 1000, (dir::Rita.Dir, n::Rita.Dir, n1::Float64, n2::Float64) -> 1.0 <= n1 <= n2 <= 5.0)
