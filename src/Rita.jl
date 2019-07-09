module Rita

#import Unitful
using StaticArrays: SVector
using LinearAlgebra: dot, normalize, norm
#using CoordinateTransformations, StaticArrays, Unitful, UnitfulAngles

export refract_dir

# location
const Pnt = SVector{3, Float64}
# difference between locations
const Diff = SVector{3, Float64}
# direction
const Dir = SVector{3, Float64}
Dir(a, b, c)::Dir = normalize([a, b, c])
Dir(alpha, beta)::Dir = Dir([ cos(alpha)*sin(beta),
                                                sin(alpha)*sin(beta),
                                                cos(beta)])

# Inspired heavily by JuliaGeometry/RayTraceEllipsoid

struct Ray
    orig::Pnt
    dir::Dir
    Ray(o::Pnt, d::Dir) = new(o, d)
end

# ray along z from origin
Ray() = Ray(Vec(0,0,0), Vec(1,0,0))

# dot(in, n) should be >0 in most cases
mirror_dir(in::Dir, n::Dir)::Dir = in - 2*dot(in, n)*n
# Snellius law
# Does also handle total internal reflection!
function refract_dir(in::Dir, n::Dir, n1::Float64, n2::Float64)::Dir
    in_parallel = dot(in, n) * n
    in_orthogonal = in - in_parallel
    q = n1/n2
    # println(in, " ", n, " ", n1, " ", n2)
    out = try
        q * in_orthogonal + sqrt(1 - q^2 * dot(in_orthogonal, in_orthogonal)) * in_parallel
    catch DomainError
        # total internal reflection
        mirror_dir(in, n)
    end
    normalize(out)
end

# Code excerpt from RayTraceEllipsoid.jl
function intersect_unit_sph(ray::Ray)::Tuple{Float64, Float64}
    v = - ray.orig
    b = dot(v, ray.dir)
    disc = b*b - dot(v, v) + 1
    if disc >= 0
        d = sqrt(disc)
        t2 = b + d
        if t2 >= 0
            t1 = b - d
            return t1 > 0 ? (t1, t2) : (Inf, t2)
        end
    end
    return (Inf, Inf)
end

end # module
