"""
  module Multroot

Calculate mutiple roots of a real or complex polymomial exploiting the
multiplicity structure of the root.
Using methods develped by Zhang et. al.
"""
module Multroot

using Printf
using LinearAlgebra

export multroot, gcdroot, pejroot
export sylves, sylves1, sylmat, cauchymt, scalsq
export hqrt, hessqr, forsub, backsub
export zminsv, zminsv1
export polyzeros, polymult, polytransform

include("polyzeros.jl")

include("sylves.jl")
include("sylmat.jl")
include("polymult.jl")
include("zminsv.jl")
include("hessqr.jl")
include("hqrt.jl")
include("cauchymt.jl")
include("scalsq.jl")
include("forsub.jl")
include("backsub.jl")
include("gcdgn.jl")
include("multroot.jl")
include("pejroot.jl")
include("gcdroot.jl")

include("polyscale.jl")
include("horner.jl")

end # module


