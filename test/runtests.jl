using Multroot
using Base.Filesystem
using Polynomials
using LinearAlgebra

import Multroot: diffmult, diffzeros

using Test

# assuming pwd = package/test
const ts = "testsuit"

function runtest(fun::Function)
    println("testing multroot with ", fun)
    @testset "$fun" begin
        p, z, job_exp, bkerr_exp, mult_exp, znorm_exp = wrap_fun(fun)
        zz, bkerr, pjcond, job = multroot(p)
        z = sort(z)
        zz = sort(z)
        job == 0 && return
        @test job == job_exp
        @test bkerr <= bkerr_exp
        @test length(zz.mult) == length(z.mult)
        lc = length(zz.mult) == length(z.mult)
        @test lc && diffmult(zz, z) <= mult_exp
        @test lc && diffzeros(zz, z) <= znorm_exp
    end
end

function wrap_fun(fun::Function, args...)
    p, z, = fun(args...)
    pp  = p isa Poly ? reverse(p.a) : p
    job_exp = 1
    bkerr_exp = Inf
    mult_exp = 0
    znorm_exp = norm(z.z, Inf) * 1e-8
    pp, z, job_exp, bkerr_exp, mult_exp, znorm_exp
end

functpath, io = mktemp(cleanup = true)
println(io, "const functions = Function[]")
for fname in filter(x -> endswith(x, "jl"), readdir(ts))
    fstring = read(abspath(ts, fname), String)
    if endswith(fname, ".m.jl")
        println(io, "fun = ")
        write(io, fstring)
        println(io)
        println(io, "push!(functions, fun)")
    else
        println(io)
        write(io, fstring)
        println(io)
    end
end
close(io)
include(functpath)

@testset "testsuite" begin
    for fun in functions
        runtest(fun)
    end
end

