using Multroot
using Base.Filesystem
using Polynomials

# assuming pwd = package/test
const ts = "testsuit"
println(pwd())

function runtest(fun::Function)
    println("testing multroot with ", fun)
    p, z = fun()
    p = p isa Poly ? reverse(p.a) : p
    zz, bkerr, job = multroot(p)
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

for fun in functions
    runtest(fun)
end    
