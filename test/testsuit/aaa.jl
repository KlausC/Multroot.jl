import Multroot: convolute

conv(a, b) = convolute(a, b)
conv(a::Poly, b::Poly) = conv(p.a, b.a)
