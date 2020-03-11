function large03()
#
#  test polynomial suggested by Z Zeng
#
    p, z = large02()
    p = conv(p, p)
    z.mult .*= 2
    p, z
end
