function large04()
#
#  test polynomial suggested by Z Zeng
#
    p, z = large03()
    p = conv(p, p)
    z.mult .*= 2
    p, z
end
