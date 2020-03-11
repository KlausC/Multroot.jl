function large05()
#
#  test polynomial suggested by Z Zeng
#
    p, z = large04()
    p = conv(p, p)
    z.mult .*= 2
    p, z
end
