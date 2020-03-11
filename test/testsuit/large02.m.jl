function large02()
#
#  test polynomial suggested by Z Zeng
#

    p, z = large01()
    p = conv(p, p)
    z.mult .*= 2
    p, z
end
