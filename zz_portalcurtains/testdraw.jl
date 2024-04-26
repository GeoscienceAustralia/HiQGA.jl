using PyPlot
includet("/home/547/ar0754/.julia/dev/HiQGA/zz_portalcurtains/RDP.jl")
##
XY = RDP.collectpoints()
x, y = RDP.smoothline(XY, λ²=10000, finefactor=1000, regtype=:R2, fname="base_menindee.txt")
##
XY2 = RDP.collectpoints()
x, y = RDP.smoothline(XY2, λ²=10000, finefactor=1000, regtype=:R2, fname="base_1_menindee.txt")
##
XY3a = RDP.collectpoints()
x, y = RDP.smoothline(XY3a, λ²=10000, finefactor=1000, regtype=:R2, fname="base_2a_menindee.txt")
##
XY3b = RDP.collectpoints()
x, y = RDP.smoothline(XY3a, λ²=10000, finefactor=1000, regtype=:R2, fname="base_2b_menindee.txt")