using PyPlot
includet("/home/547/ar0754/.julia/dev/HiQGA/zz_portalcurtains/RDP.jl")
##
XY = RDP.collectpoints()
x, y = RDP.smoothline(XY, λ²=10000, finefactor=1000, regtype=:R2, fname="base_menindee.txt")
##
XY1 = RDP.collectpoints()
x, y = RDP.smoothline(XY1, λ²=10000, finefactor=1000, regtype=:R2, fname="base_1_menindee.txt")
##
XY2a = RDP.collectpoints()
x, y = RDP.smoothline(XY2a, λ²=10000, finefactor=1000, regtype=:R2, fname="base_2a_menindee.txt")
##
XY2b = RDP.collectpoints()
x, y = RDP.smoothline(XY2b, λ²=10000, finefactor=1000, regtype=:R2, fname="base_2b_menindee.txt")
## 
XY3a = RDP.collectpoints()
x, y = RDP.smoothline(XY1a, λ²=10000, finefactor=1000, regtype=:R2, fname="base_3a_menindee.txt")
## 
XY3b = RDP.collectpoints()
x, y = RDP.smoothline(XY1b, λ²=10000, finefactor=1000, regtype=:R2, fname="base_3b_menindee.txt")
## read them all
x, y = [[xy[i] for xy in RDP.readpoints(["base_menindee.txt", "base_1_menindee.txt", "base_2a_menindee.txt", "base_2b_menindee.txt", 
                    "base_3a_menindee.txt", "base_3b_menindee.txt"])] for i in 1:2]
map(axp) do ax
    for a in ax[1:8]
        for (xx, yy) in zip(x, y)
            a.plot(xx, yy, "--k", linewidth=0.5)
        end
    end
end