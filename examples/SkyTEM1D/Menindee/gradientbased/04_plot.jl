dr, dz = 50, 2 # interpolation in m
lnames = [200613] # array of lines wanted in survey to plot
idx = [[1, 5]] # array of index arrays with indices per line in lnames 
transD_GP.plotconvandlast(soundings, dr, dz; zall, 
    figsize=(15,7), vmin=lo, vmax=hi, postfix="_β²_$(β²)_$(regtype)_bg_$(round(10. ^σ0, sigdigits=4))Spm", 
    prefix=zipsaveprefix, lnames, idx, plotforward=true, aem_in=aem,
    preferEright=true, showplot=true, logscale=true)
@info "done"
