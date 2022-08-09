dxplot, dzplot = 5, 2
transD_GP.SkyTEM1DInversion.splitlineconvandlast(soundings, dxplot, dzplot; zstart=zstart, extendfrac=extendfrac, dz=dz, nlayers=nlayers, 
    figsize=(8.5,5), vmin=lo, vmax=hi, postfix="_β²_$(β²)_$(regtype)_bg_$(10. ^σ0)Spm", 
    preferEright=true, showplot=true, logscale=true, fontsize=9, yl=[-250,80])