transD_GP.plotconvandlast(soundings, 50, 2; zstart=zstart, extendfrac=extendfrac, dz=dz, nlayers=nlayers, 
    figsize=(14,7), vmin=lo, vmax=hi, postfix="_$(regtype)_β²_$(β²)_bg_$(10. ^σ0)Spm", preferEright=true, 
    showplot=true, logscale=true, prefix=zipsaveprefix)
