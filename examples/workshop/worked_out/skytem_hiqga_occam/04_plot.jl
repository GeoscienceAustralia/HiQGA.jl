transD_GP.SkyTEM1DInversion.splitlineconvandlast(soundings, 5, 2; zstart=zstart, extendfrac=extendfrac, dz=dz, nlayers=nlayers, 
    figsize=(12,7), vmin=lo, vmax=hi, postfix="_β²_$(β²)_$(regtype)_bg_$(10. ^σ0)Spm", 
    preferEright=true, showplot=true, logscale=true)
@info "writing PDF ..."
run(`./pdfscript.sh`)
@info "done"