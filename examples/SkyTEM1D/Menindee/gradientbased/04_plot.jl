transD_GP.plotconvandlast(soundings, 50, 2; zall, 
    figsize=(15,7), vmin=lo, vmax=hi, postfix="_β²_$(β²)_$(regtype)_bg_$(round(10. ^σ0, sigdigits=4))Spm", prefix=zipsaveprefix,
    preferEright=true, showplot=true, logscale=true)
@info "writing PDF ..."
run(`./pdfscript.sh`)
@info "done"
