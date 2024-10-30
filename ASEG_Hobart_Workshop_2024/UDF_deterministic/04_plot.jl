lnames = [100401, 100502]
idx = [[50, 200],[300]]
dr = 12.5
transD_GP.plotconvandlast(soundings, dr, dz; zall, lnames, idx, 
  plotforward=true, aem_in=aem, prefix=zipsaveprefix,
  figsize=(14,5), vmin=lo, vmax=hi, postfix=zipsaveprefix, yl=[-150, 100],
  preferEright=true, showplot=true, logscale=true, saveplot=true)

