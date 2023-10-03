pro meplot, x, y, z, range, xtitle, ytitle, title

xran = [min(x), max(x)]
yran = [min(y), max(y)]
layout, 1, pos, cpos

!p.position = pos
plot, fltarr(3), /nodata, xran=xran, yran=yran, /xsty, /ysty, xtitle=xtitle, ytitle=ytitle, title=title
plot_panel, float(x), float(y), float(z), range=range

!p.position = cpos
cbar, range=range

end
