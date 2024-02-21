pro stuff

; Continuum 550-600
; Background 250-300

continuum_dseq = dseq[550:600, *, *]

dseq = reform(total(dseq, 1), [1, 512, 512])

plot1 = plot(wl_3, continuum_sp_3, title='Continuum vs Bkg, 3Jan2020, Ha panel', xtitle='Ang', name='continuum sp ha')
plot2 = plot(wl_3, background_sp_3, /overplot, color='red', name='background sp Ha')
leg = legend(target=[plot1,plot2], /auto_text_color )

end
