pro trim_blueshifts, intensity, blueshift, cutoff, trimmed
; Remove all the blueshifts where there isn't any proton aurora
; Inputs
;   intensity FLOAT array : Intensity values for the duration
;   blueshift "           : Blueshift          " 
;   cutoff FLOAT          : The divider between proton aurora and no proton aurora, based on the below formula or chosen manually
; Outputs
;   trimmed FLOAT array   : The blueshift array with any values below the cutoff replaced with NaNs

ind = where(intensity le cutoff)

trimmed = blueshift
trimmed[ind] = !VALUES.F_NAN

end

pro whats_cutoff_precious, intensity, cutoff
; Returns the cutoff at 20% over the no aurora intensity

; Lowest value
low = min(intensity)

; cutoff at 20% over the lowest intensity (roughly no aurora)
cutoff = low * 1.2

end

pro plot_all_in_time, intensity, blueshift, temperature

whats_cutoff_precious, intensity, cutoff
trim_blueshifts, intensity, blueshift, cutoff, trimmed

dims = size(intensity)
dim = dims[-1]

time = indgen(dim)

!x.margin=[8, 9]
   ; plot with no right axis
plot, time, intensity, ystyle=8, ytitle='_____ Intensity'

plot, time, temperature, /noerase, xrange=!x.crange, line=2, ystyle=4, xstyle=1

   ; make an axis on the right for the last plot
AXIS, YAXIS=1, YRANGE=!y.crange, ytitle='----- Temperature'

plot, time, trimmed, /noerase, xrange=!x.crange, line=2, ystyle=4, xstyle=1

end