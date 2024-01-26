pro get_timings, fname, years, months, days, hours, mins, secs, msecs, lens

p = read_csv(fname)
si = size(p.(0))
len = si[-1]

years = p.(0)
months = p.(1)
days = p.(2)
hours = p.(3)
mins = p.(4)
secs = p.(5)
msecs = p.(6)

lens = list([10], length = len)

end


pro spectrathing, years, months, days, hours, mins, secs, msecs, lens, spectrum

totaldseq = 0
totalchonks = 0
a = size(years)
len = a[1]

for i = 0, len-10 do begin
	getchonks, years[i], months[i], days[i], hours[i], mins[i], secs[i], msecs[i], 10, lens[i], startimes
	b = size(startimes)
	leng = b[1]
	for j=0,leng-1 do begin
		read_tim, startimes[j], 10/60., mjs0, time, dseq, icount, /nophot, tadd=10*60
		dseq = reform(total(dseq, 1), [1, 512, 512])
		totaldseq += dseq
		totalchonks += 1
	endfor
endfor

spectra, 3, mjs0, time, totaldseq/totalchonks, spectrum

end


; Can I just add dseq? does it make a difference?
; Examples 10/12/21 07:00:00 and 12/12/21 07:10:00
