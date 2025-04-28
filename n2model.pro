FUNCTION wl, wv_vk, v_lower, v_upper
    a = wv_VK[0, v_lower, v_upper, *]
    RETURN, reform(a, n_elements(a))
END

FUNCTION int, int_vk, v_lower, v_upper
    a = int_VK[0, v_lower, v_upper, *]
    RETURN, reform(a, n_elements(a))
END

FUNCTION to_nm, arr
    RETURN, arr/10.
END

PRO run_n2vk_model, Ts, wl1, wl2, wl3, sp1, sp2, sp3

    fwhm = 0.06d
    
    prog_n2_VKband_v1, Ts, wv_VK, int_VK

    ; Format
    ;  wv_VK[Tindex, v_lower, v_upper, Jindex]
    ; int_VK[Tindex, v_lower, v_upper, Jindex]

    read_tim, '30/12/2021 08:00:00', 1/60., mjs0, time, dseq, icount, /nophot
    get_w, mjs0, 3, wl3
    get_w, mjs0, 2, wl2
    get_w, mjs0, 1, wl1

    wl3 = to_nm(wl3)
    wl2 = to_nm(wl2)
    wl1 = to_nm(wl1)

    ; For Halpha panel
    ;21 6
    convolve_sp, wl(wv_vk,21,6), int(int_vk,21,6), fwhm, wl3, sp6_21
    ;20 5
    convolve_sp, wl(wv_vk,20,5), int(int_vk,20,5), fwhm, wl3, sp5_20
    ;18 2
    convolve_sp, wl(wv_vk,18,2), int(int_vk,18,2), fwhm, wl3, sp2_18
    ;17 1
    convolve_sp, wl(wv_vk,17,1), int(int_vk,17,1), fwhm, wl3, sp1_17
    ;16 0
    convolve_sp, wl(wv_vk,16,0), int(int_vk,16,0), fwhm, wl3, sp0_16

    sp3 = [[sp6_21], [sp5_20], [sp2_18], [sp1_17], [sp0_16]]

    ; For Oplus panel
    ;21 6
    convolve_sp, wl(wv_vk,21,5), int(int_vk,21,5), fwhm, wl2, sp5_21
    ;20 4
    convolve_sp, wl(wv_vk,20,4), int(int_vk,20,4), fwhm, wl2, sp4_20
    ;18 1
    convolve_sp, wl(wv_vk,18,1), int(int_vk,18,1), fwhm, wl2, sp1_18

    sp2 = [[sp5_21], [sp4_20], [sp1_18]]

    ; For OH panel
    ;21 4
    convolve_sp, wl(wv_vk,21,4), int(int_vk,21,4), fwhm, wl1, sp4_21
    ;20 3
    convolve_sp, wl(wv_vk,20,3), int(int_vk,20,3), fwhm, wl1, sp3_20
    ;17 0
    convolve_sp, wl(wv_vk,17,0), int(int_vk,17,0), fwhm, wl1, sp0_17

    sp1 = [[sp4_21], [sp3_20], [sp0_17]]
END