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

PRO run_n2vk_model, Ts, wl3, sp3
    
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
    convolve_sp, wl(wv_vk,6,21), int(int_vk,6,21), 0.6d, wl3, sp6_21
    ;20 5
    convolve_sp, wl(wv_vk,5,20), int(int_vk,5,20), 0.6d, wl3, sp5_20
    ;18 2
    convolve_sp, wl(wv_vk,2,18), int(int_vk,2,18), 0.6d, wl3, sp2_18
    ;17 1
    convolve_sp, wl(wv_vk,1,17), int(int_vk,1,17), 0.6d, wl3, sp1_17
    ;16 0
    convolve_sp, wl(wv_vk,0,16), int(int_vk,0,16), 0.6d, wl3, sp0_16

    sp3 = [[sp6_21], [sp5_20], [sp2_18], [sp1_17], [sp0_16]]
END