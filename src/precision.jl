function get_delay_hrs(NBITS, b, db, rs, prec, el=true)
    setprecision(NBITS)
    prec = 1.0e-14

    massfunc(r) = tnfw_m(0.1, 500.0, 20.0, r);
    rhofunc(r) = tnfw_ρ(0.1, 500.0, 20.0, r);
    tg = b/rs
    up_ex, up_ey, up_eu, up_edu, up_eia, up_eoa, up_eda, up_eθ =
    el? numeric_el(massfunc, rhofunc, 0.0, rs, tg, 100000; reltol=prec, abstol=prec) :
        numeric(massfunc, 0.0, rs, tg, 100000; reltol=prec, abstol=prec);
    b += db
    tg = b/rs
    down_ex, down_ey, down_eu, down_edu, down_eia, down_eoa, down_eda, down_eθ =
    el? numeric_el(massfunc, rhofunc, 0.0, rs, tg, 100000; reltol=prec, abstol=prec) :
        numeric(massfunc, 0.0, rs, tg, 100000; reltol=prec, abstol=prec);

    uex,uey,dex,dey = extrapolate_to_meeting(up_ex,up_ey,down_ex,down_ey);

    uet = time(massfunc, uex, uey);
    det = time(massfunc, dex, dey);
    delta_e = abs(uet-det);
    display(string(NBITS)*", "*string(db)*", "*string(delta_e))
    delta_e*1e9*365*24;
end
