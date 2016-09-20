

#############################
#
#   Sharply truncated ISOTHERMAL (@ R200)
#
###################################################################

function tiso_m(a, M200, r)
    R=R200(a, M200)
    r = min(R, r)
    M200.*r./R
end

function tiso_ρ(a, M200, r)
    R=R200(a, M200)
    if r>R
        return zero(r)
    end
    M200./4./π./R./r./r
end

function tiso_Ψ(a, M200, rp)
    R = R200(a, M200)
    if rp>R
        return zero(rp)
    end
    M200./2./π./R./rp.*atan(sqrt(R.*R./rp./rp-1))
end

function tiso_mp(a, M200, rp)
    R = R200(a, M200)
    rp = min(R, rp)
    @show R, rp
    quadgk(rp->2.*π.*rp.*tiso_Ψ(a, M200, rp), 0.0, rp)[1]
end



