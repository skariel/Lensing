#############################
#
#   Some general functions so we furtehr only need 3D density
#
###################################################################

# a general function for cumulative projected mass
function projected_mass(rhofunc, rp, RMAX)
    rp = min(rp, RMAX)
    rp = max(rp, 1.0e-37)
    function circ_rho_func(r)
        θ = r<rp ? π/2 : asin(rp/r)
        rhofunc(r)*4*π*r*r*(1.0-cos(θ))
    end
    quadgk(circ_rho_func, 1.0e-37, RMAX)[1]
end

# a general function for cumulative mass
function cumulative_mass(rhofunc, r, RMAX)
    r = min(r, RMAX)
    r = max(r, 1.0e-37)
    quadgk(r->4.0*π*r*r*rhofunc(r), 1.0e-37, r)[1]
end

#############################
#
#   Sharply truncated ISOTHERMAL (@ R200)
#
###################################################################

function tiso_ρ(a, M200, r)
    R=R200(a, M200)
    if r>R
        return zero(r)
    end
    M200./4./π./R./r./r
end

tiso_m(a, M200, r) =
    cumulative_mass(r->tiso_ρ(a, M200, r), r, R200(a, M200))

tiso_mp(a, M200, rp) =
    projected_mass(r->tiso_ρ(a, M200, r), rp, R200(a, M200))

#############################
#
#   Sharply truncated NFW (@ R200)
#
###################################################################

function tnfw_ρ(a, M200, RS, r)
    R=R200(a, M200)
    if r>R
        return zero(r)
    end
    rr200 = RS+R
    k = M200./(log(rr200./RS)-R./rr200)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2)
end

tnfw_m(a, M200, RS, r) =
    cumulative_mass(r->tnfw_ρ(a, M200, RS, r), r, R200(a, M200))

tnfw_mp(a, M200, RS, rp) =
    projected_mass(r->tnfw_ρ(a, M200, RS, r), rp, R200(a, M200))

#############################
#
#   Exponentially truncated NFW (@ R200)
#
###################################################################

function lnfw_ρ(a, M200, RS, r)
    R=R200(a, M200)
    rr200 = RS+R
    k = M200./(log(rr200./RS)-R./rr200)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2) ./ (1.0+exp((r-R)./RS))
end

lnfw_m(a, M200, RS, r) =
    cumulative_mass(r->lnfw_ρ(a, M200, RS, r), r, 100*R200(a, M200))

lnfw_mp(a, M200, RS, rp) =
    projected_mass(r->lnfw_ρ(a, M200, RS, r), rp, 100*R200(a, M200))

#############################
#
#   Void surrounded NFW (@ R200)
#
###################################################################

function vnfw_ρ(a, M200, RS, r; vstart=2.0, vend=3.0)
    R=R200(a, M200)
    rr200 = RS+R
    k = M200./(log(rr200./RS)-R./rr200)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2) ./ (1.0+exp((r-R)./RS)) + RHO_CRIT*OM(a).*(-1./(1.0+exp(-(r-vstart*R)./RS)) ./ (1.0+exp((r-vend*R)./RS))) 
end

vnfw_m(a, M200, RS, r; vstart=2.0, vend=10.0) =
    cumulative_mass(r->vnfw_ρ(a, M200, RS, r; vstart=vstart, vend=vend), r, 100*R200(a, M200))

vnfw_mp(a, M200, RS, rp; vstart=2.0, vend=10.0) =
    projected_mass(r->vnfw_ρ(a, M200, RS, r; vstart=vstart, vend=vend), rp, 100*R200(a, M200))

