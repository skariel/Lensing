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

# a general function for potential
function potential(massfunc, r, RMAX)
    r = max(r, 1.0e-37)
    quadgk(r->G*massfunc(r)./r./r, 1.0e37, r)[1]
end

#############################
#
#   Sharply truncated ISOTHERMAL (@ R200)
#
###################################################################

function tiso_ρ(MX, r, RMAX)
    if r>RMAX
        return zero(r)
    end
    MX./4./π./RMAX./r./r
end

function tiso_m(MX, r, RMAX)
    r = min(r, RMAX)
    MX./RMAX.*r
end

tiso_mp(MX, rp, RMAX) =
    projected_mass(r->tiso_ρ(MX, r, RMAX), rp, RMAX)

#############################
#
#   Sharply truncated NFW (@ RX)
#
###################################################################

function tnfw_ρ(MX, RS, r, RMAX)
    RMAX
    if r>RMAX
        return zero(r)
    end
    rrX = RS+RMAX
    k = MX./(log(rrX./RS)-RMAX./rrX)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2)
end

function tnfw_m(MX, RS, r, RMAX)
    if r<1.0e-17
        return zero(r)
    end 
    rrX = RS+RMAX
    k = MX./(log(rrX./RS)-RMAX./rrX)
    r = min(RMAX, r)
    rrs = RS+r
    (log(rrs./RS)-r./rrs).*k
end

tnfw_potential(MX, RS, r, RMAX) =
    potential(r->tnfw_m(MX,RS,r,RMAX), r, RMAX)

tnfw_mp(MX, RS, rp, RMAX) =
    projected_mass(r->tnfw_ρ(MX, RS, r, RMAX), rp, RMAX)

#############################
#
#   Exponentially truncated NFW (@ RX)
#
###################################################################

function lnfw_ρ(MX, RS, r, RMAX)
    rrX = RS+RMAX
    k = MX./(log(rrX./RS)-RMAX./rrX)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2) ./ (1.0+exp((r-RMAX)./RS))
end

lnfw_m(MX, RS, r, RMAX) =
    cumulative_mass(r->lnfw_ρ(MX, RS, r, RMAX), r, 100*RMAX)

lnfw_mp(MX, RS, rp, RMAX) =
    projected_mass(r->lnfw_ρ(MX, RS, r, RMAX), rp, 100*RMAX)

#############################
#
#   Void surrounded NFW
#
###################################################################

function vnfw_ρ(MX, RS, r, RMAX; vstart=2.0, vend=3.0)
    rrX = RS+RMAX
    k = MX./(log(rrX./RS)-RMAX./rrX)
    ρ0 = k/4/π/RS/RS/RS
    x = r./RS
    ρ0./(x.*(1.0+x).^2) ./ (1.0+exp((r-RMAX)./RS)) + RHO_CRIT*OM(a).*(-1./(1.0+exp(-(r-vstart*RMAX)./RS)) ./ (1.0+exp((r-vend*RMAX)./RS))) 
end

vnfw_m(MX, RS, r, RMAX; vstart=2.0, vend=10.0) =
    cumulative_mass(r->vnfw_ρ(MX, RS, r, RMAX; vstart=vstart, vend=vend), r, 100*RMAX)

vnfw_mp(MX, RS, rp, RMAX; vstart=2.0, vend=10.0) =
    projected_mass(r->vnfw_ρ(MX, RS, r, RMAX; vstart=vstart, vend=vend), rp, 100*RMAX)

