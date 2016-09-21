

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
    quadgk(rp->2.*π.*rp.*tiso_Ψ(a, M200, rp), 0.0, rp)[1]
end

#############################
#
#   Sharply truncated NFW (@ R200)
#
###################################################################

function tnfw_m(a, M200, RS, r)
    if r<=0
        r=1.0e-3
    end
    R=R200(a, M200)
    rr200 = RS+R
    k=0.0
    try
        k = M200./(log(rr200./RS)-R./rr200)
    catch
        @show rr200, RS, R
        error("xxxxxxxx")
    end
    r = min(R, r)
    rrs = RS+r
    (log(rrs./RS)-r./rrs).*k
end

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

#############################
#
#   Projecting a general density
#
###################################################################

function sphere_pick()
    u = rand()
    v = rand()
    θ = 2*π*u
    ϕ = acos(2*v-1)
    x = cos(θ)*sin(ϕ)
    y = sin(θ)*sin(ϕ)
    z = cos(θ)
    x,y,z
end

function prepare_proj_mass(massfunc, a, N=10000)
    M200 = massfunc(1.0e30)
    R = R200(a, M200)
    pm = M200/N # infinity ha!
    cm = 0.0
    rp = Float64[]
    for i in 1:N-1
        # find next radius for particle
        cm += pm
        nr = fzero(r->massfunc(r)-cm, 1.0e-10,R)
        nx, ny, nz = sphere_pick()
        push!(rp, nr*sqrt(nx*nx + ny*ny))
    end
    m = pm.*collect(1:N-1)
    m = m*M200/m[end]
    rp = sort!(rp)
    rp = rp*R/rp[end]
    Spline1D(rp, m; k=5, bc="nearest"), rp, m
end

function project_mass(spl, rp)
    evaluate(spl, rp)
end


