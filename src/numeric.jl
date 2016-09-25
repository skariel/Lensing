

function _numeric_f(integrand, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17, use78=false)
    tgα = -tgα
    ϕm=π
    ϕl = mylinspace(θ, π, N)    
    u0 = 1/rs
    du0 = dudθ(θ, tgα, rs)
    start = [u0; du0];
    
    ϕ, y = use78?
        ode78(integrand, start, ϕl, abstol=abstol, reltol=reltol) :
        ode45(integrand, start, ϕl, abstol=abstol, reltol=reltol);
    u = map(y -> y[1], y);
    du = map(y -> y[2], y);
    ixs = find(u.>0)
    u=u[ixs]
    du=du[ixs]
    ϕ=ϕ[ixs]
    
    x = [v for v in -cos(ϕ)./u]
    y = [v for v in sin(ϕ)./u]

    out_angle =  ray_angle(ϕ[end], 1/u[end], du[end])
    in_angle =  atan(tgα)

    x, y, u, du, -in_angle, out_angle, -out_angle-in_angle, ϕ
end

function numeric(massfunc, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17)
    function f(ϕ, y)
        (u, v) = y
        u_prime = v
        v_prime = 3*rg(massfunc(1./u)).*u.*u-u
        [u_prime; v_prime]    
    end
    _numeric_f(f, θ, rs, tgα, N; abstol=abstol, reltol=reltol)
end

function numeric_el(massfunc, rhofunc, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17)
    M200 = massfunc(1.0e30) # infinity, ha!
    du = dudθ(θ, tgα, rs)
    u = 1/rs
    _rg = rg(M200)
    _e = 1/(1-2*_rg*u)
    e2l2 = du*du + u*u*_e
    function integrand(ϕ, y)
        (u, v) = y
        _rg = rg(massfunc(1./u))
        u_prime = v
        v_prime = 3*_rg.*u.*u-u+e2l2*4*π*G*rhofunc(1./u)/C/C/u/u/u
        [u_prime; v_prime]    
    end
    _numeric_f(integrand, θ, rs, tgα, N; abstol=abstol, reltol=reltol, use78=true)
end

# some convenience functions:

numeric_tiso200b(z, M200b, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tiso_m(M200b, r, R200b(z, M200b)), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tiso200c(z, M200c, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tiso_m(M200c, r, R200c(z, M200c)), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tiso_el200b(z, M200b, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tiso_m(M200b, r, R200b(z, M200b)), r->tiso_ρ(M200b, r, R200b(z, M200b)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tiso_el200c(z, M200c, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tiso_m(M200c, r, R200c(z, M200c)), r->tiso_ρ(M200c, r, R200c(z, M200c)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tisoVir(z, Mvir, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tiso_m(Mvir, r, RVIR(z, MVir)), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tiso_elVir(z, Mvir, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tiso_m(Mvir, r, RVIR(z, MVir)), r->tiso_ρ(Mvir, r, RVIR(z, MVir)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 

numeric_tnfw200b(z, M200b, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tnfw_m(M200b, Rcalse, r, R200b(z, M200b)), θ, rs, tgα, N; abstol=abstol, reltol=reltol)
numeric_tnfw200c(z, M200c, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tnfw_m(M200c, Rscale, r, R200c(z, M200c)), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tnfw_el200b(z, M200b, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tnfw_m(M200b, Rscale, r, R200b(z, M200b)), r->tnfw_ρ(M200b, Rscale, r, R200b(z, M200b)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tnfw_el200c(z, M200c, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tnfw_m(M200c, Rscale, r, R200c(z, M200c)), r->tnfw_ρ(M200c, Rscale, r, R200c(z, M200c)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tnfwVir(z, Mvir, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric(r->tnfw_m(Mvir, Rscale, r, RVIR(z, MVir)), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
numeric_tnfw_elVir(z, Mvir, Rscale, θ, rs, tgα, N=100; abstol=1.0e-16, reltol=1.0e-17) =
    numeric_el(r->tnfw_m(Mvir, Rscale, r, RVIR(z, MVir)), r->tnfw_ρ(Mvir, Rscale, r, RVIR(z, MVir)),
                θ, rs, tgα, N; abstol=abstol, reltol=reltol) 
