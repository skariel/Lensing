
function numeric(massfunc, θ, rs, tgα, N=100; abstol=1.0e-13, reltol=1.0e-14)
    tgα = -tgα
    ϕm=π
    ϕl = linspace(θ, π, N)    
    function f(ϕ, y)
        (u, v) = y
        u_prime = v
        v_prime = 3*rg(massfunc(1./u)).*u.*u-u
        [u_prime; v_prime]    
    end
    u0 = 1/rs
    du0 = dudθ(θ, tgα, rs)
    start = [u0; du0];
    
    ϕ, y = ode45(f, start, ϕl, abstol=abstol, reltol=reltol);
    u = map(y -> y[1], y);
    du = map(y -> y[2], y);
    ixs = find(u.>0)
    u=u[ixs]
    du=du[ixs]
    ϕ=ϕ[ixs]
    
    x = [Float64(v) for v in -cos(ϕ)./u]
    y = [Float64(v) for v in sin(ϕ)./u]

    out_angle =  ray_angle( ϕl[end], 1/u[end], du[end])
    in_angle =  atan(tgα)

    x, y, u, du, in_angle, -out_angle, out_angle+in_angle
end

numeric_tiso(a, M200, θ, rs, tgα, N=100; abstol=1.0e-13, reltol=1.0e-14) =
    numeric(r->tiso_m(a, M200, r), θ, rs, tgα, N; abstol=abstol, reltol=reltol) 