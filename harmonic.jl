
function harmonic(a, M200, θ0, tgα0, N=100)
    tgα0 = -tgα0
    R=R200(a,M200)
    ω=sqrt(1.0-3*M200*G/C/C/R)
    
    dudθ0 = dudθ(θ0, tgα0, R)
    u0 = 1.0/R
    
    ϕ = atan(-dudθ0/ω/u0)-ω*θ0
    A = u0/cos(ω*θ0+ϕ)
    
    ϕm=π/2-ω*ϕ*2.0
    
    ϕl = linspace(θ0, ϕm, N)
    u = A*cos(ω*ϕl+ϕ)
    du = -A*ω*sin(ω*ϕl+ϕ)

    ixs = find(u.>=u0)
    u=u[ixs]
    du = du[ixs]
    ϕl = ϕl[ixs]
    
    x = [v for v in -cos(ϕl)./u]
    y = [v for v in sin(ϕl)./u]
        
    out_angle =  ray_angle( ϕl[end], 1/u[end], du[end])
    in_angle =  atan(tgα0)

    x, y, u, du, -in_angle, out_angle, -out_angle-in_angle, ϕl
end

function harmonic_e(a, M200, rmin, tgα, N=100)
    tgα =- tgα
    R=R200(a,M200)
    m0 = M200*G/C/C
    ω=sqrt(1.0-3*m0/R)
    
    ϕ = linspace(-π/2, π/2, N)
    α = ω.*ϕ
    cα = cos(α)
    sα = sin(α)

    u = +cα/rmin + m0/R/rmin/ω/ω * (α.*sα+cα.*log(cα))

    ixs = find(u.>((1./R)))
    u=u[ixs]
    ϕ=ϕ[ixs]

    a1 = atan(tgα)
    dx = -cos(ϕ[2]+π/2)./u[2] + cos(ϕ[1]+π/2)./u[1]
    dy =  sin(ϕ[2]+π/2)./u[2] - sin(ϕ[1]+π/2)./u[1]
    a2 = -atan2(dy,dx)

    @show dα = π/2+a1-a2
    θ = ϕ+dα
    @show dα-π/2

    x = [v for v in -cos(θ)./u]
    y = [v for v in sin(θ)./u]

    # TODO: du, outgoing_angle, deflection_angle
    in_angle = atan(tgα)
    x, y, u, θ, -in_angle
end

