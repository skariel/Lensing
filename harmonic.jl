
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
    
    x = [Float64(v) for v in -cos(ϕl)./u]
    y = [Float64(v) for v in sin(ϕl)./u]
        
    out_angle =  ray_angle( ϕl[end], 1/u[end], du[end])
    in_angle =  atan(tgα0)

    x, y, u, du, -in_angle, out_angle, -out_angle-in_angle, ϕl
end

function harmonic_e(a, M200, rmin, θ0, N=100)
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

    dx=hex[2]-hex[1]
    dy=hey[2]-hey[1]
    a1 = atan2(dy,dx)
    dx=ex[2]-ex[1]
    dy=ey[2]-ey[1]
    a2 = atan2(dy,dx)

    @show dα = π/2+a1-a2
    θ = ϕ+dα

    x = [Float64(v) for v in -cos(θ)./u]
    y = [Float64(v) for v in sin(θ)./u]

    x, y, u, θ
end

