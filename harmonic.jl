
function harmonic(a, M200, θ0, tgα0, N=100)
    R=R200(a,M200)
    ω=sqrt(1.0-3M200*G/C/C/R)
    
    dudθ0 = dudθ(θ0, tgα0, R)
    u0 = 1.0/R
    
    ϕ = atan(-dudθ0/ω/u0)-ω*θ0
    A = u0/cos(ω*θ0+ϕ)
    
    ϕm=π/2-ω*ϕ*2.0
    
    ϕl = linspace(θ0, ϕm, N)
    u = A*cos(ω*ϕl+ϕ)
    du = -A*ω*sin(ω*ϕl+ϕ)

    ixs = find(u.>u0)
    u=u[ixs]
    ϕl = ϕl[ixs]
    
    x = [Float64(v) for v in -cos(ϕl)./u]
    y = [Float64(v) for v in sin(ϕl)./u]
        
    out_angle =  angle( ϕl[end], 1/u[end], du[end])
    in_angle =  atan(tgα0)

    x, y, u, du, in_angle, -out_angle, out_angle+in_angle
end

