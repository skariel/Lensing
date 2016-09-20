function thin(proj_massfunc, θ, rs, tgα, N=100)
    const cθ = cos(θ)
    const sθ = sin(θ)    
    const x0 = -rs*cθ
    @show const y0 = rs*sθ
    const b = y0-tgα*x0
    const α=atan(tgα)
    const M0p = proj_massfunc(b)
    const deflection_angle = 4*M0p*G/C/C/b
    const nα = α-deflection_angle
    const tgnα = tan(nα)
    nx = -b/tgnα
    if nα >0.0
        nx = rs
    end
    
    lx = linspace(0.0,1.0,div(N,2)).^3
    x = [reverse(lx).*x0; lx.*nx]
    y = zeros(length(x))
    for i in 1:length(x)
        y[i] = if x[i]<0.0
            y0+tgα*(x[i]-x0)
        else
            b+tgnα*x[i]
        end
    end

    in_angle =  atan(tgα)
    out_angle =  nα

    u = 1./(x.*x+y.*x)
    θ = atan2(y, x)
    du = dudθ(θ, tgα, 1./u) # TODO: test this du is probably wrong, and not currently used

    x, y, u, du, in_angle, out_angle, deflection_angle
end

thin_tiso(a, M200, θ, rs, tgα, N=100) =
    thin(r->tiso_mp(a, M200, r), θ, rs, tgα, N) 