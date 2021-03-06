
function thin(proj_massfunc, θ, rs, tgα, N=100; ext=1.0)
    const cθ = cos(θ)
    const sθ = sin(θ)    
    const x0 = -rs*cθ
    const y0 = rs*sθ
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
    x = [reverse(lx).*x0; lx.*nx.*ext]
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

    x, y, u, du, in_angle, out_angle, deflection_angle, π-atan2(y,x)
end

# some convenience functions:

thin_tiso200b(z, M200b, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tiso_mp(M200b, r, R200b(z, M200b)), θ, rs, tgα, N; ext=ext)
thin_tiso200c(z, M200c, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tiso_mp(M200c, r, R200c(z, M200c)), θ, rs, tgα, N; ext=ext)
thin_tisoVir(z, Mvir, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tiso_mp(Mvir, r, RVIR(z, MVir)), θ, rs, tgα, N; ext=ext)

thin_tnfw200b(z, M200b, Rscale, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tnfw_mp(M200b, Rscale, r, R200b(z, M200b)), θ, rs, tgα, N; ext=ext)
thin_tnfw200c(z, M200c, Rscale, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tnfw_mp(M200c, Rscale, r, R200c(z, M200c)), θ, rs, tgα, N; ext=ext)
thin_tnfwVir(z, Mvir, Rscale, θ, rs, tgα, N=100; ext=1.0) =
    thin(r->tnfw_mp(Mvir, Rscale, r, RVIR(z, Mvir)), θ, rs, tgα, N; ext=ext)
