rg(M) = M.*G./C./C

function dudθ(θ, tgα, R)
    cθ = cos(θ)
    sθ = sin(θ)
    (cθ+sθ.*tgα)./R./(sθ-cθ.*tgα) 
end

function ray_angle(θ, r, dudθ)
    cθ = cos(θ)
    sθ = sin(θ)
    atan((-r.*dudθ.*sθ+cθ)./(-r.*dudθ.*cθ-sθ))
end
