rg(M) = M.*G./C./C

function dudθ(θ, tgα, R)
    cθ = cos(θ)
    sθ = sin(θ)
    (cθ+sθ.*tgα)./R./(sθ-cθ.*tgα) 
end

function dduddθ(θ, tgα, R)
    cθ = cos(θ)
    sθ = sin(θ)

#        V                 W
#    ------------      -----------
#    (cθ+sθ.*tgα)./R./(sθ-cθ.*tgα) 

    v = cθ+sθ.*tgα
    w = sθ-cθ.*tgα

    dv = -sθ+cθ.*tgα
    dw =  cθ-sθ.*tgα

#   1/R*{  V*W^(-1) }'
    1/R*( dv/w-v/w/w*dw)
end

function ray_angle(θ, r, dudθ)
    cθ = cos(θ)
    sθ = sin(θ)
    atan((-r.*dudθ.*sθ+cθ)./(-r.*dudθ.*cθ-sθ))
end
