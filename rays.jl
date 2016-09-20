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
    -atan((-r.*dudθ.*sθ+cθ)./(-r.*dudθ.*cθ-sθ))
end


function mylinspace(i,f,N)
    l = linspace(0,1,div(N,2)).^0.001
    ll = π/2-i
    lr = f-π/2
    [(1-reverse(l)).*ll+i; l.*lr+π/2]
end