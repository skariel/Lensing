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

function extrapolate_to_meeting(x1,y1,x2,y2)
    dx1 = x1[end]-x1[end-10]
    dy1 = y1[end]-y1[end-10]
    tg1 = dy1/dx1
    dx2 = x2[end]-x2[end-10]
    dy2 = -(y2[end]-y2[end-10])
    tg2 = dy2/dx2

    DX2 = (-y1[end-10]-x2[end-10]*tg1+x1[end-10]*tg1-y2[end-10])./(tg1-tg2)
    X = x2[end-10] + DX2
    Y = -y2[end-10] + DX2*tg2
    
    xx1 = x1[x1.<X]
    yy1 = y1[x1.<X]
    xx2 = x2[x2.<X]
    yy2 = y2[x2.<X]
    push!(xx1, X)
    push!(xx2, X)

    push!(yy1, Y)
    push!(yy2, -Y)

    xx1,yy1, xx2,yy2
end

function extrapolate_to_meeting_same_side(x1,y1,x2,y2)
    dx1 = x1[end]-x1[end-10]
    dy1 = y1[end]-y1[end-10]
    tg1 = dy1/dx1
    dx2 = x2[end]-x2[end-10]
    dy2 = y2[end]-y2[end-10]
    tg2 = dy2/dx2

    DX2 = (-y1[end-10]-x2[end-10]*tg1+x1[end-10]*tg1+y2[end-10])./(tg1-tg2)
    X = x2[end-10] + DX2
    Y = y2[end-10] + DX2*tg2
    
    xx1 = x1[x1.<X]
    yy1 = y1[x1.<X]
    xx2 = x2[x2.<X]
    yy2 = y2[x2.<X]
    push!(xx1, X)
    push!(xx2, X)

    push!(yy1, Y)
    push!(yy2, Y)

    xx1,yy1, xx2,yy2
end

function extrapolate_to_x2(x1,y1,x2)
    dx1 = x1[end]-x1[end-10]
    dy1 = y1[end]-y1[end-10]
    tg1 = dy1/dx1

    xx1 = x1[x1.<x2]
    yy1 = y1[x1.<x2]

    DX1 = x2-x1[end]
    DY1 = DX1*tg1

    X = x2
    Y = y1[end] + DY1

    push!(xx1, X)
    push!(yy1, Y)

    xx1,yy1
end

function to_small(v)
    Float64[float(vi) for vi in v]
end