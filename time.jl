
function time(massfunc, x, y)
    dx = x[2:end] - x[1:(end-1)]
    dy = y[2:end] - y[1:(end-1)]
    dpath = sqrt(dx.*dx + dy.*dy)
    r = sqrt(x.*x + y.*y)
    rm = 0.5*(r[2:end] + r[1:(end-1)])
    m = [massfunc(ri) for ri in rm]
    time_dilation = sqrt(1.0-2.*G.*m./rm./C./C)
    sum(dpath./time_dilation) ./ C
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