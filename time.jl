
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
