
function time(phifunc, x, y)
    dx = x[2:end] - x[1:(end-1)]
    dy = y[2:end] - y[1:(end-1)]
    dpath = sqrt(dx.*dx + dy.*dy)
    phi = zeros(length(dx))
    C2 = C*C
    @inbounds for i in 1:length(dx)
        xc = x[i] + 0.5*dx[i] 
        yc = y[i] + 0.5*dy[i] 
        phi[i] = phifunc(xc,yc)./C2
    end
    time_dilation = 1.0-2.*phi
    sum(dpath./time_dilation) ./ C
end
