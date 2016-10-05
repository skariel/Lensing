function fermat(phifunc_xy, forcefunc_xy, rs, tgα, N=100; reltol=1e-21, abstol=1e-21, ext=1.0)
    const C2 = C*C
    function f(x,yvec)
        (y, v) = yvec
        y_prime = v
        v_prime = 0.0

        phi = phifunc_xy(x,y)/C/C
        dphidx, dphidy = forcefunc_xy(x,y)
        dphidx /= C2
        dphidy /= C2
        v2 = v*v

        v_prime = 2*(dphidx*(v2*v+v)-dphidy*(v2+1)^2) / (1-2*phi)
        [y_prime, v_prime]
    end
    _x = linspace(0.0,1,div(N,2)).^2
    _l = [-reverse(_x)*rs ; _x*rs.*ext]
    x, u = ode45(f, [0.0; tgα], _l; reltol=reltol, abstol=abstol)
    y = map(u -> u[1], u);
    dy = map(u -> u[2], u);    
    x,y
end

function fermat_y(forcefunc_y, rs, tgα, N=100; reltol=1e-21, abstol=1e-21, ext=1.0)
    const C2 = C*C
    function f(x,yvec)
        (y, v) = yvec
        y_prime = v
        v_prime = 0.0

        v_prime = -2*forcefunc_y(x,y)/C2
        [y_prime, v_prime]
    end
    _x = linspace(0.0,1,div(N,2)).^2
    _l = [-reverse(_x)*rs ; _x*rs.*ext]
    x, u = ode45(f, [0.0; tgα], _l; reltol=reltol, abstol=abstol)
    y = map(u -> u[1], u);
    dy = map(u -> u[2], u);    
    x,y
end
