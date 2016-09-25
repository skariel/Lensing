# using Gadget2 standard units:
#     velocity is canonical momentum, in km/s
#     distance is in kpc/h
#     mass is in 1e10 Ms/h
const G = 43007.1
const H0 = 0.1
const C = 299792.458

# using Dark Sky cosmology:
const ΩΛ = 0.7048
const Ω0 = 1 - ΩΛ
const h = 0.6881
const σ8 = 0.8344

const RHO_CRIT_TODAY = 3*H0^2/8/π/G

# Hubble as function of scale factor a, or z
H(z) = H0*sqrt(Ω0.*(z+1.0).^3+ΩΛ)

# see here: http://mnras.oxfordjournals.org/content/322/2/419.full
OM(z) = Ω0 / (Ω0 + ΩΛ./(z+1.0).^3)
OL(z) = ΩΛ / (Ω0.*(z+1.0).^3 + ΩΛ)

function unnormalizedDa(z)
    o0 = OMa(z)
    ol = OLa(z)
    up = 5o0
    down = o0^(4/7) - ol + (1+o0/2)*(1+ol/70)
    a*up/2down
end
D(z) = unnormalizedDa(z) / unnormalizedDa(1.0)
D2(z) = -3/7*D(z)^2

F(z) = OMa(z)^(5/9)
F2(z) = 2*OMa(z)^(6/11)

# mean/critical density
ρc(z) = ρ0(z)./OM(z)

A(z) = 1.0./(z+1.0)

ρ0(z) = Ω0.*RHO_CRIT_TODAY.*(1+z).^3 # mean rho
function RXb(z, m, X)
    vX = m./ρ0(z)./X
    (vX.*3./4./π).^(1/3) 
end
R200b(z, m) = RXb(z, m, 200.0)
function RXc(z, m, X)
    vX = m./ρc(z)./X
    (vX.*3./4./π).^(1/3) 
end
R200c(z, m) = RXc(z, m, 200.0)

RVIR(z, m) = RXc(z, m, Δcvir(z))

# line of sight comoving distance
# (see: https://arxiv.org/pdf/astro-ph/9905116v4.pdf)
E(z) = sqrt(Ω0.*(1+z).^3+ΩΛ)
const DH = C/H0
Dc(z) = DH.*quadgk(z->1.0./E(z), 0, z)[1]

# angular distance, in arcseconds
rad_to_as(Θ) = Θ*206265
as_to_rad(Θ) = Θ/206265
Da_rad(z, dθ) = Dc(z)./(1+z).*dθ
Da_as(z, dθ) = Da_rad(z, as_to_rad(dθ))

# 

function Δcvir(z)
    x = OM(z)-1.0
    18.*π.*π+82.*x-39.*x.*x
end

