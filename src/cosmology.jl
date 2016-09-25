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

# using Dark-Sky halos box 1600 Mpc/h:
const SIDE_LEN = 1.6e6
const RHO_CRIT = 3*H0^2/8/π/G
const BOX_VOLUME = SIDE_LEN^3
const MASS_IN_BOX = RHO_CRIT*BOX_VOLUME*Ω0
const BOX_N = 128
const SMTH = 3000.0
const MEAN_SEP_FULL_REALIZATION = 3010.314648627409
const PARTICLE_NUM_FULL_REALIZATION = 150149632
const MEAN_SEP_SMALL_REALIZATION = 3806.0162249593727
const PARTICLE_NUM_SMALL_REALIZATION = 74293025
const DX = SIDE_LEN / BOX_N
const xl = linspace(0,SIDE_LEN,BOX_N)

# Hubble as function of scale factor a
Ha(a) = H0*sqrt(Ω0/a/a/a+ΩΛ)

# drift factor
FD(a1,a2) = quadgk(a->1./Ha(a)/a/a/a,a1,a2)[1]

# kick factor
FK(a1,a2) = quadgk(a->1./Ha(a)/a/a,a1,a2)[1]

# for convenience...
const FAC1 = 0.5*H0*H0*Ω0

# see here: http://mnras.oxfordjournals.org/content/322/2/419.full
OM(a) = Ω0 / (Ω0 + ΩΛ*a*a*a )
OL(a) = ΩΛ / (Ω0/a/a/a + ΩΛ )

function unnormalizedD(a)
    o0 = OM(a)
    ol = OL(a)
    up = 5o0
    down = o0^(4/7) - ol + (1+o0/2)*(1+ol/70)
    a*up/2down
end
D(a) = unnormalizedD(a) / unnormalizedD(1.0)
D2(a) = -3/7*D(a)^2

F(a) = OM(a)^(5/9)
F2(a) = 2*OM(a)^(6/11)

# mean density
Z(a) = 1./a-1
ρ0(a) = OM(a).*RHO_CRIT
function R200(a, m)
    v200 = m./ρ0(a)./200
    (v200.*3./4./π).^(1/3) 
end

# line of sight comoving distance
# (see: https://arxiv.org/pdf/astro-ph/9905116v4.pdf)
E(z) = sqrt(Ω0.*(1+z).^3+ΩΛ)
const DH = C/H0
Dc(z) = DH.*quadgk(z->1.0./E(z), 0, z)[1]

# angular distance, in arcseconds
rad_to_as(Θ) = Θ*206265
as_to_rad(Θ) = Θ/206265
Da_radians(z, dθ) = Dc(z)./(1+z).*dθ
Da_arcseconds(z, dθ) = Da_radians(z, as_to_rad(dθ) )


