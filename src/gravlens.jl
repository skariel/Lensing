using PyPlot
using ODE
using Dierckx
using Roots

setprecision(128);

include("cosmology.jl")
include("rays.jl")
include("mass.jl")
include("harmonic.jl")
include("numeric.jl");
include("thin.jl")
include("time.jl")
include("precision.jl")
include("defaultmodels.jl")
include("fermat.jl")
