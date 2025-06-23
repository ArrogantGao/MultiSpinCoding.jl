module MultiSpinCoding

using Graphs, LinearAlgebra
using BitBasis

export SpinGlass, Scheduler
export random_J, random_state, bit_J
export sa!, cal_energies

include("utils.jl")
include("types.jl")
include("sa.jl")

end
