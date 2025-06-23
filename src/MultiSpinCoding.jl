module MultiSpinCoding

using Random
using Graphs, LinearAlgebra
using BitBasis

export SpinGlass, Scheduler
export random_J, random_state, bit_J
export load_instance
export sa!, cal_energies

include("utils.jl")
include("types.jl")
include("sa.jl")

end
