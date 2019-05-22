module Experiments

include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")

using .PartitionsGen
using .QuadraticPartitions
using SymPy

function O_plusplus_vis(D)
  # first index is a, second is b
  A = zeros(Int,12,12)
  for i = 1:12
    A[:,i] = length.(quad_partitions.(1:12,i,D,true))
  end
  A'
  # heatmap(A') 
  #   with Plots.jl
end

function O_plus_vis(D)
  size = 6
  A = zeros(size*2 + 1, size*2 + 1) # the 1 is for zero
  for i = -6:6
    A[:,i+size+1] = length.(quad_partitions.(-6:6,i,D,false))
  end
  A'
end
function symbolic_gen_O_plusplus(k)
  n = Vector{SymPy.Sym}(k)
  x = Sym("x")
  expr = Sym(1)
  for i in 1:k
    ni = Sym("n$i")
    expr *= (1 - x^(ni + Sym(i)*√Sym(2)))
  end
  expr
end

end#module
