include("Partitions.jl")
include("QuadraticPartitions.jl")

module Experiments

using QuadraticPartitions
using SymPy

function O_plusplus_vis(D)
  # first index is a, second is b
  A = zeros(Int,12,12)
  for i = 1:12
    A[:,i] = length.(partitions_allpositive.(1:12,i,D))
  end
  A'
  # heatmap(A') 
  #   with Plots.jl
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
