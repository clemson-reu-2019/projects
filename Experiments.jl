module Experiments

include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")

using .PartitionsGen
using .QuadraticPartitions
using SymPy
using AbstractAlgebra

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

function abstract_symbolic_gen(k,D)
  R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
  expr = 1
  for i = 1:k 
    expr *= (1 - x^i)
  end 
  for i in 1:k
	for nᵢ = convert(Int,ceil(i*√D)):k
		expr *= (1 - (x^(nᵢ)) * y^i )
	end
  end
  expr
end

function symbolic_gen_O_plusplus(k,D)
  symbolic_gen_euler(k)*symbolic_gen_sqareroot(k,k,D)
end

function symbolic_gen_sqareroot(k,j,D)
  x = Sym("x")
  expr = Sym(1)
  for i in 1:k
	for nᵢ = convert(Int,ceil(i*√D)):j
      expr *= (1 - x^(nᵢ + Sym(i)*√Sym(2)))
	end
  end
  expr
end

function symbolic_gen_euler(k)
  x = Sym("x")
  expr = Sym(1)
  for i = 1:k
	expr *= (1 - x^i)
  end
  expr
end

function parts(x)
  (a,b) = (x.args[1],x.args[2])
  b = b/√Sym(3)
  (a,b)
end

function get_distinguished_element(a,b,D=3)
  D != 3 && return -1
  unit = 2 - √Sym(3)
  x₀ = Sym(a) + Sym(b)*√Sym(D)
  (a_best, b_best) = parts(x₀)
  i = 0
  xᵢ = x₀
  while true
	i += 1
	xᵢ₊₁ = simplify(xᵢ * expand(unit^i))
	(a, b) = parts(xᵢ₊₁)
	if is_wholly_positive(a,b,D) && 
		abs(a^2 + b^2) < abs(a_best^2 + b_best^2) # euclidean norm
      (a_best, b_best) = (a,b)
	  xᵢ = xᵢ₊₁
  	else
	  break
  	end
  end
  (a_best,b_best)
end



end#module
