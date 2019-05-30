module Experiments

include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")
#include("./EulerCoefficients.jl")

using .PartitionsGen
using .QuadraticPartitions
#using .EulerCoefficients
using SymPy
using AbstractAlgebra
using DelimitedFiles
using Plots
using PlotUtils

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

function O_plus_vis(D,size=6)
  A = zeros(Int, size*2 + 1, size*2 + 1) # the 1 is for zero
  for i = -size:size
    A[:,i+size+1] = length.(quad_partitions.(-size:size,i,D,false))
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

"""
take the output from expander.m2 and process it into a form
that can be copied and pasted into julia
"""
function process_macaulay2_polynomial(filename)
  replace_at(s,i,new) = s[1:i-1] * new * s[i+1:end]
  
  m2out = read(filename, String)
  # remove all ---- lines
  lines = split(m2out, "\n")
  lines = filter(line -> !occursin("-----", line), lines)
  nodashes = join(lines)
  # remove all whitespace
  nowhsp = filter(x -> !isspace(x), nodashes)
  # replace first { with [
  i = findfirst("{",nowhsp).start
  opensquare = replace_at(nowhsp,i,"[")
  # replace last } with ]
  i = findlast("}",opensquare).start
  clsdsquare = replace_at(opensquare,i,"]")
  # replace all { with ( and } with )
  noopenbrak = replace(clsdsquare, "{" => "(")
  noclsdbrak = replace(noopenbrak, "}" => ")")
  open(filename * "jl", "w") do f
	write(f, noclsdbrak)
  end
end

function process_macaulay2_data(tuples)
  sparse = zeros(Int, length(tuples), 3)
  for i = 1:length(tuples)
    ((a,b),c) = tuples[i]
    sparse[i,:] = [a,b,c]
  end
  maxA = maximum(sparse[:,1])
  maxB = maximum(sparse[:,2])
  recurGrid = zeros(Int, maxA+1, maxB+1)
  for i = 1:length(tuples)
    ((a,b),c) = tuples[i]
	0 <= b && (recurGrid[a+1,b+1] = c)
  end
  recurGrid
end

function generate_norm_array(N,D)
  NORMS = zeros(Int, N+1,N+1)
  for a = 0:N
    for b = 0:N
      if Experiments.QuadraticPartitions.is_wholly_positive(a,b,D)
        NORMS[a+1,b+1] = a^2 - D*b^2
      end
    end
  end
  NORMS
end

function nonzerocoefs(M,level=-10000)
  A = deepcopy(M)
  A[A .== 0] .= level
  heatmap(A, yflip=true, c=cgrad([:red,:blue],[0.1,1]))
  gui()
end

function plotrow(M,row)
  G = deepcopy(M[row,:])
  plot(G)
  gui()
end

function logwithzeros(M,level=0.1)
  N = deepcopy(M)
  N[N .== 0] .= level
  log.(N)
end

getpointsequal(M,n) = foldl((A,x) -> [A; x[1] x[2]], findall(M .== n), init=zeros(Int,0,2))

function timesunit(a,b,n,U =3 + 2*√Sym(2))
  (a + b*√Sym(2)) * (U^n)
end

end#module
