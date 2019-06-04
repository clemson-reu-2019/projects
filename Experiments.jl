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
using LsqFit

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
  evalfile(filename * "jl")
end

"""
Process the data from macaulay2 into an array which can be
used to compute partition numbers
"""
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

function process_new_coefficients(filename,D,allpositive)
  tuples = process_macaulay2_polynomial(filename)
  coefs = process_macaulay2_data(tuples)
  QuadraticPartitions.incorporate_coefficients(coefs,D,allpositive)
end

function generate_norm_array(N,D)
  M = ceil(Int, N/√D)
  NORMS = zeros(Int, N+1,M+1)
  for a = 0:N
    for b = 0:M
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

quadratic_mult((a,b),(c,d),D) = (a*c + D*b*d, b*c + a*d)

function quadratic_power((a,b),N,D)
  prod = (1,0)
  for i = 1:N
    prod = quadratic_mult(prod,(a,b),D)
  end
  prod
end

function modelexpsin(ydata,c₀)
  @. model(x,c) = c[1]*exp(c[2]*x)*sin(c[3]*√(x+ c[4]))
  curve_fit(model, 1:length(ydata),ydata,c₀)
end

function linreg(x,y)
  @. model(x,c) = c[1]*x + c[2]
  curve_fit(model,x,y,[1.0,0.0])
end
#  -0.9734146840101862 
#   0.18881479886532654
#   9.13793675064356   
# 104.04780167727091

function modelexpsqrt(ydata,c₀)
  @. model(x,c) = c[1]*exp(c[2]*√x)
  curve_fit(model,1:length(ydata),ydata,c₀)
end

function modelsqrt(ydata,c₀)
  @. model(x,c) = c[1]*√x + c[2]
  curve_fit(model,1:length(ydata),ydata,c₀)
end

function highestBFor(a,D)
  floor(Int, a / √D)
end

# goes one member beyond limitA right now
function pellClassFor(a,b,D=2,unit=(3,2),limitA=100)
  class = [(a,b)]#Array{Tuple{typeof(a),typeof(b)}}(undef,0)
  current = (a,b)
  while current[1] <= limitA
    current = quadratic_mult(current,unit,D)
    class = [class; current]
  end

  conjunit = (unit[1], -unit[2])
  current = (a,b)
  while current[1] <= limitA
    current = quadratic_mult(current,conjunit,D)
    class = [class; current]
  end
  class
end

samepellclass(cl1,cl2) = (Set(cl1) == Set(cl2))

function addconjugates(pellclass)
  class = Set(pellclass)
  for (a,b) in pellclass
    push!(class,(a,-b))
  end
  collect(class)
end

const PellClass{T} = Array{Tuple{T,T}}

function findpellclasses(N,D=2,unit=(3,2))
  norms = generate_norm_array(N,D)
  classes = Dict{Int,Array{PellClass{Int}}}()
  boundOnFundSol = floor(Int,N/√(unit[1] + unit[2]*√D))

  pushclass(c,i) = haskey(classes,i) ? push!(classes[i], c) : classes[i] = [c]

  for i = 1:boundOnFundSol
    matches = [(ind[1],ind[2]) for ind in findall(norms .== i)]
    matches = map(m -> m .- 1, matches)

    while 0 < length(matches)
      firstclass = pellClassFor(matches[1]...,D,unit,N)
      pushclass(firstclass, i)

      accumSet = addconjugates(firstclass)
      matches = filter(m -> !(m in accumSet), matches)
    end
  end
  classes
end

function comparepellclasses(N,D=2,unit=(3,2),allclasses=nothing)
  allclasses == nothing && (allclasses = Experiments.findpellclasses(N))
  println("Evaluating all pell classes for counterexample")
  outliers = filter(pair -> 1 < length(pair[2]),allclasses)
  for (norm,classes) in outliers
    # minimum a
    part_nums = zeros(Int,length(classes))
    for i = 1:length(classes)
      minind = argmin(map(t -> t[1],classes[i]))
      part_nums[i] = partition_number(classes[i][minind]...,D)
    end
    #println("$part_nums")
    if !allunique(part_nums)
      println("Counterexample found! $norm")
    end
  end
  outliers
end

function pellclasses_to_thetafn(pellclasses)
  θ = zeros(Int, maximum(keys(pellclasses)))
  for (norm,classes) in pellclasses
    θ[norm] = length(classes)
  end
  θ
end
    
    


end#module
