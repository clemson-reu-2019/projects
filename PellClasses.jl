module PellClasses

include("./EulerCoefficients.jl")
using .EulerCoefficients

using Primes

const PellClass{T} = Array{Tuple{T,T}}

qmult((a,b),(c,d),D) = (a*c + D*b*d, b*c + a*d)

"""
Generate an 2D array which contains all of the norms
for wholly positive numbers a + b√d which have 
a ≦ N
"""
function generate_norm_array(N,D)
  M = ceil(Int, N/√D)
  NORMS = zeros(Int, N+1,M+1)
  for a = 0:N
    for b = 0:M
      if is_wholly_positive(a,b,D)
        NORMS[a+1,b+1] = a^2 - D*b^2
      end
    end
  end
  NORMS
end
export generate_norm_array

# goes one member beyond limitA right now
function pellClassFor(a,b,D=2,unit=(3,2),limitA=100)
  class = [(a,b)]#Array{Tuple{typeof(a),typeof(b)}}(undef,0)
  current = (a,b)
  while current[1] <= limitA
    current = qmult(current,unit,D)
    class = [class; current]
  end

  conjunit = (unit[1], -unit[2])
  current = (a,b)
  while current[1] <= limitA
    current = qmult(current,conjunit,D)
    class = [class; current]
  end
  class
end
export pellClassFor

samepellclass(cl1,cl2) = (Set(cl1) == Set(cl2))

function minelement(class)
  minind = argmin(map(t -> t[1],class))
  class[minind]
end
export minelement

function addconjugates(pellclass)
  class = Set(pellclass)
  for (a,b) in pellclass
    push!(class,(a,-b))
  end
  collect(class)
end

"""
Find all pell classes which have a primitive
solution which is bounded by A
"""
function findpellclasses(A,D=2,unit=(3,2),ignoreconj=true)
  norms = generate_norm_array(A,D)
  classes = Dict{Int,Array{PellClass{Int}}}()
  #boundOnFundSol = floor(Int,N/√(unit[1] + unit[2]*√D))
  boundOnFundSol = floor(Int,A^2/(unit[1] + unit[2]*√D))

  pushclass(c,i) = haskey(classes,i) ? push!(classes[i], c) : classes[i] = [c]

  for i = 1:boundOnFundSol
    matches = [(ind[1],ind[2]) for ind in findall(norms .== i)]
    matches = map(m -> m .- 1, matches)

    while 0 < length(matches)
      firstclass = pellClassFor(matches[1]...,D,unit,A)
      pushclass(firstclass, i)

      if !ignoreconj
        conj = (matches[1][1], -matches[1][2])
        if !(conj in firstclass)
          conjclass = pellClassFor(conj...,D,unit,A)
          pushclass(conjclass, i)
        end
      end

      accumSet = addconjugates(firstclass)
      matches = filter(m -> !(m in accumSet), matches)
    end
  end
  classes
end
export findpellclasses

"""
Compute the pell theta function 
from a the output of findpellclasses()
"""
function pellthetafn(pellclasses)
  θ = zeros(Int, maximum(keys(pellclasses)))
  for (norm,classes) in pellclasses
    θ[norm] = length(classes)
  end
  θ
end
export pellthetafn

function pellthetafn(N,D,unit,ignoreconj)
  pellclasses_to_thetafn(findpellclasses(N,D,unit,ignoreconj))
end

"""
Decide whether a prime is a solvable prime

Does not consider the case when p is not prime. 
Be sure that p is a prime!!!
"""
function is_solvable(p,D)
  D == 2 && return (p % 8 == 1 || p % 8 == 7)
  D == 5 && return (p % 10 == 1 || p % 10 == 9)
  D == 10 && return (is_solvable(p,2) && is_solvable(p,5))
  false
end

function is_symmetric(n,D)
  D == 2 && return n == 2
  false
end

function theta_coef(n,D)
  !(D == 2 || D == 5) && throw(ArgumentError("not supported"))
  n == 1 && return 1
  factzn = factor(Dict,n)
  SlvP = Dict(1 => 0)

  for (pᵢ,eᵢ) in factzn
    if is_solvable(pᵢ,D)
      SlvP[pᵢ] = eᵢ
    else
      if !is_symmetric(pᵢ,D) && (eᵢ % 2 == 1)
        return 0
      end
    end
  end

  prod(collect(values(SlvP)) .+ 1)
end
export theta_coef

end#module
