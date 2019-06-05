module PellClasses

include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")
using .QuadraticPartitions

const PellClass{T} = Array{Tuple{T,T}}

qmult((a,b),(c,d),D) = (a*c + D*b*d, b*c + a*d)

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

samepellclass(cl1,cl2) = (Set(cl1) == Set(cl2))

function addconjugates(pellclass)
  class = Set(pellclass)
  for (a,b) in pellclass
    push!(class,(a,-b))
  end
  collect(class)
end


function findpellclasses(N,D=2,unit=(3,2),ignoreconj=true)
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

      accumSet = ignoreconj ? addconjugates(firstclass) : firstclass
      matches = filter(m -> !(m in accumSet), matches)
    end
  end
  classes
end

function pellclasses_to_thetafn(pellclasses)
  θ = zeros(Int, maximum(keys(pellclasses)))
  for (norm,classes) in pellclasses
    θ[norm] = length(classes)
  end
  θ
end

function pellthetafn(N,D=2,unit=(3,2),ignoreconj=true)
  pellclasses_to_thetafn(findpellclasses(N,D,unit,ignoreconj))
end

end#module
