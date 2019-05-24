module QuadraticPartitions

using Combinatorics
using ..PartitionsGen

"""
Returns true if x + y√D is wholly positive
"""
function is_wholly_positive(x,y,D)
  0 < x && D*y^2 < x^2
end
export is_wholly_positive

"""
Takes a partition p of some integer a, and 
a partition q of some integer b, and an injective
function (in the form of a permutation) from the
parts of q to the parts of p which defines a
quadratic partitions.

Returns a 2-D array where the first row is p and the
second row is q, but ordered by the assignment of p.

Example:
p = [1,1,1,1,1] a partition of 5
q = [1,2,3] a partition of 6
perm = [5,1,2]
Output:
[ 1 2
  1 3
  1 0
  1 0
  1 1 ]
"""
function quad_partition(p,q,perm)
  sqrow = zeros(eltype(q), length(p))
  sqrow[perm] = q
  [p sqrow]
end

order(parti) = sortslices(parti, dims=1, by=x->(x[1],x[2]))

"""
Takes a partition in 2-D array form (see quad_partition)
and checks if every entry is wholly positive
"""
function all_wholly_positive(parti,D)
  for i in 1:size(parti,1)
    if parti[i,1] - abs(parti[i,2])*√D < 0
      return false
    end
  end
  return true
end

"""
Generate all partitions of a + b√D in Q(√D)
which reside in the set O₊₊ 
"""
function quad_partitions(a,b,D,allpositive=false)
  bound = a - b*floor(Int,√D)
  if allpositive
    bparts_func = l -> partitions_lessterms(b,minimum((l,bound)))
  else
    bparts_func = l -> generalized_partitions(b,bound,minimum((l,bound)))
  end
  ps = Array{Matrix{Int}, 1}()

  !is_wholly_positive(a,b,D) && return ps

  # manually add in partitions with no generalized term for 0
  if !allpositive && b == 0
    for p in partitions(a)
      push!(ps,quad_partition(p,zeros(Int,0),[]))
    end
  end

  for pₐ in partitions(a)
	  lₐ = length(pₐ)
    for pᵦ in bparts_func(lₐ)
	  pᵦ == [] && continue
      # skip a few examples that would really waste time
      if !is_wholly_positive(maximum(pₐ),minimum(pᵦ),D)
        continue
      end

      for perm in permutations([1:lₐ;],length(pᵦ))
        p = order(quad_partition(pₐ,pᵦ,perm))
        if all_wholly_positive(p,D)
          push!(ps,p)
        end
      end
    end
  end
  unique(ps)
end
export quad_partitions

"""
Generate all partitions of a + b√D in O₊₊ which satisfy
the predicate pred
"""
function partitions_allpositive(a,b,D,pred::Function,allpositive=true)
  filter(pred,partitions_allpositive(a,b,D,allpositive))
end

"""
Generate all of the generalized partitions (that is, the
partitions with negative numbers included) which have
parts all less than maxpart, and no more than naxNumParts
parts.
"""
function generalized_partitions(b,maxNumParts,maxpart)
  maxpart == 1 && return Array{Vector{Int}, 1}([[b]])
  maxpart*maxNumParts < b && return Array{Vector{Int}, 1}([[]])

  ps = Array{Vector{Int}, 1}()
  plists = Array{Array{Vector{Int}, 1}}(undef, maximum((maxpart*maxNumParts,b)))
  for i = 1:maxpart*maxNumParts
    plists[i] = partitions_lessterms_leq(i,maxNumParts,maxpart)
  end
  
  if b != 0 
    0 < b && append!(ps,plists[b])
    b < 0 && append!(ps,-1 .* reverse.(plists[-b]))
  end
  for i = abs(b)+1:maxpart*maxNumParts
    for pᵢ in plists[i]
      l = length(pᵢ)
      compat = filter(p -> length(p) <= maxNumParts-l, plists[i-abs(b)])
      for pm in compat
        gen_p = [-reverse(pm); pᵢ]
        b < 0 && (gen_p = -reverse(gen_p))
        push!(ps,gen_p)
      end
    end
  end
  ps
end
export generalized_partitions

end#module
