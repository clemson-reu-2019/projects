module QuadraticPartitions

using Combinatorics
using Partitions

"""
Returns true if x + y√D is wholly positive
"""
function is_wholly_positive(x,y,D)
  0 < x && 3y^2 < x^2
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

order(parti) = sortrows(parti, by=x->(x[1],x[2]))

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
function partitions_allpositive(a,b,D)
  ps = Array{Matrix{Int}, 1}()
  !is_wholly_positive(a,b,D) && return ps
  for pₐ in partitions(a)
	  lₐ = length(pₐ) # - b *floor(sqrt(D))
    for pᵦ in partitions_lessterms(b,lₐ)
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
export partitions_allpositive

"""
Generate all partitions of a + b√D in O₊₊ which satisfy
the predicate pred
"""
function partitions_allpositive(a,b,D,pred::Function)
  filter(pred,partitions_allpositive(a,b,D))
end

end#module
