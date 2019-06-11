module PartitionsGen

using Memoize

function no_parts_exceeding(partition,n)
  for part in partition
    if n < part
      return false
    end
  end
  true
end
export no_parts_exceeding

# Generating Partitions
# 
# The following algorithms are used to generate
# partitions themselves. They all perform computations
# and are memoized with the @memoize macro

#"""
#Generate the partitions of n which have no elements lower than l
#https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
#
#Note: this seems to get bad at around 65 (with memoization)
# because of too many allocations.
#"""
@memoize function partitions_lowest(n,l)
  ps = Array{Vector{Int}, 1}()
  push!(ps,[n])
  for i in l:(n ÷ 2) #integer division
    for P in partitions_lowest(n-i,i)
      push!(ps,[i,P...])
    end
  end
  ps
end

#"""
#Generate partitions of n which have less terms than l
# This algorithm is rule_asc from 
# http://jeromekelleher.net/category/combinatorics.html
#"""
@memoize function rule_asc(n,l)
  n == 1 && return [[1]]
  ps = Array{Vector{Int}, 1}()
  A = zeros(Int, n)
  k = 2
  A[2] = n

  while k != 1
    x = A[k-1] + 1
    y = A[k] - 1
    k -= 1
    while x <= y && k < l
      A[k] = x
      y -= x
      k += 1
    end
    A[k] = x + y
    push!(ps, deepcopy(A[1:k]))
  end
  ps
end

# Partition Interface
# 
# The following functions provide an interface to 
# easily get partitions. They don't necessarily 
# do computation, and they aren't optimized for
# any sort of efficieny

"""
    partitions_of(n)

Generate all partitions of n
"""
function partitions_of(n)
  rule_asc(n,n)
end
export partitions_of


"""
    partitions_of(n,nParts,bound=n)

Generate all partitions of n which have nParts parts with all
parts less than bound
"""
function partitions_of(n,nParts,bound=n)
  if n <= bound
    ps = partitions_of(n)
  else
    ps = partitions_leq(n,bound)
  end
  filter(p -> nParts == length(p), ps)
end

"""
Generate all partitions of n which satisfy predicate pred
"""
function partitions_of(n,pred::Function)
  filter(pred,partitions_of(n))
end

# a more efficient algorithm for this might exist
"""
Generate all partitions of n which are less than or equal to l
"""
partitions_leq(n,l) = partitions_of(n, p -> no_parts_exceeding(p,l))
export partitions_leq

"""
Generate all partitions of n which have m or less terms
"""
partitions_lessterms(n,m) = rule_asc(n,m)
export partitions_lessterms

"""
Generate all partitions of n which have m or less terms and which
have each part less than or equal to l
"""
partitions_lessterms_leq(n,m,l) = filter(p -> no_parts_exceeding(p,l), partitions_lessterms(n,m))
export partitions_lessterms_leq

end#module
