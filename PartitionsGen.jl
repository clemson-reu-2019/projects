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

#"""
#Generate the partitions of n which have no elements lower than l
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

"""
    partitions_of(n)

Generate all partitions of n
"""
function partitions_of(n)
  partitions_lowest(n,1)
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

"""
Generate all partitions of n which are less than or equal to m
"""
partitions_leq(n,m) = partitions_of(n, p -> no_parts_exceeding(p,m))
export partitions_leq

"""
Generate all partitions of n which have m or less terms
"""
partitions_lessterms(n,m) = partitions_of(n,p -> length(p) <= m)
export partitions_lessterms

# TODO???: write an efficient method for this or the previous method
#          and memoize
"""
Generate all partitions of n which have m or less terms and which
have each part less than or equal to l
"""
partitions_lessterms_leq(n,m,l) = partitions_of(n,p -> no_parts_exceeding(p,l) && length(p) <= m)
export partitions_lessterms_leq

end#module
