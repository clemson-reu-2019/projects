module QuadraticPartitions

using Combinatorics
using Memoize
using DelimitedFiles
using ..PartitionsGen

"""
Returns true if x + y√D is wholly positive
"""
function is_wholly_positive(x,y,D)
  0 < x && D*y^2 < x^2
end
export is_wholly_positive

sqrt2_80_dat = convert(Matrix{Int}, readdlm("80.dat"))
sqrt2_plus_50_dat = convert(Matrix{Int}, readdlm("50plus.dat"))

"""
Returns the euler coeffieicent corresponding
to a certain wholly positive number.
"""
function euler_coef(a,b,D,allpositive)
  D != 2 && (print("2 only supported!"); return 0) # not implemented
  !is_wholly_positive(a,b,D) && return 0

  if allpositive 
    sqrt2_80_dat[a+1,b+1]
  else
    sqrt2_plus_50_dat[a+1,convert(Int, abs(b))+1]
  end
end

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

If allpositive is true, then search for partitions
in O₊₊

Otherwise, if allpositive is false, 
find partitions in O₊ 
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
function quad_partitions(a,b,D,pred::Function,allpositive=true)
  filter(pred,quad_partitions(a,b,D,allpositive))
end

"""
Generate all of the generalized partitions (that is, the
partitions with negative numbers included) which have
parts all less than maxpart, and no more than naxNumParts
parts.
"""
function generalized_partitions(b,maxNumParts,maxpart)
  maxpart == 1 && return Array{Vector{Int}, 1}([[b]])
  #maxpart*maxNumParts < b && return Array{Vector{Int}, 1}([[]])

  # maybe this belongs here for type safety?
  b = convert(Int, abs(b))

  ps = Array{Vector{Int}, 1}()
  N = maximum((maxpart*maxNumParts,b))
  plists = Array{Array{Vector{Int}, 1}}(undef, N)
  for i = 1:N
    plists[i] = partitions_lessterms_leq(i,maxNumParts,maxpart)
  end
  
  if b != 0 
	0 < b && append!(ps,partitions_lessterms(b,maxNumParts))
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

#"""
#Generate the partition number in a brute-force manner
#"""
@memoize function partition_number_brute(a,b,D,allpositive=false)
  (a == 0 && b == 0) && return 1
  length(quad_partitions(a,b,D,allpositive))
end

"""
Give the biggest value of a and b where the 
brute force algorithm is performant, tested manually.
"""
function p_num_brute_force_bound(D,allpositive=false)
  return 3
  (D == 2 && allpositive) && return 10
  (D == 2 && !allpositive) && return 7
  (D == 3 && allpositive) && return 9
  (D == 3 && !allpositive) && return 7
  (D == 5 && allpositive) && return 11
  (D == 5 && !allpositive) && return 7
  5 # a decent bet that this will be fast
end

#"""
#Calculate the partition number p(n) of n = a + b√D
#using the recursive euler expansion
#
#if allpositive is true, calculates p₊(n)
#"""
@memoize function partition_number(a,b,D,allpositive=false)
  #println("Starting $a + $b√2")
  0 < b && return partition_number(a,-b,D,allpositive)

  if a == 0 && b == 0
    return 1
  end
  # b is necessarily less than a
  #if a <= p_num_brute_force_bound(D,allpositive) 
  #  print("$a found: brute\n")
	#  return partition_number_brute(a,b,D,allpositive)
  #end

  #p = partition_number(a-1,b,D,allpositive)
  #print("$a-1:$b, $p\n")

  p = 0
  for i = 0:a-1
    if allpositive
      start = 0
      last = b
    else
      bound = ceil(Int, i*√D)
      start = -bound#- ceil(Int, i*√D)
      last = +bound
    end
    #println("START: $start for $i")

    for j = start:last
      #println("$i,$j")
      eul = euler_coef(a-i,b-j,D,allpositive)
      #println("  ...found $eul")
      eul == 0 && continue
      #println("nonzero coef found!")
      #pn = partition_number(i,j,D,allpositive)
      p -= eul * partition_number(i,j,D,allpositive)
      #print("running total $a + $b√2: $p\n")
      #println("term added to $a + $b√2 from e_$(a-i),$(b-j): - $eul * $pn")
    end
  end
  #println("done computing for $a + $b√2: got $p")
  p
end
export partition_number

"""
Generate a 2-D grid of values of p(n) of size N
by brute-forcing it.

If allpositive is true, generate p₊(n) instead
"""
function partitions_grid_brute(N,D,allpositive=false)
  A = zeros(Int,N+1,N+1)
  for i = 0:N
    A[:,i+1] = length.(quad_partitions.(0:N,i,D,allpositive))
  end
  A'
end

"""
Generate a 2-D grid of values of p(n) of size N
by using the recursive Euler algorithm.

If allpositive is true, generate p₊(n) instead
"""
function partitions_grid(N,D,allpositive=false)
  A = zeros(Int,N+1,N+1)
  for i = 0:N
    for j = 0:N
      if is_wholly_positive(i,j,D)
        A[i+1,j+1] = partition_number(i,j,D,allpositive)
      end
    end
  end
  A'
end
export partitions_grid

end#module
