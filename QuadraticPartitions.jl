module QuadraticPartitions

include("./EulerCoefficients.jl")
include("./TotallyLess.jl")
include("./OnePart.jl")
include("./PellClasses.jl")

import Base.Threads.@spawn
import Base.Threads.@sync

using Combinatorics
using Memoize
using Primes
using ..PartitionsGen
using .EulerCoefficients
using .TotallyLess
using .OnePart
using .PellClasses

# BRUTE FORCE ALGORITHM

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
  nPartsBound = floor(Int,√(a^2 - D*b^2))
  if allpositive
    bparts_func = l -> partitions_lessterms(b,minimum((l,bound,nPartsBound)))
  else
    bparts_func = l -> generalized_partitions(b,minimum((l,bound,nPartsBound)),bound)
  end
  ps = Array{Matrix{Int}, 1}()

  !is_wholly_positive(a,b,D) && return ps

  # manually add in partitions with no generalized term for 0
  if !allpositive && b == 0
    for p in partitions(a)
      push!(ps,quad_partition(p,zeros(Int,0),[]))
    end
  end

  for pₐ in partitions_lessterms(a, nPartsBound)
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

#"""
#Give the biggest value of a and b where the
#brute force algorithm is performant, tested manually.
#"""
#function p_num_brute_force_bound(D,allpositive=false)
#  (D == 2 && allpositive) && return 10
#  (D == 2 && !allpositive) && return 7
#  (D == 3 && allpositive) && return 9
#  (D == 3 && !allpositive) && return 7
#  (D == 5 && allpositive) && return 11
#  (D == 5 && !allpositive) && return 7
#  5 # a decent bet that this will be fast
#end

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

# BRUTE FORCE ALGORITHM USING DECOMPOSITIONS

# there is a bug if b \neq 0
function quad_partitions_decomp(a,b,d,allpositive=false)
  allpositive && throw(ArgumentError("Ok++ not supported!"))
  ps = Set([])

  for decomp in OnePart.QuickDecompositions(a,b,d)
    if all(iszero, decomp .|> x -> x[2])
      # decomposition is integral, thus b = 0
      for partition in partitions_of(a)
        push!(ps, [(x, 0) for x in partition])
      end
    else
      for setpart in partitions(decomp)
        parti = map(set -> foldl((x,y) -> x .+ y, set), setpart)
        sorted = sort(sort(parti, by = x -> x[2]), by = x -> x[1])
        push!(ps, sorted)
      end
    end
  end

  collect(ps)
end

# RECURSIVE ALGORITHM USING EULER PRODUCT EXPANSION

#"""
#Calculate the partition number p(n) of n = a + b√D
#using the recursive euler expansion
#
#if allpositive is true, calculates p₊(n)
#"""
@memoize function partition_number_euler(a,b,D,allpositive=false)
  #println("Starting $a + $b√2")
  b < 0 && return partition_number_euler(a,-b,D,allpositive)

  if a == 0 && b == 0
    return 1
  end

  # So it turns out that the brute force is still much slower
  # than this algorithm even at very small values...
  # honestly not that surprising, but for that reaason
  # the fallback is not necessary:
  #
  # b is necessarily less than a
  #if a <= p_num_brute_force_bound(D,allpositive)
  #  #print("$a found: brute\n")
  #  return partition_number_brute(a,b,D,allpositive)
  #end

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
      #pn = partition_number_euler(i,j,D,allpositive)
      p -= eul * partition_number_euler(i,j,D,allpositive)
      #print("running total $a + $b√2: $p\n")
      #println("term added to $a + $b√2 from e_$(a-i),$(b-j): - $eul * $pn")
    end
  end
  #println("done computing for $a + $b√2: got $p")
  p
end
export partition_number_euler


"""
Generate a 2-D grid of values of p(n) of size N
by using the recursive Euler algorithm.

If allpositive is true, generate p₊(n) instead
"""
function partitions_grid(t::Type,N,D,allpositive=false,alg=nothing,data=nothing)
  alg == nothing && (alg = partition_number)
  maxB = floor(Int, N / √D)
  A = zeros(t,N+1,maxB+1)
  for i = zero(t):N
    @sync for j = zero(t):maxB
      if is_wholly_positive(i,j,D)
        if data == nothing
          A[i+1,j+1] = alg(i,j,D,allpositive)
        else
          @spawn A[i+1,j+1] = alg(i,j,D,allpositive,data)[2]
        end
      end
    end
  end
  A'
end
export partitions_grid

partitions_grid(N,D,allpositive=false,alg=nothing) = partitions_grid(Int,N,D,allpositive,alg)

# RECURSIVE ALGORITHM USING FINITE PRODUCT EXPANSION

#"""
#Calculate the integer partition number pᵣ(n)
#which gives the number of partitions of n
#with less than or equal to r parts
#"""
@memoize function partition_number_leq(t::Type,n,r)
  r == 1 && return one(t)
  n == 0 && return one(t)
  sum(partition_number_leq.(t,n:-r:0,r-1))
end
export partition_number_leq

function partition_number_leq_array(t::Type,n,r)
  partitiongrid=fill(zero(t),n+1,r)
  for n₀ in 0:n
    partitiongrid[n₀+1,1]=one(t)
    for r₀ in 2:min(r,n₀)
      partitiongrid[n₀+1,r₀]=partitiongrid[n₀,r₀-1]+partitiongrid[n₀+1-r₀,r₀]
    end
  end
  return sum(partitiongrid[n+1,1:r])
end
"""
Returns an array of tuples of all of the
wholly positive integers whith a < maxA
"""
function all_whpstvi(maxA,D)
  WP = [(1,0)]
  for a = 2:maxA
    for b = 0:floor(Int, a / √D)
      push!(WP,(a,b))
    end
  end
  WP
end

"""
Returns true of (a,b) does not exceed (c,d) in the field
Q(√D)
"""
function ds_not_exceed((a,b),(c,d),D)
  t = a + b*√D
  t̄ = a - b*√D
  r = c + d*√D
  t <= r && t̄ <= r
end
export ds_not_exceed

# the following three algorithms are rather brute-force
# TODO: figure out a smarter way?

"""
Returns an array of all wholly positive numbers which do not exceed
(c,d)
"""
function all_ds_not_exceed((c,d),D)
  bound = ceil(Int, c + d*√D)
  filter(((x,y),) -> ds_not_exceed((x,y),(c,d),D), all_whpstvi(bound,D))
end
export all_ds_not_exceed

"""
Given a wholly positive integer, gives the next wholly positive integer
with respect to the enumeration that uses "does not exceed"
"""
function next_whpstvi((a,b),D)
  rord((x,y)) = x + y*√D
  bar((x,y)) = (x,-y)
  r = rord((a,b))
  nextreal = b == 0 ? a+1 : ceil(Int, r)

  eligible = (s,) -> rord(bar(s)) < r && r < rord(s)

  ELG = filter(eligible, all_whpstvi(a,D))
  push!(ELG,(nextreal,0))

  minind = argmin(rord.(ELG))
  ELG[minind]
end
export next_whpstvi

"""
gives the greatest wholly positive integer that is less than (a,b)
and does not exceed (a,b)
"""
function largest_whpstvi_lt((a,b),D)
  (a,b)==(1,0) && throw(ArgumentException("No wholly positive numbers less than 1"))

  rord((x,y)) = x + y*√D
  r = rord((a,b))

  eligible = (s,) -> rord(s) < r && ds_not_exceed(s,(a,b),D)
  ELG = filter(eligible, all_whpstvi(2*a,D))
  maxind = argmax(rord.(ELG))
  ELG[maxind]
end

"""
Returns an iterator of all the wholly positive numbers in the
lattice { n - (k*r + q*r̄) }
"""
function whpstv_lattice_lt(n,r,D)
  bar((a,b)) = (a,-b)
  add((a,b),(c,d)) = (a+c,b+d)

  pts = Set([n])

  #dont have a formula for the bounds of this for loop
  k = 0
  while true
    q = 0
    while true
      newN = add(n, -1 .* (add(k .* r, q .* bar(r))))
      if newN == (0,0) || is_wholly_positive(newN...,D)
        push!(pts, newN)
      else
        break
      end
      q += 1
    end
    q == 0 && break
    k += 1
  end

  pts
end

#"""
#Recursively compute the partition number using the recursive
#formula for the number of partitions with parts which are
#less than (c,d)
#"""
@memoize function partition_number_leq(n,r,D)
  #println(" "^(spac), "starting p$n, r=$r")
  #global spac += 1

  # the order of these base cases matters
  n == (0,0) && return 1#(spac -= 1; println(" "^(spac), "p(0), got 1"); return 1)
  !is_wholly_positive(n...,D) && return 0#(spac -= 1; return 0)
  if r == (1,0)
    (_,b) = n
    if b == 0
      return 1
      #spac -= 1; println(" "^(spac), "p$n, 1 part, got 1"); return 1
    else
      return 0
      #spac -= 1; return 0
    end
  end

  total = 0

  s = largest_whpstvi_lt(r,D)
  for nminusl in whpstv_lattice_lt(n,r,D)
    total += partition_number_leq(nminusl,s,D)
  end

  #global spac -= 1
  #println(" "^(spac), "p$n, r=$r got $sum")
  total
end

"""
Gives the partition number for a + b√D using the recursive recurrence formula
"""
function partition_number_recurs(a,b,D,allpositive=false)
  allpositive && throw(ArgumentException("Ok++ not supported by this algorithm"))
  partition_number_leq((a,b),(a,b),D)
end

# RECURSIVE ALGORITHM USING THE GENERATINGFUNCTIONOLOGY METHOD

#### The following code was copied from RosettaCode.com
#### https://rosettacode.org/wiki/Proper_divisors#Julia
####
#### ... and it was adapted for the purposes of this project.
####
#### This code thus uses the GNU Free Documentation Liscense
#### https://www.gnu.org/licenses/old-licenses/fdl-1.2.html

function properdivisors(n)
    0 < n || throw(ArgumentError("number to be factored must be ≥ 0, got $n"))
    1 < n || return Integer[]
    !isprime(n) || return [1, n]
    f = factor(n)
    d = [1]
    for (k, v) in f
        c = [k^i for i in 0:v]
        d = d*c'
        d = reshape(d, length(d))
    end
    sort!(d)
    return d[1:end]
end

#### This ends the code licesned under the GNU Free Documentation Liscense

function σ₁(n)
  n == 1 && return 1
  sum(properdivisors(n))
end

#"""
#Gets the partition number in using the gfology-style
#generating function
#"""
@memoize Dict function partition_number_gfology(a,b,D,allpositive=false,data=nothing,verbose=false)
  allpositive && throw(ArgumentException("Ok++ not yet supported"))

  (a,b) == (0,0) && return 1
  !is_wholly_positive(a,b,D) && return 0
  (a,b) == (1,0) && return 1
  b < 0 && return partition_number_gfology(a,-b,D,allpositive,data)

  if data != nothing
    a <= size(data,2) - 1 && ( return data[b+1,a+1] )
  end

  if D == 2 # remove this restriction later; we need the unit
    (a₁,b₁) = minelement(pellClassFor(a,b,D,(3,2),a))
    (a₁,b₁) != (a,b) && ( return partition_number_gfology(a₁,b₁,D,allpositive,data) )
  end

  # auxillary functions because I'm too lazy to make a datatype
  norm((a,b)) = a^2 - D*(b^2)
  bar((a,b)) = (a,-b)
  times((a,b),(c,d)) = (a*c + b*d*D, b*c + a*d)
  gcid((a,b)) = gcd(a,b)

  n = (a,b)
  whollyLessThanN = TotallyLess.listPoints(a,b,D,true)
  total = (0,0)

  verbose && ( l = length(whollyLessThanN); println("$l possible parts") )

  for m in reverse(whollyLessThanN)
    p_nminusm = partition_number_gfology((n .- m)...,D,allpositive,data)
    g = gcid(m)

    total = @. total + div(m,g) * σ₁(g) * p_nminusm

    (verbose) && ( i = n .- m; println("$i ... tot=$total") )
  end

  verbose && ( println("grand total is... $total") )

  div.(times(bar(n), total), norm(n))[1]
end
partition_number = partition_number_gfology
export partition_number

#"""
#Gets the partition number in using the gfology-style
#generating function
#
# this one uses rational numbers in the hopes that they will require
# less big-numbered calculations, but it turns out that it doesn't
# really help according to timing tests I've done
#"""
@memoize function partition_number_gfology_ratnl(a,b,D,allpositive=false)
  allpositive && throw(ArgumentException("Ok++ not yet supported"))

  (a,b) == (0,0) && return 1
  !is_wholly_positive(a,b,D) && return 0
  (a,b) == (1,0) && return 1
  b < 0 && return partition_number_gfology_ratnl(a,-b,D,allpositive)

  # auxillary functions because I'm too lazy to make a datatype
  norm((a,b)) = a^2 - D*(b^2)
  bar((a,b)) = (a,-b)
  times((a,b),(c,d)) = (a*c + b*d*D, b*c + a*d)
  gcid((a,b)) = gcd(a,b)

  n = (a,b)
  whollyLessThanN = TotallyLess.listPoints(a,b,D,true)
  total = (0//1,0//1)
  norminv = 1//norm(n)

  for m in whollyLessThanN
    p_nminusm = partition_number_gfology_ratnl((n .- m)...,D,allpositive)
    g = gcid(m)

    total = @. total + norminv*(m * 1//g * σ₁(g) * p_nminusm)

  end

  convert(Int, times(bar(n), total)[1])
end

# DYNAMIC ALGORITHM USING GENERATING FUNCTIONOLOGY

function calc_partition_num(a,b,D,recurs,allTtllyLess)
  (a,b) == (0,0) && return 1
  !is_wholly_positive(a,b,D) && return 0
  (a,b) == (1,0) && return 1

  b = abs(b)

  if D == 2 # remove this restriction later; we need the unit
    (a₁,b₁) = minelement(pellClassFor(a,b,D,(3,2),a))
    (a₁,b₁) != (a,b) && ( return recurs(a₁,b₁) )
  end

  # auxillary functions because I'm too lazy to make a datatype
  norm((a,b)) = a^2 - D*(b^2)
  bar((a,b)) = (a,-b)
  times((a,b),(c,d)) = (a*c + b*d*D, b*c + a*d)
  gcid((a,b)) = gcd(a,b)

  n = (a,b)
  total = (0,0)

  for m in reverse(allTtllyLess(a,b,D,true))
    p_nminusm = recurs((n .- m)...)
    #l = n .- m
    g = gcid(m)

    total = @. total + div(m,g) * σ₁(g) * p_nminusm
  end

  div.(times(bar(n), total), norm(n))[1]
end

"""
To start anew, data must be an array
full of missings but which has data[1,1] == 1
at the intex (1,1).
"""
function partition_number_dynamic(a,b,D,allpositive=false,data=nothing,
                                  verbose=false,allTtllyLess=TotallyLess.listPoints)
  allpositive && ( throw(ArgumentError("Unsupported!")) )
  data == nothing && ( throw(ArgumentError("This algorithm requires data.")) )
  # assumes p(0,0) is 1!!!
  p(a,b) = data[abs(b)+1,a+1]

  val = p(a,b)
  !ismissing(val) && ( return (data,val) )

  for (c,d) in allTtllyLess(a,b,D,true)
    !ismissing(p(c,d)) && ( continue )

    # calculate the partition number for (c,d)
    data[abs(d)+1,c+1] = calc_partition_num(c,d,D,p,allTtllyLess)
    verbose && ( nm = p(c,d) ; println("($c,$d): $nm thread $(Threads.threadid())"))
  end

  (data,p(a,b))
end

end#module
