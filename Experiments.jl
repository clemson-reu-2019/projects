module Experiments

include("./PartitionsGen.jl")
include("./EulerCoefficients.jl")
include("./QuadraticPartitions.jl")
include("./PellClasses.jl")

using .PartitionsGen
using .EulerCoefficients
using .QuadraticPartitions
using .PellClasses
#using .EulerCoefficients
using SymPy
using AbstractAlgebra
using DelimitedFiles
using Plots
using PlotUtils
using LsqFit
using Primes
using OffsetArrays
using BitIntegers
using HomotopyContinuation

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
  xᵢ₊₁ = simplify(xᵢ * Sympy.expand(unit^i))
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
  lines =
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

"""
(a,b) / (c,d)
"""
function quadratic_div((a,b),(c,d),D)
  cd_norm = c^2 - D*d^2
  cd_inv = (c // cd_norm,-d//cd_norm)
  quadratic_mult((a,b),cd_inv,D)
end

"""
quadratic_power((a,b),N,D)

raises (a,b) to the Nth power in Q(√D)

Note that this is implemented naively, so that it only
works for N an integer.
"""
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

function modelratnlpower(ydata)
  @. model(x,c) = c[1]*x^(c[2])
  curve_fit(model,1:length(ydata),ydata,[0,0.5])
end

function modelratnlpowernoconst(ydata)
  @. model(x,c) = x^(c[1])
  curve_fit(model,1:length(ydata),ydata,[0.5])
end

function modelpartitionquotient(xrange,ydata)
  @. model(x,c) = c[1]*x^(2 / (2 + c[2])) + c[3]
  curve_fit(model,xrange,ydata,[1.0,1.0,0])
end

function modeltwoplusϵ(ϵ,xrange,ydata)
  @. model(x,c) = c[1]*x^(2 / (2 + ϵ))
  curve_fit(model,xrange,ydata,[1.0])
end

function fitpartitionquotient(range,data)
  fit = modelpartitionquotient(range,log.(data[range]))
  c = coef(fit)
  f(x) = c[1]*x^(2 / (2 + c[2])) + c[3]
  (f,c[2])
end

function highestBFor(a,D)
  floor(Int, a / √D)
end

function partition_nums_for_classes(classes,D,p)
  part_nums = zeros(typeof(p(1,0,D)),length(classes))
  for i = 1:length(classes)
    minind = argmin(map(t -> t[1],classes[i]))
    part_nums[i] = p(classes[i][minind]...,D)
  end
  part_nums
end

"""
comparepellclasses(N,D,unit,allclasses,partition_number)
This method finds counterexamples to the Pell Class conjecture

Example:

comparepellclasses(100,3,(2,1),nothing,p)

calculates the counterexamples which have norm less than 100,
which reside in the field Q(√3), where the partition function
p is used to calculate partition numbers. 
The allclasses parameter is 
to put a pre-calculated list of pell classes, if you want.
Putting "nothing" make the function calculate its own pell classes.

The quantity N is a norm, and this will search for counterexamples
for classes which have norm less than N.
"""
function comparepellclasses(N,D=2,unit=(3,2),allclasses=nothing,partition_number=QuadraticPartitions.partition_number)
  allclasses == nothing && (allclasses = findpellclasses(N,D,unit))
  println("Evaluating all pell classes for counterexample")
  outliers = filter(pair -> 1 < length(pair[2]),allclasses)

  counterexes = []
  for (norm,classes) in outliers
    # minimum a
    part_nums = partition_nums_for_classes(classes,D,partition_number)
    #part_nums = zeros(typeof(N),length(classes))
    #for i = 1:length(classes)
    #  minind = argmin(map(t -> t[1],classes[i]))
    #  part_nums[i] = partition_number(classes[i][minind]...,D)
    #end
    #println("$part_nums")
    if !allunique(part_nums)
      println("Counterexample to PC conjecture found! $norm")
      push!(counterexes, norm)
    end
  end
  eval = x -> outliers[x]
  (counterexes, eval.(counterexes))
end

function all_PC_counterexamples(N,A,partition_number)
  ALL = []
  for a in A
    println()
    println("Starting for d=$a")
    PCs = comparepellclasses(N,a,fund_unit(a),nothing,partition_number)
    push!(ALL,PCs)
  end
  ALL
end

function find_norminv_counterex(N,D,unit,allclasses,partition_number=QuadraticPartitions.partition_number)
  allclasses == nothing && (allclasses = findpellclasses(N,D,unit))
  println("Evaluating all pell classes for counterexample")
  pnums = fill(Array{typeof(partition_number(1,0,D)),1}(), N)

  for (norm,classes) in allclasses
    pnums[norm] = partition_nums_for_classes(classes,D,partition_number)
  end

  println("counterexamples are of the form (norm1,norm2,partition_number)")


  counterexes = []
  #maybe not the most efficient but oh well.

  # this should be O(n^2), because the 
  #  number of pell classes is so small in comparison 
  #  to the number of possible norms.

  # also this will print any duplicates from the same
  # norm twice, a quirk that makes it technically incorrect but
  # possibly useful for identifying this class of examples.
  for i in 1:length(pnums)
    for j in i:length(pnums)
      for i_pnum in 1:length(pnums[i])
        for j_pnum in 1:length(pnums[j])
          (i,i_pnum) == (j,j_pnum) && continue
          p = pnums[i][i_pnum]
          q = pnums[j][j_pnum]
          if p == q 
            print("Counterexample Found!")
            println(" ($i,$j,$p)")
            push!(counterexes, (i,j,p))
          end
        end
      end
    end
  end

  (counterexes, allclasses)
end

function fund_unit(d)
  F = Dict(
           2 => (3,2),
           3 => (2,1),
           5 => (9,4),
           6 => (5,2),
           7 => (8,3),
           10 => (19,6),
           11 => (10,3),
           13 => (649,180),
           14 => (15,4),
           15 => (4,1),
           17 => (33,8),
           19 => (170,39),
           21 => (55,12),
           22 => (197,42),
           23 => (24,5),
           26 => (51,10),
           #29 => (9801,1820),
           30 => (11,2),
           31 => (1520,273),
           33 => (23,4),
           34 => (35,6),
           35 => (6,1),
           70 => (251,30),
          )
  F[d]
end

function num_pn_norminv_counterexes(d,p)
  CE = Experiments.find_norminv_counterex(1000,d,fund_unit(d),nothing,p)[1]
  p3((x,y,z)) = z
  length(unique(p3.(CE)))
end

all_unique_keys_sorted(dict) = sort(unique(collect(keys(dict))))
  

"""
gives the numbers in X which are in the set of primes from 1 to N, 
mod n
"""
function primesmod(N,n,X)
	return filter!((x)->mod(x,n) in X,Primes.primes(N))
end

function testmod(n,X)
    N=unique([1:n;].^2 .%n)
	for x in X
		unique2=true
		uniquenegative2=true
		for n₁ in N
			for n₂ in N
				if mod(n₁-x*n₂,n)==2 && unique2==true
					println("Possible solution 2 for $x mod $n")
					unique2=false
				end
				if mod(n₁-x*n₂,n)==(n-2) && uniquenegative2==true
					println("Possible solution -2 for $x mod $n")
					uniquenegative2=false
				end
			end
		end
		if unique2==true ||  uniquenegative2==true
			println("At Most One Case Above")
    end
  end
end

function partitions_parts_distr(N)
  DIST = zeros(Int,N,N)
  for n = 1:N
    for r = 1:N
      DIST[n,r] = partition_number_leq(n,r)
    end
  end
  DIST
end

function all_not_exceed_bool(a,b,D,N=0)
  N == 0 && (N = a+1)
  A = zeros(Int, N+1, highestBFor(N,D)+1)
  for (x,y) in QuadraticPartitions.all_whpstvi(N,D)
    if ds_not_exceed((x,y),(a,b),D)
      A[x+1,y+1] = 2
    else
      A[x+1,y+1] = 1
    end
  end
  A[a+1,b+1] = 3
  A'
end

function find_congruences(p,base,mults,mods,offsets,N)
  length(mults) != length(mods) && ( throw(ArgumentError("there must be the same number of multaplicative factors and modulos")) )

  for i in 1:length(mults)
    m = mults[i]
    md = mods[i]
    for offset in offsets
      scales = m .* collect(0:N) .+ offset
      (a,b) = base
      res = zeros(Int,length(scales))
      for i in 1:length(scales)
        res[i] = convert(Int, p((scales[i] .* base)...) .% md)
      end
      if length(unique(res)) == 1
        ans = res[1]
        println("Found congruence: p($m*k + $offset) = $ans (mod $md)")
      end
    end
  end
end

find_congruences(p,mults,mods,offsets,N) = find_congruences(p,(1,0),mults,mods,offsets,N)

function norm_representations(N,D)
  norms = [a^2 - D*b^2 for a=0:N,b=0:N]
  norms[norms .< 0] .= 0
  sort!(unique(norms))
end

function sums_of_norms(N,D,unit,num=2,allowneg=false)
  norms = [a^2 - D*b^2 for a=0:N,b=0:N]
  if !allowneg
    norms[norms .< 0] .= 0
  end
  norms = unique(norms)

  range = minimum(norms):maximum(norms)
  sums = zeros(Int, fill(range,num)...)
  sumset = Set{Int}()

  for n1 in norms
    for n2 in norms
      if num == 3
        for n3 in norms
          sums[n1,n2,n3] = n1 + n2 + n3
          push!(sumset,n1 + n2 + n3)
        end
      else
        sums[n1,n2] = n1 + n2
        push!(sumset,n1 + n2)
      end
    end
  end

  if !allowneg
    bound = floor(Int, N^2 / (unit[1] + unit[2]*√D))
    println("Found all sums at least to $bound")
  end
  (sums,sort(collect(sumset)))
end

function parts_over_partitions_distr(n,p=nothing)
  p == nothing && (p = partition_number_leq)
  leqs = p.(n,1:n)
  distr = zeros(Int,n)
  distr[1] = leqs[1]
  for i = 2:n
    distr[i] = leqs[i] - leqs[i-1]
  end
  distr
end

function maxnumparts_to(N)
  argmax.(parts_over_partitions_distr.(1:N))
end

leftside(m,a,b) = m*(m-1)*a^(m-2)*b^2

rightside(m,a,b) = (a+b)^m - 2a^m + (a-b)^m

function testineq(N,aRange,bRange,D)
  for i = 1:N
    for j = 1:highestBFor(i,D)
      for a in aRange
        for b in bRange
          a < b && continue
          m = i + j*√D
          if leftside(m,a,b) > rightside(m,a,b)
            println("Bad: m=$i + $j√$D, a=$a, b=$b")
            return
          end
        end
      end
    end
  end
end

function factorial(a,b,D)
  prod = BigInt(1)
  for (x,y) in all_ds_not_exceed((a,b),D)
    if y == 0
      prod *= x
    else
      prod *= x^2 - D*y^2
    end
  end
  prod
end

function injectivetest(a,b,d,bool=false)
  Q=QuadraticPartitions.quad_partitions(a,b,d,bool)
  N=Array{Tuple{Int,Int}}[]
  for Q₁ in Q
    N₁=Tuple{Int,Int}[]
    for i in 1:size(Q₁)[1]
      for j in i:size(Q₁)[1]
        ps=(Q₁[i,1]*Q₁[j,1],2*Q₁[i,2]*Q₁[j,2])
        push!(N₁,ps)
      end
    end
    sort!(N₁,by =x -> x[2])
    push!(N,N₁)
  end
  return length(unique(N))==length(Q)
end

function injectivecounter(n,u,d)
  PositiveUnits=Tuple[]
  NegativeUnits=Tuple[]
  counter=0
  for a in 1:u
    for b in 0:Int(floor(a/√d))
      if a^2-d*b^2==1
        push!(PositiveUnits,(a,b))
        push!(NegativeUnits,(a,-b))
      end
    end
  end
  for a₁ in 1:n
    for b₁ in 1:Int(floor(a₁/√d))
      for a₂ in 1:n
        for b₂ in 1:Int(floor(a₂/√d))
          for u₁ in PositiveUnits
            for u₂ in NegativeUnits
              if u₁[1]≠u₂[1]
                if a₁+a₂ == u₁[1]*a₁+d*u₁[2]*b₁+u₂[1]*a₂+d*u₂[2]*b₂ && b₁+b₂ == a₁*u₁[2]+b₁*u₁[1]+a₂*u₂[2]+b₂*u₂[1]
                  println((a₁+a₂,b₁+b₂))
                  println((a₁,b₁))
                  println((a₂,b₂))
                  println(u₁)
                  println(u₂)
                  println("")
                  counter=counter+1
                end
              end
            end
          end
        end
      end
    end
  end
  println("$counter Counterexamples")
end

"""
the following has become a bit of sphaghetti,
in an effort to support arrays which might have
missing values, might have overflowed values,
and might have missing values set to -1. Hopefully
this whole code won't be run much more, so this
won't need to be messed with too much.
"""
function process_grid(t::Type,filename,D)
  println("Reading $filename...")
  grid = readdlm_withmissings(t, filename)

  # assumes that if there is a negative value in
  #  the array at all, then there are no missings
  #  in between the first row and the first
  #  negative value.
  ovflow_occurs = foldl((accum,x) -> accum || (ismissing(x) ? false : x), grid .< 0)
  if ovflow_occurs
    ovflowind = findfirst(grid[1,:] .< 0)
    ovflowind == nothing && return grid
    reducedgrid = grid[1:highestBFor(ovflowind-1,D),1:(ovflowind-1)]
  else
    reducedgrid = grid
  end
  mgrid = convert(Array{Union{t,Missing}}, reducedgrid)
  # account for possibilities of missings
  neg1inds = mgrid .== t(-1)
  neg1inds[ismissing.(neg1inds)] .= false
  mgrid[convert(Array{Bool}, neg1inds)] .= missing
  #count(mgrid .< 0) != 0 && println("Found some error $D")
  
  # forcefully set p(0)
  mgrid[1,1] = 1

  mgrid
end

function process_all_files(t::Type)
  # 3 is the files, as opposed to the directories or the prefix
  filenames = collect(walkdir("./data"))[1][3]
  ds = parse.(Int,map(x -> x[1], match.(r"grid(\d*)",filenames)))

  grids = process_grid.(t,"./data/" .* filenames,ds)
  gridict = Dict()
  for i in 1:length(ds)
    d = ds[i]
    println("Merging in d = √$d")
    if d in keys(gridict)
      gridict[d] = EulerCoefficients.merge_matrices(gridict[d],grids[i])
    else
      gridict[d] = grids[i]
    end
  end
  #Dict([(ds[i],grids[i]) for i in 1:length(ds)])
  gridict
end

function testdivision(t::Type,alg,N,lim=100001)
  numers = rand(one(t):lim,N)
  denoms = rand(one(t):1000,N)
  outputs = zeros(t,N)
  for i in 1:N
    outputs[i] = alg(numers[i],denoms[i])
  end
  outputs
end

function testplusidentity(N,p)
  for i in 1:N
    p(2*i) >= p(i)^2 && ( print("$i ") )
  end
end

function plusidentitydata(N,p)
  [ (p[i]^2) / p[2*i] for i in 1:N]
end

function plusidentitydatacube(N,p)
  [ (p(i)^3) / p(3*i) for i in 1:N]
end

function testdynamicalg(a,b,D,table)
  data = Array{Union{Missing, Int1024}}(missing, a+1,a+1)
  firstN = table[1:a,1:a]
  data[axes(firstN)...] = firstN
  @time QuadraticPartitions.partition_number_dynamic(Int1024(a),Int1024(b),D,data)
end

function testdynamicalgdiffdata(a,a₀,D,table)
  data = Array{Union{Missing, Int1024}}(missing, a+1,a+1)
  firstN = deepcopy(table[1:(a₀+1),1:(a₀+1)])
  data[axes(firstN)...] = firstN
  data
end

function create_emptygrid(size,D)
  bSize = highestBFor(size,D)
  grid = Matrix{Union{Int1024,Missing}}(missing,bSize+1,size+1)
  grid[1,1] = 1
  grid
end

function best_partition_function()
  p(i,j,D,allpos,data) = QuadraticPartitions.
                  partition_number_dynamic(i,j,D,allpos,data,true) 
end

function best_partition_function(f)
  p(i,j,D,allpos,data) = QuadraticPartitions.
                  partition_number_dynamic(i,j,D,allpos,data,true,
                       QuadraticPartitions.TotallyLess.listPointsConductor(f))
end
  

# I've been running tests with this:
# mydata = @time Experiments.QuadraticPartitions.partitions_grid(Int1024,1000,D,false,p,grid)
# where p and grid are created using best_partition_function() and create_emptygrid()

#julia> A = [  2
#         3
#         5
#         6
#         7
#        10
#        11
#        13
#        14
#        15
#        17
#        19
#        21
#        22
#        23
#        26
#         31
#        33
#        34
#        35
#       70]

function compare_normcexes_to_units(grids,Ds)
  p(a,b,D) = grids[D][b+1,a+1]
  normcexes = Experiments.num_pn_norminv_counterexes.(Ds,p)

  p1((x,y)) = x
  as = p1.(Experiments.fund_unit.(Ds))

  p2((x,y)) = y
  bs = p2.(Experiments.fund_unit.(Ds))

  plot(log.((as.^2) + (bs .^ 2)), normcexes,seriestype=:scatter,
       xlabel="size of fundamental unit: log (a^2 + b^2)",
       ylabel="number of counterexamples to norm invariance")

  res = Experiments.linreg(log.((as.^2) + (bs .^ 2)),normcexes)
  plot!(x -> coef(res)[1] *x + coef(res)[2])
end

# Experiments.QuadraticPartitions.quad_partitions_decomp(8,0,2)
# filter(arr -> reduce((accum,(a,b)) -> (accum && (b % 2 == 0)), arr, init=true), eightPartitions)

function partitions_decomp_restrict(a,b,D,conductor)
  b % conductor == 0 || (return missing)
  dec_partitions = QuadraticPartitions.quad_partitions_decomp(a,b,D)
  filter(arr -> reduce((accum,(a,b)) -> (accum && (b % conductor == 0)), 
                       arr, init=true), 
         dec_partitions)
end

function partition_number_decomp_restrict(a,b,D,conductor)
  b % conductor == 0 || (return missing)
  length(partitions_decomp_restrict(a,b,D,conductor))
end 



end#module
