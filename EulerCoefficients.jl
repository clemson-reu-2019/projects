module EulerCoefficients

using DelimitedFiles 

"""
Returns true if x + y√D is wholly positive
"""
function is_wholly_positive(x,y,D)
  0 < x && D*y^2 < x^2
end
export is_wholly_positive

# PERSISTANT STORAGE OF EULER COEFFICIENTS

"""
Contains the power series coefficients which have been pre-calculated
for the product expansion of (the inverse of) the generating function
for p(n). Data is categorized by Ok++ vs Ok+ and by which square root
we choose to adjoin.

DO NOT MODIFY this dictionary. It is just a copy of data which
has been saved to files. To add more coefficients, update the files
and call load_data()
"""
COEF_DATA = Dict{Tuple{Int,Bool},Matrix{Union{Int, Missing}}}()

reset_coefs() = (global COEF_DATA = Dict{Tuple{Int,Bool},Matrix{Union{Int, Missing}}}())

"""
Loads data from sqrt*.dat files into 
the dictionary COEF_DATA.
"""
function load_data()
  # Delete current data, because it is always a copy from the files

  reset_coefs()
  frgx = r"sqrt(.*)\.dat"
  for file in filter(f -> occursin(frgx, f), readdir())
    idstr = match(frgx,file).captures[1]

    if length("plusplus") < length(idstr) && idstr[end-7:end] == "plusplus"
      allpositive = true
      idstr = idstr[1:end-8]
    else
      allpositive = false
    end

    D = parse(Int, idstr)
	M = convert(Matrix{Any}, readdlm(file))
	M[M .== "missing"] .= missing
	M = convert(Matrix{Union{Int, Missing}},M)
    COEF_DATA[(D,allpositive)] = M
  end
end

"""
Takes the arrays M and N and returns
and array which contains the values of both
M and N at the indices of M and N, and missing
everywhere else.
Thus, the arrays M and N have been "merged"

For this to work, the matrices must agree where they 
have common indices.

E.g.
M = 
[1  2
 3  4
 5  6]
N = 
[1  2  3
 3  4  5]
merge_arrays(M,N) = 
[ 1  2  3       
 3  4  5       
 5  6   missing]
"""
function merge_matrices(M,N)
  sizes = [size(A,l) for A=(M,N), l=(1,2)]

  for j = 1:minimum(sizes[:,1])
    for k = 1:minimum(sizes[:,2])
      if (M[j,k] !== missing && N[j,k] !== missing) && M[j,k] != N[j,k]
        throw(ArgumentError("""Error: cannot merge two matrices
							Reason: they have different entries at $j,$k"""))
      end
    end
  end

  m1 = maximum(sizes[:,1])
  m2 = maximum(sizes[:,2])
  MN = Union{Missing, eltype(M)}[missing for n=1:m1, m=1:m2]
  
  # the following overrwrites some values, but
  #   this code shouldn't need to be too performant
  MN[1:size(N,1),1:size(N,2)] = N
  MN[1:size(M,1),1:size(M,2)] = M
  MN
end

"""
Saves the coefficient matrix into a persistant store, 
with the following naming convention:

sqrt2.dat - holds data for Q(√2)
sqrt3.dat - holds data for Q(√3)
...
...
sqrtN.dat - holds data for Q(√N)

For files which denote coefficients for Ok++, 
the postfix 'plusplus' will be added:
e.g.
sqrt2plusplus.dat
"""
function incorporate_coefficients(M,D,allpositive)
  load_data()
  if !haskey(COEF_DATA, (D, allpositive))
    new_coefs = M
  else
    C = COEF_DATA[(D, allpositive)]
    new_coefs = merge_matrices(M,C)
  end

  filename = "sqrt$D"
  allpositive && (filename *= "plusplus")
  filename *= ".dat"

  writedlm(filename, new_coefs)

  load_data()
end

"""
Returns the euler coeffieicent corresponding
to a certain wholly positive number.
"""
function euler_coef(a,b,D,allpositive)
  !is_wholly_positive(a,b,D) && return 0

  key = (D,allpositive)

  # make sure the coeficient data exists
  if !haskey(COEF_DATA, key)
    load_data()
    if !haskey(COEF_DATA, key)
      print("Cannot find coefficient data for $D, $allpositive")
      return 0
    end
  end

  coefs = COEF_DATA[(D,allpositive)]
  eul = coefs[a+1,convert(Int,abs(b))+1]
  eul .=== missing && println("Error: Coefficient missing for $a + $b√$D")
  eul
end
export euler_coef

end#module
