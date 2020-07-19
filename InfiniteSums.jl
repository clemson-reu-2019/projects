
using Nemo

# quadratic field basics

function fund_unit(D)
  D == 2 && return (3,2)
  (0,0)
end

function norm((a,b),D)
  a^2 - D * b^2
end

function realv((a,b),D)
  a + b * √D
end

function conjv((a,b),D)
  a - b * √D
end

# some constants

function b(D)
  π / log(fund_unit(D))
end

function λ(μ,D)
  exp(i * b(D) * log( conjv(μ,D) / realv(μ,D)))
end



# basic summation term expression

function zeta_term(s,n)
  1 / (n^s)
end

#function gamma_term(s)
#  
#end

function heckezeta_term(s,n)
  
end

# main ingredients 

# Nemo provides
# zeta(s), gamma(s), lgamma(s) [log of gamma], 
#   and rgamma(s) [reciprocal of gamma]

#more parameters will probably be needed.

function frakS1_term(n,N)
  
end

function frakS2_term(n,N)
  
end

function frakS3_term(n,N)
  
end

function partialSum(f,maxN)
	sum(f.(1:maxN))
end

function frakSconstant(N,maxNumTerms=1000)
  
end


