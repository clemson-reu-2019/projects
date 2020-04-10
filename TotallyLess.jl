module TotallyLess

"""
Gives a list of the numbers which are totally
less than a + b*√D

the output is in ascending order of a, then b
"""
function listPoints(a,b,D,excludezero=false,conductor=1)
  List=[]
  b₀=floor((a-b*√D)/2)
  b₁=floor((a+b*√D)/2)
  x=0
  y=a+1
  if excludezero==true
      x=1
  end
  for a₁ in x:a
      if a₁<=b₀
          for b₂ in Int(ceil(-a₁/√D)):Int(floor(a₁/√D))
              push!(List,(a₁,b₂))
          end
      end
      if b₀<a₁<=b₁
          for b₂ in Int(ceil((a₁-a)/√D+b)):Int(floor(a₁/√D))
              push!(List,(a₁,b₂))
          end
      end
      if a₁>b₁ && a₁<y
          for b₂ in Int(ceil((a₁-a)/√D+b)):Int(floor((-a₁+a)/√D+b))
              push!(List,(a₁,b₂))
          end
      end
  end
  isInOrder((a,b)) = b == 0 || (b % conductor == 0)
  return filter!(isInOrder, List)
end

"""
A curried version of listPoints if you want to pass around
a function for a fixed conductor.
"""
function listPointsConductor(f)
  lP(a,b,D,twopart) = listPoints(a,b,D,twopart,f)
end

end
