module TotallyLess

"""
Gives a list of the numbers which are totally
less than a + b*√d
"""
function listPoints(a,b,d,twopart=false)
  # make sure output type is same as input type
  t = typeof(a)

  List=[]
  b₀=floor((a-b*√d)/2)
  b₁=floor((a+b*√d)/2)
  x=0
  y=a+1
  if twopart==true
      x=1
  end
  for a₁ in x:a
      if a₁<=b₀
          for b₂ in t(ceil(-a₁/√d)):t(floor(a₁/√d))
              push!(List,(a₁,b₂))
          end
      end
      if b₀<a₁<=b₁
          for b₂ in t(ceil((a₁-a)/√d+b)):t(floor(a₁/√d))
              push!(List,(a₁,b₂))
          end
      end
      if a₁>b₁ && a₁<y
          for b₂ in t(ceil((a₁-a)/√d+b)):t(floor((-a₁+a)/√d+b))
              push!(List,(a₁,b₂))
          end
      end
  end
  return List
end

end
