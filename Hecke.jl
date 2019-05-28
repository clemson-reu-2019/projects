module Hecke

using Primes

  function c(s,f,N,n)
    c=1
    for l in 1:N
      if mod(N,l)==0 && Primes.isprime(l)
        c=c*c2(s,f,N,n,l)
      end
    end
    return c
  end

  function c2(s,f,N,n,l)
  v=0
  b=0
    for i in 1:Int(floor(log(l,N)))
      if mod(N,l^(i))==0
        v=i
      end
    end
    for j in 1:Int(floor(log(l,f)))
      if mod(f,l^(j))==0
        b=j
      end
    end
    a=A(s,n,l,b,v)
    c=a
    if mod((s^2-4*n)/(f^2),l)==0
      b=B(s,n,l,b,v)
      c=a+b
    end
    return c
  end

  function A(s,n,l,b,v)
    L=Int[]
    K=Int[]
    if b>0
      for i in 0:(l^(v+b+1)-1)
        x=s/2+i*l^b
        if mod(x^2-s*x+n,l^(v+2b))==0
          push!(L,x)
        end
      end
    end
    if b==0
      for x in 0:(l^(v+2b)-1)
        if mod(x^2-s*x+n,l^(v+2b))==0
          push!(L,x)
        end
      end
    end
    K=copy(L)
    for j in 1:(length(L))
      for k in 1:(length(L))
        if j!=k && mod(L[j],l^(v+b))==mod(L[k],l^(v+b)) && k<=length(K)
          deleteat!(K,k)
        end
      end
    end
    return length(K)
  end

  function B(s,n,l,b,v)
    L=Int[]
    K=Int[]
    if b>0
      for i in 0:(l^(v+b+1)-1)
        x=s/2+i*l^b
        if mod(x^2-s*x+n,l^(v+2b))==0
          push!(L,x)
        end
      end
    end
    if b==0
      for x in 0:(l^(v+2b+1)-1)
        if mod(x^2-s*x+n,l^(v+2b))==0
          push!(L,x)
        end
      end
    end
    K=copy(L)
    for j in 1:(length(L))
      for k in 1:(length(L))
        if j!=k && mod(L[j],l^(v+b))==mod(L[k],l^(v+b)) && k<=length(K)
          K=deleteat!(K,k)
        end
      end
    end
    return length(K)
  end

end
