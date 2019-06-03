module Hecke

using Primes

  function multic(n)
    m=M(n)
    y=0
    for N in 1:m
      if mod(m,N)==0
        for s in 0:isqrt(4n)
          if (s^2-4n)!=0
            T=t(s,n)
            for f in 1:T
              if mod(T,f)==0
                x=c(s,f,N,n)
                if x!=0
                    y=y+1
                    println("c($s,$f,$N,$n)=$x")
                end
              end
            end
          end
        end
      end
    end
    println("$y nonzero c values")
  end

  function t(s,n)
    b=0
    for a in 1:isqrt(abs(s^2-4n))
      if mod(s^2-4n,a^2)==0
        b=a
      end
    end
    if mod((s^2-4n)/b^2,4)!=1
      b=b/2
    end
    return b
  end


  function M(n)
  m=1
  for l in 1:4n
    if Primes.isprime(l) && gcd(l,n)==1
      m=m*l^(M(n,l))
    end
  end
  return m
  end

  function M(n,l)
  m=0
  for s in 0:isqrt(4n)
    if s^2-4n!=0
      for d in 1:isqrt(n)
        if mod(abs(s^2-4n),l)==0 || mod(n/d-d,l)==0
          for s in 0:isqrt(4n)
            if s^2-4n!=0
              v=0
              for i in 1:Int(floor(log(l,abs(s^2-4n))))
                if mod(abs(s^2-4n),l^(i))==0
                  v=i
                end
              end
              if v>m
                m=v
              end
            end
          end
          for d in 1:isqrt(n)
            v=0
            for i in 1:Int(floor(log(l,abs(s^2-4n))))
              if mod(abs(s^2-4n),l^(i))==0
                v=i
              end
            end
            if v>m
              m=v
            end
          end
          return m+1
        end
      end
    end
  end
  return m
  end

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
