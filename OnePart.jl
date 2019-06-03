module OnePart

using Plots

    function specialbd(m,n)
        Maxnorm=Float64[]
        D=Int[]
        maxnorm=0
        for d in 2:m
            if d==Int(ceil(sqrt(d))^2-1)
                push!(D,d)
                maxnorm=multib(d,n)
                push!(Maxnorm,multib(d,n))
                if maxnorm≠1
                    println("Counterexample")
                end
            end
        end
        println(D,Maxnorm)
    end

    function multibd(m,n)
        Maxnorm=Float64[]
        D=Int[]
        for d in 2:m
            if sqrt(d)≠floor(sqrt(d))
                push!(D,d)
                push!(Maxnorm,multib(d,n))
            end
        end
        plot(D,Maxnorm)
    end


    function multib(d,n)
        maxnorm=0
        maxa=0
        maxb=0
        ptest=true
        N=0
        for b in 0:n
            a=Int(ceil(b*√d))
            ptest=true
            for x in 1:floor(b/2)
                if a≥Int(ceil((b-x)*√d)+ceil(x*√d))
                    ptest=false
                end
                N=Norm(a,b,d)
            end
            if ptest==true && N>maxnorm
                maxnorm=N
                maxa=a
                maxb=b
            end
        end
        return maxnorm
    end

    function Norm(a,b,d)
        return a^2-d*b^2
    end

end
