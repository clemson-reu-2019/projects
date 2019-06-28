module OnePart
using Plots
using Combinatorics

    function smallb(a,b,d,u)
        a₀=a
        b₀=b
        u₀=(1,0)
        if b>0
            while abs(b₀*u[1]-a₀*u[2])<b₀
                a₀,b₀=(a₀*u[1]-b₀*u[2]*d),(u[1]*b₀-u[2]*a₀)
                u₀=(u[1]*u₀[1]+d*u[2]*u₀[2],u[1]*u₀[2]+u₀[1]*u[2])
            end
        end
        if b<0
            while abs(b₀*u[1]+a₀*u[2])<b₀
                a₀,b₀=(a₀*u[1]+b₀*u[2]*d),(u[1]*b₀+u[2]*a₀)
                u₀=(u[1]*u₀[1]-d*u[2]*u₀[2],u[1]*u₀[2]-u₀[1]*u[2])
            end
        end
        return(((a₀,b₀),u₀))
    end

    function DecompositionLength(n,d,u)
        counter=false
        for b₁ in 0:Int(floor(n/√d))
            for a₁ in Int(ceil(b₁*√d)):n
                if LongDecomposition(a₁,b₁,d,u)≠true
                    println("Counterexample for $a₁+$b₁√$d")
                    counter=true
                end
            end
        end
        if ~counter
            println("No counterexamples")
        end
    end

    function LongDecomposition(a,b,d,u,output=false)
        N=a^2-d*b^2
        decomposition=Tuple[]
        Blocks=Tuple{Int,Int}[]
        um=smallb(a,b,d,u)
        a=um[1][1]
        b=um[1][2]
        u=um[2]
        for a₀ in Int(ceil(√d)):a
            b₀=Int(floor(a₀/√d))
            ptest=true
            if a₀>ceil(b₀*√d)
                ptest=false
            end
            for x in 1:Int(floor(b₀/2))
                if a₀≥Int(ceil((b₀-x)*√d)+ceil(x*√d))
                    ptest=false
                end
            end
            if ptest
                if b>0 && a-a₀≥abs((b-b₀)*√d)
                    push!(Blocks,(a₀,b₀))
                end
                if b<0 && a-a₀≥abs((b+b₀)*√d)
                    push!(Blocks,(a₀,-b₀))
                end
            end
        end
        push!(Blocks,(1,0))
        a₀=a
        b₀=b
        for block in Blocks
            while (a₀-block[1])-abs(b₀-block[2])*√d≥0 && abs(b₀-block[2])*abs(b)==(b₀-block[2])*b
                a₀=a₀-block[1]
                b₀=b₀-block[2]
                push!(decomposition,((block[1]*u[1]+d*block[2]*u[2]),(block[2]*u[1]+block[1]*u[2])))
            end
        end
        if output
            return decomposition
        end
        if sqrt(N/d)>length(decomposition)
            return false
        end
        if sqrt(N/d)≤length(decomposition)
            return true
        end
    end

    function QuickDecompositions(a,b,d)
        Blocks=Tuple{Int,Int}[]
        Decompositions=Array{Tuple{Int,Int}}[]
        push!(Blocks,(1,0))
        ptest=true
        for a₀ in Int(ceil(√d)):a
            b₀=Int(floor(a₀/√d))
            ptest=true
            if a₀>ceil(b₀*√d)
                ptest=false
            end
            for x in 1:Int(floor(b₀/2))
                if a₀≥Int(ceil((b₀-x)*√d)+ceil(x*√d))
                    ptest=false
                end
            end
            if ptest
                if a-a₀>abs((b-b₀)*√d)
                    push!(Blocks,(a₀,b₀))
                    push!(Blocks,(a₀,-b₀))
                end
            end
        end
        n=Int(floor(sqrt(a^2-d*b^2)))
        B=length(Blocks)
        K=Array{Int64,1}[]
        push!(K,digits(0, base = B, pad=n).+1)
        for i in 1:B^n
            push!(K,digits(i, base = B, pad=n).+1)
        end
        for k in K
            decomposition=Tuple[]
            suma=0
            sumb=0
            for j in 1:n
                suma=Blocks[k[j]][1]+suma
                sumb=Blocks[k[j]][2]+sumb
                push!(decomposition,Blocks[k[j]])
                if suma==a && sumb==b
                    sort!(decomposition, by =x -> x[2])
                    push!(Decompositions,decomposition)
                end
            end
        end
        return unique(Decompositions)
    end

    function Decompositions(a,b,d)
        Blocks=Tuple{Int,Int}[]
        Decompositions=Array{Tuple{Int,Int}}[]
        push!(Blocks,(1,0))
        ptest=true
        for a₀ in Int(ceil(√d)):a
            b₀=Int(floor(a₀/√d))
            ptest=true
            if a₀>ceil(b₀*√d)
                ptest=false
            end
            for x in 1:Int(floor(b₀/2))
                if a₀≥Int(ceil((b₀-x)*√d)+ceil(x*√d))
                    ptest=false
                end
            end
            if ptest
                push!(Blocks,(a₀,b₀))
                push!(Blocks,(a₀,-b₀))
            end
        end
        Orderings=permutations(Blocks)
        n=Int(floor(sqrt(a^2-d*b^2)))
        B=length(Blocks)
        K=Array{Int64,1}[]
        push!(K,digits(0, base = n, pad=B).+1)
        for i in 1:n^B
            push!(K,digits(i, base = n, pad=B).+1)
        end
        OrderingsTimes=Tuple[]
        for order in Orderings
            for k in K
                push!(OrderingsTimes,(order,k))
            end
        end
        for bt in OrderingsTimes
            a₀=a
            b₀=b
            decomposition=Tuple[]
            for i in 1:B
                for j in 1:bt[2][i]
                    if (a₀-bt[1][i][1])-abs(b₀-bt[1][i][2])*√d≥0
                        a₀=a₀-bt[1][i][1]
                        b₀=b₀-bt[1][i][2]
                        push!(decomposition,bt[1][i])
                    end
                end
            end
            if a₀==0
                sort!(decomposition, by =x -> x[2])
                push!(Decompositions,decomposition)
            end
        end
        return unique(Decompositions)
    end

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
        maxnorm=0
        D=Int[]
        for d in 2:m
            if sqrt(d)≠floor(sqrt(d))
                push!(D,d)
                maxnorm=multib(d,n)
                push!(Maxnorm,maxnorm)
                if sqrt(d-1)==floor(sqrt(d-1)) && maxnorm≠d
                    print("counterexample")
                end
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
            if ptest && N>maxnorm
                maxnorm=N
                maxa=a
                maxb=b
            end
        end
        return maxnorm
    end

    function FindOneParts(d,n)
        A=Tuple[]
        for b in 0:n
            a=Int(ceil(b*√d))
            ptest=true
            for x in 1:floor(b/2)
                if a≥Int(ceil((b-x)*√d)+ceil(x*√d))
                    ptest=false
                end
                N=Norm(a,b,d)
            end
            if ptest
                push!(A,(a,b))
            end
        end
        return A
    end

    function Norm(a,b,d)
        return a^2-d*b^2
    end

end
