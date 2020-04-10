module Lattice
using Plots
using Primes
include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")

function conjugatetest(d,n)
    for b₁ in 0:Int(floor(n/√d))
        for a₁ in Int(ceil(b₁*√d)):n
            for b₂ in 1:Int(floor(n/√d))
                for a₂ in Int(ceil(b₂*√d)):n
                    p₁=Int(Points(a₁*a₂+b₁*b₂*d,b₁*a₂+a₁*b₂,d))
                    p₂=Int(Points(a₁*a₂-b₁*b₂*d,abs(a₁*b₂-b₁*a₂),d))
                    if p₁≠p₂
                        println("For c₁=$a₁+$b₁√$d and c₂=$a₂+$b₂√$d, $p₁ points for c₁c₂ and $p₂ points for c₁̄c₂")
                        if QuadraticPartitions.partition_number(a₁*a₂+b₁*b₂*d,b₁*a₂+a₁*b₂,d,false) == QuadraticPartitions.partition_number(a₁*a₂-b₁*b₂*d,abs(a₁*b₂-b₁*a₂),d,false)
                            println("c₁c₂ and  c₁̄c₂ have the same partition number")
                        end
                        if abs(p₁-p₂)≠2
                            println("Difference not equal to 2")
                        end
                    end
                end
            end
        end
    end
end

"""

"""
function differencetest(d,n,D)
    A=Array{Int64}[]
    for b₁ in 0:Int(floor(n/√d))
        for a₁ in Int(ceil(b₁*√d)):n
            for b₂ in 1:Int(floor(n/√d))
                for a₂ in Int(ceil(b₂*√d)):n
                    p₁=Int(Points(a₁*a₂+b₁*b₂*d,b₁*a₂+a₁*b₂,d))
                    p₂=Int(Points(a₁*a₂-b₁*b₂*d,abs(a₁*b₂-b₁*a₂),d))
                    if abs(p₁-p₂)==D
                        x=abs(p₁-p₂)
                    # println("For c₁c₂=$(a₁*a₂+b₁*b₂*d)+$(b₁*a₂+a₁*b₂)√$d and c₁̄c₂=$(a₁*a₂-b₁*b₂*d)+$(abs(a₁*b₂-b₁*a₂))√$d, difference is $x")
                    # println("For norm $((a₁*a₂+b₁*b₂*d)^2-d*(b₁*a₂+a₁*b₂)^2), difference is $x")
                        push!(A,Primes.factor(Vector, (a₁*a₂+b₁*b₂*d)^2-d*(b₁*a₂+a₁*b₂)^2))
                    end
                end
            end
        end
    end
    return unique(A)
end


    function multiab(d,n)
        A=Float64[]
        B=Float64[]
        m=0
        a₀=0
        b₀=0
        N₀=0
        for b in 0:floor(n/√d)
            for a in Int(ceil(b*√d)):n
                E=Error(a,b,d)
                push!(B,a^2-d*b^2)
                push!(A,Error(a,b,d))
                if abs(E)>abs(m)
                    m=E
                    a₀=a
                    N₀=a^2-d*b^2
                    b₀=b
                end
            end
        end
        println("Calculated error for all a+b√$d with a ranging from 1 to $n and b ranging from 0 to $(floor(n/√d))")
        println("Maximum error $m at n=$a₀+$b₀√$d, N(n)=$N₀")
        plot(B,A)
    end

    function multia(b,d,n)
        A=Float64[]
        B=Float64[]
        m=0
        x=0
        N=0
        for a in Int(ceil(b*√d)):n
            E=Error(a,b,d)
            push!(B,a^2-d*b^2)
            push!(A,Error(a,b,d))
            if abs(E)>abs(m)
                m=E
                x=a
                N=a^2-d*b^2
            end
        end
        println("Calculated error for all a+$b√$d with a ranging from $(Int(ceil(b*√d))) to $n")
        println("Maximum error $m at a=$x, N(n)=$N")
        plot(B,A)
    end

    function Error(a,b,d)
        return Area(a,b,d)-Points(a,b,d)
    end

    function Norm(a,b,d)
        return a^2-d*b^2
    end

    function Area(a,b,d)
        return Norm(a,b,d)/(2*√d)
    end

    function Points(a,b,d)
        sum=0
        b₀=floor((a-b*√d)/2)
        b₁=floor((a+b*√d)/2)
        for i in 1:a
            if i<=b₀
                sum = sum + floor(i/√d)-floor(-i/√d)
            end
            if b₀<i<=b₁
                sum = sum + floor(i/√d)-floor((i-a)/√d+b)
            end
            if i>b₁
                sum = sum + floor((-i+a)/√d+b)-floor((i-a)/√d+b)
            end
        end
        return sum
    end
end
