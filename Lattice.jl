module Lattice

using Plots

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
