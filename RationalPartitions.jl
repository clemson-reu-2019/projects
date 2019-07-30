module RationalPartitions

include("./PartitionsGen.jl")
include("./QuadraticPartitions.jl")
using .PartitionsGen
using .QuadraticPartitions
using BitIntegers

    function multiplicativetest2()
        bound=2324
        p₂=2
        counterbool=false
        for n in 2:bound
            pₙ=QuadraticPartitions.partition_number_leq_array(Int1024,Int1024(n),Int1024(n))
            pₓ=QuadraticPartitions.partition_number_leq_array(Int1024,Int1024(2*n),Int1024(2*n))
            if (p₂*pₙ)≤pₓ
                println("Checked m=2 and n=$n: NO COUNTEREXAMPLE")
            end
            if (p₂*pₙ)>pₓ
                println("COUNTEREXAMPLE for m=2 and n=$n")
                counterbool=trueinc
            end
        end
        if !counterbool
            print("No Counterexamples")
        end
    end

    function multiplicativetest()
        bound=293
        counterbool=false
        for m in 3:bound
            pₘ=QuadraticPartitions.partition_number_leq_array(Int1024,Int1024(m),Int1024(m))
            for n in 3:Int(floor(asymptoticbound(m)))
                pₙ=QuadraticPartitions.partition_number_leq_array(Int1024,Int1024(n),Int1024(n))
                pₓ=QuadraticPartitions.partition_number_leq_array(Int1024,Int1024(m*n),Int1024(m*n))
                if (pₘ*pₙ)≤pₓ
                    println("Checked m=$m and n=$n: NO COUNTEREXAMPLE")
                end
                if (pₘ*pₙ)>pₓ
                    println("COUNTEREXAMPLE for m=$m and n=$n")
                    counterbool=true
                end
            end
        end
        if !counterbool
            print("No Counterexamples")
        end
    end

    function asymptoticbound(m)
        return 1/((m-2)^2)*(10*√2*√(m*(m+3)^2/((m-2)^4))*m^2+2*m^2-40*√2*√(m*(m+3)^2/((m-2)^4))*m+37*m+40*√2*√(m*(m+3)^2/((m-2)^4))+18)
    end
end
