using Printf
using Plots

function compute_pi(λ;Z=512, maxiter=1024)
    setprecision(Z)

    function Pochhammer(x, n)
        result = big(1)
        for i in 0:n-1
            result *= x + big(i)
        end
        result
    end

    function f(n, λ)
        x = ((2n + 1) ^ 2 / (4 * (n + λ))) - n
        (1 / (n + λ) - 4 / (2n + 1)) * Pochhammer(x, n - 1) / factorial(n)
    end

    real(4 + sum(f(big(k), λ) for k in 1:maxiter))
end

function print_diff(λ;Z=512, maxiter=1024)
    setprecision(Z)
    str_cpi = "$(compute_pi(λ;Z=Z, maxiter=maxiter))"
    str_pi = "$(BigFloat(π))"

    println("   j | 1 2")
    println("-----------")
    for j in 1:length(min(str_cpi, str_pi))
        flag = str_cpi[j] != str_pi[j] ? "*" : ""
        head = @sprintf("%4d", j)
        println("$(head) | $(str_cpi[j]) $(str_pi[j]) $flag")
    end
end

function diff_maxiter(;Z=2048, maxiter=128, iterK=6)
    R = 10:10:300
    
    f = plot(size=(500, 300), dpi=150)

    for k in 1:iterK
        mik = maxiter * k
        diff = setprecision(Z) do
            K = zeros(BigFloat, length(R))
            for (ir, λ) in enumerate(R)
                pi_λ = compute_pi(λ; Z=Z, maxiter=mik)
                K[ir] = pi_λ
            end
            abs.(BigFloat(π) .- K) .+ BigFloat(1e-128)
        end
        plot!(f, R, diff, label="loss(Z=$Z, mi=$mik)", yscale=:log10)
    end

    savefig(f, "figures/loss_diff_mi_Z$(Z)_K$(iterK).png")
end