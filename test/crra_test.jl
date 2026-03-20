function crra_test()
    @testset "CRRA" begin
        println("CRRA_test")
        for pSigma in [1.0, 2.0]
            crraS = CRRA(pSigma, true);
            cM = collect(range(2.0, 3.0, length = 4)) * collect(range(1.5, 2.5, length = 5))';
            uM = utility(crraS, cM);
            u5 = utility(crraS, cM[5]);
            @test u5 ≈ uM[5]
            muM = marginal_utility(crraS, cM);
            @test size(cM) == size(uM)
            @test all(muM .> 0)

            dc = 1e-5;
            u2M = utility(crraS, cM .+ dc);
            mu2M = (u2M .- uM) ./ dc;
            @test muM ≈ mu2M  atol = 1e-5
        end
    end
end

function pv_test()
    @testset "PV consumption" begin
        for pSigma in [1.0, 2.0]
            crraS = CRRA(pSigma, true)
            R = 1.07;
            beta = 0.97;
            betaR = beta * R;
            T = 3;
            pv = pv_consumption(crraS, beta, R, T);
            g = cons_growth(crraS, betaR);
            # println("Consumption growth $g")
            pvTrue = 1 + g/R + (g/R)^2;
            @test pv ≈ pvTrue

            ltIncome = 17.3;
            cV = cons_path(crraS, beta, R, T, ltIncome);
            @test length(cV) == T  &&  all(cV .> 0.0)
            cGrowthV = cV[2 : T] ./ cV[1 : (T-1)];
            @test all(abs.(cGrowthV .- g) .< 1e-4)
            discV = (1/R) .^ (0 : (T-1));
            pv = sum(cV .* discV);
            @test pv ≈ ltIncome

            devV = euler_deviation(crraS, beta, R, cV);
            @test all(abs.(devV) .< 1e-4)

            utilV = utility(crraS, cV);
            # uGrowth = util_growth(crraS, g);
            # uGrowthV = utilV[2 : T] ./ utilV[1 : (T-1)];
            # @test all(abs.(uGrowthV .- uGrowth) .< 1e-4)
            betaV = beta .^ (0 : (T-1));
            ltUtil = sum(utilV .* betaV);
            @test ltUtil ≈ lifetime_utility(crraS, beta, R, T, ltIncome)

            # Two ways of computing lifetime utility from utility
            util0 = -1.23;
            u1 = lifetime_utility(util0, beta, T);
            u2 = lifetime_utility(fill(util0, T), beta, T);
            @test u1 ≈ u2

            # Marginal utility of wealth
            ltUtil ≈ lifetime_utility(crraS, beta, R, T, ltIncome);
            dy = 1e-6;
            ltUtil2 = lifetime_utility(crraS, beta, R, T, ltIncome + dy);
            mu = (ltUtil2 - ltUtil) / dy;
            mu2 = EconLH.mu_wealth(crraS, beta, R, T, ltIncome);
            @test isapprox(mu, mu2, atol = 1e-5)
        end
    end
end


function euler_test()
    @testset "Euler" begin
        for pSigma in [1.0, 2.0]
            crraS = CRRA(pSigma, true)
            R = 1.07;
            beta = 0.97;
            T = 5;
            g = cons_growth(crraS, beta .* R);
            devV = euler_deviation(crraS, beta, R, g .^ (1:T));
            @test all(abs.(devV) .< 1e-4)
        end
    end
end

@testset "CRRA" begin
    crra_test()
    pv_test()
    euler_test()
end

# --------------