@testset "CRRA" begin
   println("CRRA_test")
   crraS = CRRA(2.0, true);
   cM = collect(range(2.0, 3.0, length = 4)) * collect(range(1.5, 2.5, length = 5))';
   uM = utility(crraS, cM);
   muM = marginal_utility(crraS, cM);
   @test size(cM) == size(uM)
   @test all(muM .> 0)

   dc = 1e-4;
   u2M = utility(crraS, cM .+ dc);
   mu2M = (u2M .- uM) ./ dc;
   @test muM ≈ mu2M  atol = 1e-5
end


@testset "PV consumption" begin
    crraS = CRRA(2.0, true)
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
    discV = (1/R) .^ (0 : (T-1));
    pv = sum(cV .* discV);
    @test pv ≈ ltIncome

    devV = euler_deviation(crraS, beta, R, cV);
    @test all(abs.(devV) .< 1e-4)

    utilV = utility(crraS, cV);
    betaV = beta .^ (0 : (T-1));
    ltUtil = sum(utilV .* betaV);
    @test ltUtil ≈ lifetime_utility(crraS, beta, R, T, ltIncome)
end


@testset "Euler" begin
    crraS = CRRA(2.0, true)
    R = 1.07;
    beta = 0.97;
    T = 5;
    g = cons_growth(crraS, beta .* R);
    devV = euler_deviation(crraS, beta, R, g .^ (1:T));
    @test all(abs.(devV) .< 1e-4)
end
