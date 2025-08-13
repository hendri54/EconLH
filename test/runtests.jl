using Distributions, Random, Test
using EconLH


function gini_test()
	@testset "Gini" begin
		rng = MersenneTwister(43);
		n = 50_000;
		xV = rand(rng, n);
		wtV = 1.0 .+ rand(rng, n);
		g = gini(xV, wtV);
		# Gini for uniform distribution is about 1/3
		@test isapprox(g, 0.333; atol = 0.005);

		g1 = gini(10 .* xV, wtV);
		@test isapprox(g, g1; atol = 0.005);

		g2 = gini(1000 .+ xV, wtV);
		@test g2 < 0.05;

		xV[1] = 10_000.0;
		g3 = gini(xV, wtV);
		@test g3 > 0.5;
	end
end


function pv_test()
	@testset "present value" begin
		R = 1.05;
		T = 3;

		function f(x)
			return 0.7 .* x
		end
		pv = EconLH.present_value(f, R, T);
		@test isa(pv, Float64)
		@test pv ≈ f(1) + f(2) ./ R + f(3) ./ R^2
	end
end

@testset "econLH" begin
	gini_test();
	pv_test();
	include("extreme_value_decision_test.jl");
	include("prod_fct_test.jl");
	include("ces_test.jl");
	include("logistic_test.jl");
end

# ----------
