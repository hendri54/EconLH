using Distributions, Random, Test
using EconLH.ProductionFunctionsLH
using EconLH

function pv_test()
	@testset "present value" begin
		R = 1.05;
		T = 3;

		function f(x)
			return 0.7 .* x
		end
		pv = EconLH.present_value(f, R, T);
		@test length(pv) == 1
		@test isa(pv, Float64)
		@test pv ≈ f(1) + f(2) ./ R + f(3) ./ R^2
	end
end

@testset "econLH" begin
	pv_test();
	include("extreme_value_decision_test.jl");
	include("prod_fct_test.jl");
	include("ces_test.jl")
end

# ----------
