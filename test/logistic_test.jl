using EconLH, Test;

function logistic_test()
    @testset "Logistic" begin
        n = 15;
        xV = range(-2.0, 3.0, length = n);

        lb = -2.5;
        ub = 1.5;
        f0 = 3.0;
        yV = logistic(xV; lb, ub, f0);
        @test isa(yV, Vector{Float64})
        @test size(yV) == size(xV)
        @test all(yV .> lb)  &&  all(yV .< ub);
        @test all(diff(yV) .> 0.0);

        glf = GeneralizedLogistic{Float64}(; lb, ub, f0);
        y2V = logistic(glf, xV);
        @test yV ≈ y2V
    end
end


@testset "Logistic" begin
    logistic_test();
end

# ----------