using EconLH, Test

uf = EconLH;

dsum(pBeta, xV) = sum((pBeta .^ (0 : (length(xV)-1))) .* xV);


function discounted_sum_test()
    @testset "Discounted sum" begin
        pBeta = 0.9;
        xV = 3.0
        @test discounted_sum(pBeta, xV) == xV

        xV = collect(range(-3.0, 2.0, length = 5));
        @test isapprox(dsum(pBeta, xV), discounted_sum(pBeta, xV), atol = 1e-7)
    end
end


function constant_growth_test()
    @testset "Constant growth" begin
        T = 7;
        x1 = 3.0;
        for gFactor in (0.8, 1.0, 1.3)
            v = uf.const_growth_vector(x1, gFactor, T);
            @test v[1] == x1
            @test length(v) == T
            for j = 2 : T
                @test isapprox(v[j] / v[j-1], gFactor)
            end

            vSum = uf.const_growth_sum(gFactor, T);
            @test isa(vSum, typeof(x1))
            @test isapprox(vSum * x1, sum(v))
        end
    end
end


@testset "Helpers" begin
    discounted_sum_test();
    constant_growth_test();
end

# ------------