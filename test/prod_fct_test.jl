using Test
using EconLH.ProductionFunctionsLH

pf = ProductionFunctionsLH;

function prod_fct_test(fS)
    println(fS)
    n = n_inputs(fS);
    T = length(productivities(fS));
    xM = LinRange(0.5, 10.5, T) * LinRange(1.2, 0.4, n)';

    yV = output(fS, xM);
    @test size(yV) == (T,)
    @test all(yV .> 0.0)

    if constant_returns(fS)
        y2V = output(fS, xM .* 2.0);
        @test all(isapprox(yV .* 2.0, y2V, atol = 1e-4));
        for t = 1 : T
            xMin, xMax = extrema(xM[t,:]);
            @test yV[t] / productivities(fS)[t] > xMin
            @test yV[t] / productivities(fS)[t] < xMax
        end
    end

    # Marginal products
    mpM = mproducts(fS, xM);
    println(mpM)
    dx = 0.00001;
    for j = 1 : n
      x2M = copy(xM);
      x2M[:,j] = x2M[:,j] .+ dx;
      y2V = output(fS, x2M);
      derivV = (y2V .- yV) ./ dx;
      @test  all(abs.(derivV .- mpM[:,j]) ./ max.(0.1, mpM[:,j]) .< 1e-3)

      if !all(abs.(derivV .- mpM[:,j]) ./ max.(0.1, mpM[:,j]) .< 1e-3)
        println(derivV)
        println(mpM[:, j])
      end
    end
end


@testset "Production Functions" begin
    n = 3;
    T = 5;
    substElastV = (0.5, 0.99, 1.0, 1.01, 3.0);
    for substElast in substElastV
        fS = pf.make_test_ces(T, n, substElast);
        prod_fct_test(fS);
    end
end

# -------------