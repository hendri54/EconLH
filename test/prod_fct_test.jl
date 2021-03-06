using Test
using EconLH.ProductionFunctionsLH

pf = ProductionFunctionsLH;

function prod_fct_test(fS)
    println(fS)
    n = n_inputs(fS);
    T = length(productivities(fS));
    if T == 1
        xM = LinRange(1.2, 0.4, n)';
    else
        xM = LinRange(0.5, 10.5, T) * LinRange(1.2, 0.4, n)';
    end

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


function output_direct_test(substElast)
    @testset "Direct output computation" begin
        n = 3;
        T = 4;
        alphaV = (1 : n) ./ n;
        xM = (1:T) .+ (1:n)' .+ 0.5;
        yV = ces_output(xM, alphaV, substElast);
        @test all(yV .> 0.0);
        @test size(yV) == (T,);

        # One set of inputs
        y = ces_output(xM[2,:], alphaV, substElast);
        @test isapprox(y, yV[2]);
    end
end


@testset "Production Functions" begin
    n = 3;
    substElastV = (0.5, 0.99, 1.0, 1.01, 3.0);
    TV = (1, 5);
    for T in TV
        for substElast in substElastV
            fS = pf.make_test_ces(T, n, substElast);
            prod_fct_test(fS);
            output_direct_test(substElast);
        end
    end

end

# -------------