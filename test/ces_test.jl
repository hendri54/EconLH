using EconLH.ProductionFunctionsLH
using Random, Test
# Tests that output is continuous around elasticity 1
function elast_close_to_one_test()
  n = 3;
  T = 5;
  xM = LinRange(0.5, 5.0, T) * LinRange(2.0, 0.1, n)';

  elastV = LinRange(0.95, 1.05, 50);
  y1V = zeros(length(elastV));
  for (j, substElast) in enumerate(elastV)
    fS = pf.make_test_ces(T, n, substElast);
    yV = output(fS, xM);
    y1V[j] = yV[1]
  end
  @test all(abs.(diff(y1V)) .< 0.01)

  @show y1V
end

@testset "CES" begin
  elast_close_to_one_test();  
end

# ------------------
