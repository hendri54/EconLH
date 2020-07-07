using Distributions, Random, Test

## Helper: check a decision by simulation
function check_by_sim(d :: ExtremeValueDecision{T},  value_ixM :: Array{T,2},
    prob_ixM :: Array{T,2},  eVal_iV :: Array{T,1})    where T <: AbstractFloat

    rng = MersenneTwister(123);
    nTypes, n = size(value_ixM);
    success = true;

    for iType = 1 : nTypes
        success1 = check_one_by_sim(d, vec(value_ixM[iType,:]),  
            vec(prob_ixM[iType,:]),  eVal_iV[iType],  rng);
        success = success && success1;
    end

    return success
end


## Check one type by simulation
# Returns a Bool (success or failure)
function check_one_by_sim(d :: ExtremeValueDecision{F},  valueV :: Vector{F},
    probV :: Vector{F},  eVal :: F,  rng :: AbstractRNG)  where  F <: AbstractFloat

    # Random vars for all alternatives
    n = length(valueV);
    nSim = Int64(5e5);
    simProbV, simVal = sim_extreme_value(d, valueV, nSim, rng);

    success = true;
    if !isapprox(simProbV, probV, atol = 5e-3)
       println("True probs:  ", probV);
       println("Sim probs:   ", simProbV);
       success = false;
    end

    # Check expected value
    if abs(simVal ./ eVal - 1.0) > 5e-3
       @warn "Discrepancy in values (demeaned = $(demeaned(d))): $simVal,  $eVal"
       success = false;
    end

    return success
end

function sim_extreme_value(d :: ExtremeValueDecision{F},  valueV :: Vector{F},
    nSim :: Integer, rng :: AbstractRNG) where F

    # Random vars for all alternatives
    n = length(valueV);
    maxV = zeros(F, nSim);
    countV = zeros(Int, n);
    for i1 = 1 : nSim
        randV = draw_gumbel_shocks(rng, pref_scale(d), n; demeaned = demeaned(d));
        # Shocks must be added, not subtracted.
        maxV[i1], maxIdx = findmax(valueV .+ randV);
        countV[maxIdx] += 1;
    end
    simProbV = countV ./ nSim;
    simVal = mean(maxV);
    return simProbV, simVal
end


function one_option_test()
    @testset "One option" begin
        prefScale = 1.4;
        nTypes = 3;
        value_iV = collect(range(32.0, 48.0, length = nTypes));

        # Demeaned
        d = ExtremeValueDecision(prefScale, true, true);
        prob_iV, eVal_iV = extreme_value_decision(d, value_iV);
        @test prob_iV == ones(Float64, nTypes)
        @test isapprox(eVal_iV, value_iV, rtol = 0.02)

        # Not demeaned
        d2 = ExtremeValueDecision(prefScale, false, true);
        prob_iV, eVal2_iV = extreme_value_decision(d2, value_iV);
        @test prob_iV == ones(Float64, nTypes)
        # If Gumbel shocks are subtracted (they have positive mean), demeaned expected values should be higher
        @test isapprox(eVal2_iV .- prefScale * EconLH.EulerConst, eVal_iV)
    end
end


function check_by_sim_test()
    @testset "Check by simulation" begin
        for demeaned in [true, false]
            n = 4;
            prefScale = 1.4;
            d = ExtremeValueDecision(prefScale, demeaned, true);

            # Row vector of values (as matrix)
            value_xV = collect(range(50.0, 55.0, length = n)');

            for nTypes in [1, 3]
                if nTypes == 1
                    value_ixM = value_xV;
                else
                    value_ixM = collect(range(1.0, 0.9, length = nTypes)) *
                        collect(range(50.0, 55.0, length = n)');
                end

                prob_ixM, eVal_iV = EconLH.extreme_value_decision(d, value_ixM);
                # This makes sure that there is something meaningful to simulate
                @test all(prob_ixM .< 0.95)
                @test check_by_sim(d, value_ixM, prob_ixM, eVal_iV);
            end
        end
    end
end

function draw_gumbel_test()
    @testset "Draw gumbel" begin
        rng = MersenneTwister(123);
        prefScale = 0.4;
        nSim = Int(1e3);

        shockV = draw_gumbel_shocks(rng, prefScale, nSim; demeaned = true);
        @test isapprox(mean(shockV), 0.0, atol = 0.05)
        @test size(shockV) == (nSim,)

        shockM = draw_gumbel_shocks(rng, prefScale, (nSim, 2); demeaned = false);
        @test isapprox(mean(shockM), prefScale * EconLH.EulerConst, atol = 0.05)
        @test size(shockM) == (nSim, 2)
    end
end

@testset "ExtremeValueDecision" begin
    one_option_test();
    check_by_sim_test();
    draw_gumbel_test();
end

# ---------------
