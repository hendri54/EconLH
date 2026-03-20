"""
	$(SIGNATURES)

Utility from consumption.
"""
function utility(uS :: CRRA, cM :: AbstractArray{F1}) where F1
    if uS.dbg
        @assert all(x -> x .> zero(F1), cM)  "Negative consumption"
    end

    if uS.sigma == 1.0
        utilM = log.(cM);
    else
        utilM = utility_crra(cM, 1.0 - uS.sigma);
    end
    return utilM
end

# utility_log(cM) = log.(cM);
utility_crra(cM, oneMinusSigma) = (cM .^ oneMinusSigma) ./ oneMinusSigma .- 1.0;


function utility(uS :: CRRA, c :: T1) where T1 <: AbstractFloat
    if uS.sigma == 1.0
        utilM = log(c);
    else
        # oneMinusSigma = 1.0 - uS.sigma;
        utilM = utility_crra(c, 1.0 - uS.sigma);
    end
    return utilM
end


"""
	$(SIGNATURES)

Marginal utility.
"""
function marginal_utility(uS :: CRRA, cM :: AbstractArray{T1}) where
    T1 <: AbstractFloat

    if uS.dbg
        @assert all(x -> x > 0.0, cM)  "Negative consumption"
    end

    if uS.sigma == 1.0
        muM = one(T1) ./ cM;
    else
        muM = cM .^ (-uS.sigma);
    end
    return muM
end

function marginal_utility(uS :: CRRA, c :: T1) where T1 <: AbstractFloat
    if uS.sigma == 1.0
        muM = one(T1) / c;
    else
        muM = c ^ (-uS.sigma);
    end
    return muM
end


"""
    $(SIGNATURES)

Euler equation deviation.
"""
function euler_deviation(uS :: UtilityOneArg, beta :: T1, R :: T1,
    cV :: Vector{T1}) where T1 <: AbstractFloat

    muV = marginal_utility(uS, cV);
    T = length(muV);
    return muV[2 : T] .* beta .* R ./ muV[1 : (T-1)] .- 1.0;
end


"""
	$(SIGNATURES)

Consumption growth factor. Constant for CRRA. Otherwise this would have to depend on consumption levels.
"""
function cons_growth(uS :: CRRA,  betaR :: T1) where
    T1 <: AbstractFloat

    return betaR ^ (1.0 / uS.sigma)
end


# """
# 	$(SIGNATURES)

# Growth rate of utility.
# Input is consumption growth factor.
# """
# util_growth(uS :: CRRA, cGrowthFactor :: T1) where T1 <: AbstractFloat = 
#     cGrowthFactor ^ (one(T1) - uS.sigma);


"""
    $(SIGNATURES)

Present value of consumption as multiple of c1.
Periods 1...T. First not discounted.
"""
function pv_consumption(uS :: CRRA,
    beta :: T1, R :: T1, T :: T2) where
    {T1 <: AbstractFloat, T2 <: Integer}

    g = cons_growth(uS, beta * R);
    pv = const_growth_sum(g/R, T);
    return pv
    # return ((g/R) ^ T - one(T1)) / (g/R - one(T1))
end


"""
	$(SIGNATURES)

Consumption at age 1 for given lifetime income.
"""
function cons_age1(uS :: CRRA, beta :: T1, R :: T1, T :: T2,
    ltIncome) where {T1 <: AbstractFloat, T2 <: Integer}

    c1 = ltIncome ./ pv_consumption(uS, beta, R, T);
    return c1
end


"""
    $(SIGNATURES)

Consumption path.
"""
function cons_path(uS :: CRRA, beta :: T1, R :: T1, T :: T2,
    ltIncome :: T1) where {T1 <: AbstractFloat, T2 <: Integer}

    c1 = cons_age1(uS, beta, R, T, ltIncome);
    g = cons_growth(uS, beta .* R);
    cV = const_growth_vector(c1, g, T);
    return cV
    # return c1 .* (g .^ (0 : (T-1)))
end


"""
    $(SIGNATURES)

Lifetime utility, given present value of income.

Note that the growth rate of utility is not constant (because of the -1 in the utility function). So lifetime utility has to be computed from the consumption or utility path.
"""
function lifetime_utility(uS :: CRRA, beta :: T1, R :: T1, T :: T2,
    ltIncome :: T1) where {T1 <: AbstractFloat, T2 <: Integer}

    if uS.dbg
        @assert ltIncome .> 0.0  "Negative income"
    end

    # cGrowthFactor = cons_growth(uS, beta * R);
    # c1 = cons_age1(uS, beta, R, T, ltIncome);

    # if uS.sigma == 1
    #   uFct = c -> log.(c);
    # else
    #   uFct = c -> utility_crra(c, 1.0 - uS.sigma);
    # end
    # util1 = uFct(c1);
    # # util1 = utility(uS, c1);

    # ltu = util1;
    # ct = c1;
    # betaFactor = one(T1);
    # for t = 2 : T
    #     ct *= cGrowthFactor;
    #     betaFactor *= beta;
    #     ltu += betaFactor * uFct(ct);
    #     # utility(uS, ct);
    # end
    # return ltu
    utilV = utility(uS, cons_path(uS, beta, R, T, ltIncome));
    return lifetime_utility(utilV, beta, T)
end


"""
	$(SIGNATURES)

Marginal utility of initial assets given present value of income.
"""
function mu_wealth(uS :: CRRA, beta :: T1, R :: T1, T :: T2,
    ltIncome) where {T1 <: AbstractFloat, T2 <: Integer}

    c1 = cons_age1(uS, beta, R, T, ltIncome);
    return marginal_utility(uS, c1)
end


"""
	$(SIGNATURES)

Lifetime utility for a given utility sequence.
"""
function lifetime_utility(utilV :: Vector{T1},  beta :: T1,  T :: Integer)  where T1 <: AbstractFloat
    return discounted_sum(beta, utilV);
end

# Same for a constant utility flow
function lifetime_utility(util :: T1, beta :: T1, T :: Integer) where T1 <: AbstractFloat

    return util * (beta ^ T - 1.0) / (beta - 1.0)
end

# -----------
