export gini;

"""
	$(SIGNATURES)

Gini coefficient with weighted data. For discrete probability distribution.
Based on https://en.wikipedia.org/wiki/Gini_coefficient#Definition
"""
function gini(xV :: AbstractVector{T}, wtV :: AbstractVector{T}) where T
    @assert all(x -> x >= 0, wtV);
    @assert all(x -> x >= 0, xV);
    
    totalWt = sum(wtV);
    idxV = sortperm(xV);
    cumValue = zero(T);
    numer = zero(T);
    for j in idxV
        wt = wtV[j] / totalWt;
        x = xV[j];
        newProd = wt * x;
        numer += wt * (2 * cumValue + newProd);
        cumValue += newProd;
    end
    g = one(T) - numer / cumValue;
    return g
end

# -----------