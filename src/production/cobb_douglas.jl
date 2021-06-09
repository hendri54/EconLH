"""
	$(SIGNATURES)

Output of a Cobb-Douglas. Returns a vector of length 1 for xM with 1 row.
"""
function output_cobb_douglas(xM :: AbstractMatrix{F1}, 
    expV :: AbstractVector{F1}) where F1 <: Real

    T, n = size(xM);
    outputV = ones(F1, T);
    for j = 1 : n
        @views outputV .*= (xM[:,j] .^ expV[j]);
    end
    return outputV
end

# Returns scalar output.
function output_cobb_douglas(xV :: AbstractVector{F1}, 
    expV :: AbstractVector{F1}) where F1 <: Real

    n = length(xV);
    output = one(F1);
    for j = 1 : n
        output *= xV[j] ^ expV[j];
    end
    return output
end


function mproducts_cobb_douglas(xM :: AbstractMatrix{F1}, 
    expV :: AbstractVector{F1}) where F1 <: Real

    yV = output_cobb_douglas(xM, expV);
    mpM = (yV * expV') ./ xM;
    return mpM
end

# -----------