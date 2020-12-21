function output_cobb_douglas(xM :: AbstractMatrix{F1}, 
    expV :: AbstractVector{F1}) where F1 <: Real

    T, n = size(xM);
    outputV = ones(F1, T);
    for j = 1 : n
        outputV .*= (xM[:,j] .^ expV[j]);
    end
    return outputV
end


function mproducts_cobb_douglas(xM :: AbstractMatrix{F1}, 
    expV :: AbstractVector{F1}) where F1 <: Real

    yV = output_cobb_douglas(xM, expV);
    mpM = (yV * expV') ./ xM;
    return mpM
end

# -----------