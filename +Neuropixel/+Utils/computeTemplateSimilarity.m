function [simScore, bestLag] = computeTemplateSimilarity(W, U)
    % covariance matrix between all templates
    WtW = Neuropixel.Utils.getMeWtW_nomex(single(W), single(U));
    
    nt0 = size(W, 1);
    lags = nt0:-1:-nt0;       
    
    % the similarity score between templates is simply the correlation,
    % taken as the max over several consecutive time delays
    [simScore, lagInd] = max(WtW, [], 3);
    simScore = min(max(simScore, 0), 1);
    bestLag = lags(lagInd);
    
    % when not similar at all, lag estimates are useless
    bestLag(simScore < 0.1) = 0;
end