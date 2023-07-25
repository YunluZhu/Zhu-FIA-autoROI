function dFoF_KSDensity = extract_dFF_ksDensity(response)
    numROIs=size(response,2);
    mu=[];
    [smoothDist,x] = ksdensity(response(:));
    for thisROI = 1:numROIs
      traceCell = response(:,thisROI);
      [smoothDist,x] = ksdensity(traceCell);
      [valuePeak,indPeak]=max(smoothDist);
      mu(thisROI) = mean(x(1:indPeak));
    end   
    FminusMu = bsxfun(@minus,response, mu);
    dFoF_KSDensity = {bsxfun(@rdivide,FminusMu, mu)};
end