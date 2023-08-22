% https://github.com/amitani/baseline_kde/blob/master/mode_kde.m

function [dFoF_KSDensity, dF_KSDensity, mu] = extract_dFF_ksDensity_advanced(response)
    numROIs=size(response,2);
    mu=[];
    % [smoothDist,x] = ksdensity(response(:));
    for thisROI = 1:numROIs
      traceCell = response(:,thisROI);
      [f,xi] = ksdensity(traceCell,'npoints',200);
        [~,ii] = max(f);
        ii_1 = max(ii-1,1);
        ii_2 = min(ii+1,length(f));
        if(ii_2-ii_1 == 2) % if two are similarly high compared to the third,
                           % let's just focus on the interval between two.
            if(f(ii_2)>f(ii_1))
                if(f(ii)-f(ii_2)<f(ii_2)-f(ii_1))
                    ii_1 = ii;
                end
            else
                if(f(ii)-f(ii_1)<f(ii_1)-f(ii_2))
                    ii_2 = ii;
                end
            end
        end

      xx = linspace(xi(ii_1),xi(ii_2),201); % x100 interpolation
      [f, xi]=ksdensity(traceCell,xx);
      [~, ii] = max(f);
      mu(thisROI) = xi(ii);

    end   
    FminusMu = bsxfun(@minus,response, mu);
    dF_KSDensity = {FminusMu};
    dFoF_KSDensity = {bsxfun(@rdivide,FminusMu, mu)};
end