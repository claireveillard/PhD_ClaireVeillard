function [alpha,d18Odol,d18Ow,Talpha] = matthewDW(T,d18Odol,d18Ow)
        
        if isnan(T)
           alpha = (d18Odol+1000)/(d18Ow+1000);
           Talpha=(sqrt(1/((1000*log(alpha)+3.24)/3060000)))-273.15;
        else
            alpha=exp(((3060000/power(T+273,2))-3.24)/1000);
            % If no d18Odol is specified, d18Odol=from alpha and d18Ow
            if isempty(d18Odol)
           d18Odol = alpha*(d18Ow+1000)-1000;
            end
            % If no d18Ow is specified, d18Ow=from alpha and d18Odol
            if isempty(d18Ow)
           d18Ow=(d18Odol+1000-(alpha*1000))/alpha;
            end
        end

end
