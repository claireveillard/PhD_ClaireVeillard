function [alpha,d13Ccalc,d13Cco2,Talpha] = romanekCalcCO2(T,d13Ccalc,d13Cco2)
         
        if isnan(T)
            alpha=1000*(((d13Ccalc+1000)/(d13Cco2+1000))-1);
            Talpha=(-alpha+11.98)/0.12;
        else
            alpha=11.98-(0.12*T);
             % If no d18Odol is specified, d18Odol=from alpha and d18Ow
            if isempty(d13Ccalc)
               d13Ccalc = (((alpha/1000)+1)*(d13Cco2+1000))-1000;
            end
             % If no d18Ow is specified, d18Ow=from alpha and d18Odol
            if isempty(d13Cco2)
               d13Cco2=1/(((alpha/1000)+1)/(d13Cclac+1000));
            end
        end
end