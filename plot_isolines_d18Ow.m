function plot_isolines_d18Ow(a,b)
        
        xaxis=[-20:1:21]+30.92; %%ROCK PDB
        zaxis=[-20:1:21]; %FLUID SMOW
        A = nan(1,size(zaxis,2));
    
        for i=1:2:size(zaxis,2)
            z(1:size(zaxis,2)) = zaxis(i);
            [~,~,~,Talpha]  = arrayfun(@matthewDW,A,xaxis,z);
            Talpha(imag(Talpha) ~= 0) = NaN;    
            
            plot(xaxis-30.92,Talpha,'--b','LineWidth',0.001); hold on;
                
            if a(1)<a(2) 
                    xlim(a)
            end
                
            if a(1)>a(2)
                xlim([a(2) a(1)])
            end
            
            if b(1)<b(2)
               ylim(b)
            end
                
            if b(1)>b(2)
                ylim([b(2) b(1)])
            end
                
            if a(1)==a(2)
                xlim([abs(a(1))-1 abs(a(2))+1])
            end
                
            if b(1)==b(2)
                ylim([abs(b(1))-10 abs(b(2))+10])
            end
        end
end 

