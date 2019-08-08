classdef cst
   properties (Constant)
     
      %concentration elements in water
       cMgf=1284;
       cCaf=413; %http://delloyd.50megs.com/moreinfo/seawater.html#Top
       FWf = 18.01528; %water
       cOf=889000;
       cCf=200;
       cSrf=480000;
       Rhof=1;
       
      %concentration elements in pure cacite
        cMgcal=0; %https://deepblue.lib.umich.edu/bitstream/handle/2027.42/30096/0000468.pdf?sequence=1&isAllowed=y
        cCacal=400000; %Banner
        
      %concentration elements in HMg cacite
        cMgHMC=50000; %https://deepblue.lib.umich.edu/bitstream/handle/2027.42/30096/0000468.pdf?sequence=1&isAllowed=y
        cCaHMC=200000;

      %concentration Mg and Ca dolomite MP
        cMgdol=91*1000;
        cCadol=204*1000;

      %concentration elements stoichio dolomite
        cMgdol5050=132000;
        cCadol5050=217000;
        cOdol=480000; %Banner and Hanson 1990
        cCdol=120000;
        cSrdol=10000; %Banner and Hanson 1990
        
        %Dolomite rock properties 
        Rhodol=2.85;
        FWs = 184.40; 

        %Calcite rock properties
        Rhocal=2.71;
        
   end
end