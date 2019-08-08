clear
clc

%% PREAMBLE
%This script computes the change in clumped, oxygen and carbon isotopes
% during  dolomite crystallization from a mixture of fluids.

%Section 1 describes the initial conditions for fluid A and fluid B.
%Three cases are tested: the relative carbon content in fluid A is 20 times
%less than in fluid B (Pc=0.05), the carbon content in fluid A is
%equal to the carbon content in fluid B and the relative carbon content 
%in fluid A is 20 times more than in fluid B (Pc=20). 

%Section 2 is the iterative calculation of isotopic abundance.
%The outputs are stored in a matrix M of dimension 3. 
%The first dimension is the number of mixtures
%The second dimension is the number of variables calculated in the model 
%(e.g.d18O, D47 etc.). 
%The third dimension is the number of Pc tested.

%Section 3 shows how to plot the results.  

%% REQUIRED FUNCTIONS AND CLASS
%write the path where the fucntions and classes are stored
addpath('/Users/cv14/Documents/MATLAB/functions')

%daviesD2T: conversion from D47 to temperature
%daviesT2D: conversion from temperature to D47
%mattewDW: calibration between temperature, oxygen isotopic composition
%(calibration of Matthew and Katz 1977)
%romanekCalcCO2: calibration between temperature, dissolved CO2 and calcite
%(calibration of Romanek 1977)

%% DEFINITION OF THE VARIABLES
%*Fluid A*
%d18Oa: oxygen isotopic composition of fluid A
%d13Ca: carbon isotopic composition of fluid A
%Ta: temperature of fluid A

%*Fluid B*
%d18Ob: oxygen isotopic composition of fluid A
%d13Cb: carbon isotopic composition of fluid A
%Tb: temperature of fluid A

%*Fluid mixture*
%Pc: the relative carbon content in fluid A relative to fluid B
%r: relative weight of fluid A divided by fluid B
%d18Owm: oxygen isotopic composition of the mixture
%d13Cwm: carbon isotopic composition of the mixture
%Twm: temperature of the fluid mixture

%*Dolomite
%D47a: clumped isotope composition of dolomite which have precipitated from
%100% of fluid A

%D47b: clumped isotope composition of dolomite which have precipitated from
%100% of fluid B


%% SECTION 1: INITIAL CONDITIONS- Fluid A and Fluid B
%Fluid A
d18Oa=0; 
d13Ca=-10;
Ta=15;
D47a=daviesT2D(Ta);

%Fluid B
d18Ob=4;
d13Cb=-5;
Tb=40;
D47b=daviesT2D(Tb);

%Select the values of PC
Pc_range=[0.05,1,20];

l=0;

%% SECTION 2: RUN THE MODEL FOR DIFFERENT CASES
for k=1:size(Pc_range,2)
    Pc=Pc_range(:,k);
    n=0;
    l=l+1;
    r=0;
    
    %% SECTION 3: ITERATIVE CALCULATION OF ISOTOPIC ABUNDANCES 
    while r<500
        r=r+0.001;
        n=n+1;
        
        %Mixing of oxygen isotopes in the fluid mixture
        d18Owm=((r*d18Oa) + d18Ob)/(1+r);

        %Mixing of carbon isotopes in the fluid mixture
        d13Cwm=((r*d13Ca)+(Pc*d13Cb))/(Pc+r);
        
        %Clumped isotopes of dolomite which precipitate from the mixture
        D47m=((r*D47a)+D47b)/(1+r);
        
        %Oxygen isotopes of dolomite which precipitate from the mixture
        [a,d18Odol,~]=matthewDW(daviesD2T(D47m),'',d18Owm);
        
        %Carbon isotopes of dolomite which precipitate from the mixture
        [ac,d13Cdol,~]=romanekCalcCO2(daviesD2T(D47m),'',d13Cwm);
        
        %Store the results
        listd18Om(n,1)=d18Owm;
        listd13Cm(n,1)=d13Cwm;
        listr(n,1)=r;
        listD47m(n,1)=D47m;
        listd18Odol(n,1)=d18Odol-30.92;
        listd13Cdol(n,1)=d13Cdol;
        
    end
   
 M(:,:,l)=[listd18Om,listd13Cm,listD47m,listr,listd18Odol,listd13Cdol];  
end

M=[M,arrayfun(@daviesD2T,M(:,3,:))];

%Find the rows of M corresponding to the values r=0.1, r=1 and r=10.
idx1=find(round(M(:,4,1),3)==0.1);
idx2=find(round(M(:,4,1),3)==1);
idx3=find(round(M(:,4,1),3)==10);

%% SECTION 4: PLOT THE RESULTS

%subplot(2,2,1): Plot the oxygen isotopic composition of the fluid mixture 
%subplot(2,2,2): Plot the oxygen isotopic composition of dolomite
%subplot(2,2,3): Plot the carbon isotopic composition of the fluid mixture 
%subplot(2,2,4): Plot the carbon isotopic composition of dolomite

col= [[0.8 0.5 0.2];[0 0 1];[1 0 1]];

%..........Plot the oxygen isotopic composition of the fluid mixture...
subplot(2,2,1)
ax1=gca;
plot(M(:,1,1),arrayfun(@daviesD2T,M(:,3,1)),'r'); hold on;
plot(d18Ob,Tb,'ko'); hold on; plot(d18Oa,Ta,'ko');
text(d18Ob+0.2,Tb+2,'B'); 
text(d18Oa-0.2,Ta-2,'A');

plot([M(idx1,1,1) M(idx2,1,1) M(idx3,1,1)],...
    [arrayfun(@daviesD2T,M(idx1,3,1)) arrayfun(@daviesD2T,M(idx2,3,1)) ...
    arrayfun(@daviesD2T,M(idx3,3,1))] ,'kx','Color','r')
text([M(idx1,1,1) M(idx2,1,1) M(idx3,1,1)],...
    [arrayfun(@daviesD2T,M(idx1,3,1)) arrayfun(@daviesD2T,M(idx2,3,1))...
    arrayfun(@daviesD2T,M(idx3,3,1))] ,...
    num2cell([M(idx1,4,1) M(idx2,4,1) M(idx3,4,1) ]),...
    'VerticalAlignment','bottom','Color','r')

%..........Plot the oxygen isotopic composition of dolomite...
subplot(2,2,2)
ax2=gca;
plot(M(:,5,1),arrayfun(@daviesD2T,M(:,3,1)),'r'); hold on;
plot(M(1,5,1),Tb,'ko'); hold on; plot(M(end,5,1),Ta,'ko');
text(M(1,5,1),Tb+2,'DB'); text(M(end,5,1),Ta-2,'DA');

plot([M(idx1,5,1) M(idx2,5,1) M(idx3,5,1)],...
    [arrayfun(@daviesD2T,M(idx1,3,1)) arrayfun(@daviesD2T,M(idx2,3,1)) ...
    arrayfun(@daviesD2T,M(idx3,3,1))] ,'kx','Color','r')
text([M(idx1,5,1) M(idx2,5,1) M(idx3,5,1)],...
    [arrayfun(@daviesD2T,M(idx1,3,1)) arrayfun(@daviesD2T,M(idx2,3,1)) ...
    arrayfun(@daviesD2T,M(idx3,3,1))] ,num2cell([M(idx1,4,1) M(idx2,4,1)...
    M(idx3,4,1) ]),'VerticalAlignment','bottom','Color','r')

%Plot the carbon isotopic composition of the fluid mixture and dolomite
lgd=[];

for kk=1:size(Pc_range,2)
    
    %..........Plot the carbon isotopic composition of the fluid mixture...
    subplot(2,2,3)
    ax3=gca;
    
    plot(M(:,2,kk),arrayfun(@daviesD2T,M(:,3,kk)),'Color',col(kk,:)); 
    hold on;
    plot(d13Cb,Tb,'ko'); hold on; plot(d13Ca,Ta,'ko');
    text(d13Cb,Tb+2,'B'); hold on; text(d13Ca,Ta-3,'A');
    
    plot([M(idx1,2,kk) M(idx2,2,kk) M(idx3,2,kk)],...
        [arrayfun(@daviesD2T,M(idx1,3,kk)) ...
        arrayfun(@daviesD2T,M(idx2,3,kk)) ...
        arrayfun(@daviesD2T,M(idx3,3,kk))] ,...
        'kx','MarkerFaceColor','k','Color',col(kk,:)); hold on;
    text([M(idx1,2,kk) M(idx2,2,kk) M(idx3,2,kk)],...
        [arrayfun(@daviesD2T,M(idx1,3,kk)) ...
        arrayfun(@daviesD2T,M(idx2,3,kk)) ...
        arrayfun(@daviesD2T,M(idx3,3,kk))] ,...
        num2cell([M(idx1,4,kk) M(idx2,4,kk) M(idx3,4,kk) ]),...
        'VerticalAlignment','bottom','Color',col(kk,:)); hold on;
   
    %..........Plot the carbon isotopic composition of dolomite............
    subplot(2,2,4)
    ax4=gca;
    
    p4=plot(M(:,6,kk),arrayfun(@daviesD2T,M(:,3,kk)),'Color',col(kk,:),...
    'DisplayName',sprintf('P_{C}=%g', Pc_range(:,kk))); hold on;
    plot([M(idx1,6,kk) M(idx2,6,kk) M(idx3,6,kk)],...
        [arrayfun(@daviesD2T,M(idx1,3,kk)) ...
        arrayfun(@daviesD2T,M(idx2,3,kk)) ...
        arrayfun(@daviesD2T,M(idx3,3,kk))] ,...
        'kx','MarkerFaceColor','k','Color',col(kk,:));
    text([M(idx1,6,kk) M(idx2,6,kk) M(idx3,6,kk)],...
        [arrayfun(@daviesD2T,M(idx1,3,kk)) ...
        arrayfun(@daviesD2T,M(idx2,3,kk)) ...
        arrayfun(@daviesD2T,M(idx3,3,kk))] ,...
        num2cell([M(idx1,4,kk) M(idx2,4,kk) M(idx3,4,kk) ]),...
        'VerticalAlignment','bottom','Color',col(kk,:)); hold on;
    
    plot(M(1,6,1),Tb,'ko'); hold on; plot(M(end,6,1),Ta,'ko');
    text(M(1,6,1),Tb+2,'DB'); text(M(end,6,1),Ta-3,'DA');
    lgd=[lgd p4];
end

legend(lgd)
ax1.LineWidth = 2; ax1.LineWidth = 2;ax1.YDir='reverse'; 
ax1.XDir='reverse';
ylim(ax1,[10 45]);xlim(ax1,[d18Oa-0.5 d18Ob+0.5])
xlabel(ax1,'\delta^{18}O_{wm}');ylabel(ax1,'T_{wm}')


ax2.LineWidth = 2; ax2.LineWidth = 2;ax2.YDir='reverse'; ax2.XDir='reverse';
ylim(ax2,[10 45]);
xlabel(ax2,'\delta^{18}O_{dol}');ylabel(ax2,'T(\Delta_{47})_{dol}')
 
ax3.LineWidth = 2; ax3.LineWidth = 2;ax3.YDir='reverse'; ax3.XDir='reverse';
ylim(ax3,[10 45]);xlim(ax3,[d13Ca-0.5 d13Cb+0.5])
xlabel(ax3,'\delta^{13}C_{CO2}');ylabel(ax3,'T_{wm}')

ax4.LineWidth = 2; ax4.LineWidth = 2;ax4.YDir='reverse'; ax4.XDir='reverse';
ylim(ax4,[10 45]);
xlabel(ax4,'\delta^{13}C_{dol}');ylabel(ax4,'T(\Delta_{47})_{dol}')

