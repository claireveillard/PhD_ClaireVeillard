
clc
clear

%% PREAMBLE
%This script computes the change in clumped, oxygen and carbon isotopes
% during  addition of dolomite cement to an initial dolomite D0. 
%The user chooses the total volume of cement, the temperature of the
%cement, the oxygen isotope composition of the fluid from which the cement
%has precipitated and the carbon isotopic composition of the cement. 

%Section 1 describes the initial conditions. The user can change the
%initial condition in the matlab file called 'ic'.

%Section 2 is the iterative calculation of isotopic abundance.
%The outputs are stored in a matrix M of dimension 2. 
%The first dimension is the number of model iterations.
%The second dimension is the number of variables calculated in the model 
%(e.g.d18O, TD47 etc.). 

%Section 5 shows how to plot the results.  
%Subplot 1: Bulk dolomite clumped isotope composition as a function of its
%oxygen isotope composition during cementation.

%Subplot 2: Bulk dolomite clumped isotope composition as a function of the 
%apparent fluid composition during cementation.

%Subplot 3: Bulk dolomite clumped isotope composition as a function of its
%carbon isotope composition during cementation.

%% REQUIRED FUNCTIONS AND CLASS
addpath('/Users/cv14/Documents/MATLAB/functions')
%daviesD2T: conversion from D47 to temperature
%daviesT2D: conversion from temperature to D47
%plot_isolines_d18Ow: calibration of Matthew and Katz (1977) to plot 
%isolines of d18Ow on a T=f(d18Odol) plot
%plot_isolines_d18Odol: calibration of Matthew and Katz (1977) to plot 
%isoline of d18Odol on a T=f(d18Ow) plot

%% DEFINITION OF THE VARIABLES

% Bulk dolomite
%Add '0' to the variable name to describe the initial conditions
%d13Cdol: Carbon isotopic composition of dolomite (permill)
%d18Odol: Oxygen isotopic composition of dolomite (permill)
%D47dol: Clumped isotopic composition of dolomite (permill)
%V0: Initial dolomite volume (from 0 to 1)
%V:Total volume of dolomite: D0 + cement

% Dolomite cement
%d13Ccm: Carbon isotopic composition of dolomite cement (permill)
%d18Ocm: Oxygen isotopic composition of dolomite cement  (permill)
%D47cm: Clumped isotopic composition of dolomite cement (permill)
%Tcm: Temperature of dolomite cement (°C)
%D47cm: Clumped isotopic composition of dolomite cement (permill)

%Fluids
%d18Ow: Oxygen isotopic composition of fluid (permill)
%d18Owapp: Apparent oxygen isotopic composition of fluid (permill)

% *Fluid-rock interaction*
%n: Model iteration (i.e pore volume of fluid) ([ ])


%% SECTION 1: INITIAL CONDITIONS: DOLOMITE (D0) AND FLUID (W0)

%D0
[~,d18Odol0,~]=matthewDW(ic.T,'',ic.d18Ow);
d18Odol=d18Odol0;
D47dol=daviesT2D(ic.T);
d13Cdol=ic.d13Cdol;
V0=1; 
V=V0;
n=0;

%W0
d18Ow=ic.d18Ow;

%% DOLOMITE CEMENT
%*The user chooses the total volume of dolomite cement*
Vtot=input(['Enter the total volume of dolomite cement. e.g. Enter 5 if'...
    '\n there is 5 times more cement than initial dolomite \n']);

%*The user chooses the dolomite cement temperature of precipitation*
Tcm=input('Enter the temperature of dolomite cement. e.g. Enter 40\n');

%*The user chooses the fluid oxygen isotopic composition 
%from which the cement precipitates (in VSMOW)*
d18Owcm=input(['Enter d18O (SMOW) from which the cement precipitates.'...
    'e.g. Enter -2 \n']);

%*The user chooses the cement carbon isotopic composition* 

d13Ccm=input(['Enter d13C (VPDB) of dolomite cement'...
    'e.g. Enter 0 \n']);

%Compute the clumped isotope composition of the cement
[acm,d18Ocm,~]=matthewDW(Tcm,'',d18Owcm);
D47cm=daviesT2D(Tcm);

%% ITERATIVE CALCULATION OF ISOTOPIC ABUNDANCES
listd18Odol=[];
listD47=[];
listVcm=[];
listV=[];
listd18Owapp=[];
listTcm=[];
M=[];

for Vcm=0:1:Vtot*V0
    n=n+1;
    
    V=Vcm+V;
    
    %Mixing of oxygen isotopes
    D47dol=((V*D47dol)+(Vcm*D47cm))/(V+Vcm);
    d18Odol=((V*d18Odol)+(Vcm*d18Ocm))/(V+Vcm);
    
    %Mixing of carbon isotopes
    d13Cdol=((V*d13Cdol)+(Vcm*d13Ccm))/(V+Vcm); 

    %FLUID
    [~,~,d18Owapp] = matthewDW(daviesD2T(D47dol),d18Odol,'');

    %RESULT
    listd18Odol(n,1)=d18Odol-30.92;
    listD47(n,1)=D47dol;
    listVcm(n,1)=Vcm;
    listV(n,1)=V;
    listd18Owapp(n,1)=d18Owapp;
    listTcm(n,1)=Tcm;
    listd13C(n,1)=d13Cdol;

end
M(:,:)=[listd18Odol,listD47,listd18Owapp,listV,listTcm,listd13C];    

%% PLOTS

figure('Renderer', 'painters', 'Position', [10 10 900 600])

%Subplot 1: Bulk dolomite clumped isotope composition as a function of its
%oxygen isotope composition during cementation.

subplot(1,3,1)

plot_isolines_d18Ow([d18Odol0-30.92 M(end,1)],[ic.T Tcm])
hold on;

plot(M(:,1),arrayfun(@daviesD2T,M(:,2)),'-d','Color','k','MarkerSize',...
        10,'MarkerFaceColor','none','LineWidth',2); hold on;
 
plot(d18Odol0-30.92,ic.T,'MarkerSize',8,'Marker','d',...
    'MarkerFaceColor','w','Color','k'); hold on;
text(d18Odol0-30.92,ic.T,'D0','FontSize',20,'FontWeight','bold',...
    'VerticalAlignment','bottom')

set(findobj(gcf,'type','axes'),'FontSize',15,'FontWeight','Normal',...
    'YDir','reverse','LineWidth', 2,'XDir','reverse',...
    'TickLength',[0.01, 0.01],'FontSize', 20)
xlabel('\delta^{18}O_{dol}','Fontsize',20);
ylabel('T(\Delta_{47 dol}) (ºC)','Fontsize',20);
box on


%Subplot 2: Bulk dolomite clumped isotope composition as a function of the 
%apparent fluid composition during cementation.

subplot(1,3,2)

plot_isolines_d18Odol([ic.d18Ow M(end,3)],[ic.T Tcm]);

plot(M(:,3),arrayfun(@daviesD2T,M(:,2)),'-d','Color','k',...
    'MarkerSize',10,'MarkerFaceColor','none','LineWidth',2); hold on;

plot(ic.d18Ow,ic.T,'MarkerSize',8,'Marker','d','MarkerFaceColor','w',...
    'Color','k'); hold on;
text(ic.d18Ow,ic.T,'D0','FontSize',20,'FontWeight','bold',...
    'VerticalAlignment','bottom')

xlabel('\delta^{18}O_{w (app)}','Fontsize',20);
ylabel('T(\Delta_{47 dol}) (ºC)','Fontsize',20);
set(findobj(gcf,'type','axes'),'FontSize',15,'FontWeight','Normal',...
    'YDir','reverse','LineWidth', 2,'XDir','reverse',...
    'TickLength',[0.01, 0.01],'FontSize', 20);
box on

%Subplot 3: Bulk dolomite clumped isotope composition as a function of its
%carbon isotope composition during cementation.

subplot(1,3,3)

plot(M(:,6),arrayfun(@daviesD2T,M(:,2)),'-d','Color','k','MarkerSize',...
        10,'MarkerFaceColor','none','LineWidth',2); hold on;

plot(ic.d13Cdol,ic.T,'MarkerSize',8,'Marker','d',...
    'MarkerFaceColor','w','Color','k');
text(ic.d13Cdol,ic.T,'D0','FontSize',20,'FontWeight','bold',...
    'VerticalAlignment','bottom')

xlabel('\delta^{13}C_{dol}','Fontsize',20);
ylabel('T(\Delta_{47 dol}) (ºC)','Fontsize',20);
set(findobj(gcf,'type','axes'),'FontSize',15,'FontWeight','Normal',...
    'YDir','reverse','LineWidth', 2,'XDir','reverse',...
    'TickLength',[0.01, 0.01],'FontSize', 20);

if Tcm<ic.T
    ylim([Tcm ic.T])
else
    ylim([ic.T Tcm])
end
legend off
box on

