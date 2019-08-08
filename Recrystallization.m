clc
clear
%%.......Author: Claire Veillard..............%%
%%......Date 25/06/2019.......................%%

%% PREAMBLE
% See Veillard, C. M., John, C. M., Krevor, S., & Najorka, J. (2019).
%Rock-buffered recrystallization of Marion Plateau dolomites at 
%low temperature evidenced by clumped isotope thermometry and 
%X-ray diffraction analysis. 
%Geochimica et Cosmochimica Acta, 252, 190-212.

%This script computes the change in clumped, oxygen and carbon isotopes
% during  dolomite recrystallization via dissolution/re-precipitation. 
%The user chooses if recrystallization is complete or incomplete. 
%If recrystallization is  incomplete, the user chooses the percentage 
%of dolomite which recrystallizes at each iteration (n).

%Section 1 describes the initial conditions. The user can change the
%initial condition in the matlab file called 'ic'.

%Section 2 is the iterative calculation of isotopic abundance.
%The outputs are stored in a matrix M of dimension 3. 
%The first dimension is the number of model iterations - or pore volumes (n).
%The second dimension is the number of variables calculated in the model 
%(e.g.d18O, TD47 etc.). The third dimension is the number of cases (e.g.
%recrystallization happens at Tsys=20ºC, Tsys=30ºC etc.).

%Section 3 describes how to find isovalues of N in M. The user can choose
%the values of N in the array 'chooseN'.

%Section 4 describes the best fit functions between T47 and d18Odol 
%for the seleted value of N. The Matlab Optimization Toolbox is required. 

%Section 5 shows how to plot the results.  
%Figure 1: Plot of T=f(d18O). 
%Figure 2: Extent of the reaction, IEQ vs cumulative fluid/rock ratio (N). 
%Figure 3: Model parameters vs cumulative water/rock ratio (N).  


%% REQUIRED FUNCTIONS AND CLASS
%write the path where the fucntions and classes are stored
addpath('/Users/cv14/Documents/MATLAB/functions')

%daviesD2T: conversion from D47 to temperature
%daviesT2D: conversion from temperature to D47
%mattewDW: calibration between temperature, oxygen isotopic composition
%of water and dolomite from Matthew and Katz, 1997.
%plot_isolines_d18Ow: calibration of Matthew and Katz (1977) to plot 
%isolines of d18Ow on a T=f(d18Odol) plot
%vector color: array containing rgb triplets

%% DEFINITION OF THE VARIABLES
%Add '0' to the variable name to describe the initial conditions
%Add 'eq' to the variable name to describe the equilibrium

% *Elements concentration in dolomite and fluid*
%cCf: Concentration of carbon in fluid (ppm)
%cCdol: Concentration of carbon in dolomite (ppm)
%cOf: Concentration of oxygen in fluid (ppm)
%cOdol: Concentration of oxygen in dolomite (ppm)
%Rhodol: Dolomite density (kg/m^3)
%Rhof: Fluid density (kg/m^3)

% *Oxygen and carbon isotopic composition of dolomite and fluid*
%d13Cdol: Carbon isotopic composition of dolomite (permill)
%d13Cw: Total dissolved carbon isotopic composition of the fluid (permill)
%d18Odol: Oxygen isotopic composition of dolomite (permill)
%d18Ow: Oxygen isotopic composition of fluid (permill)
%d18Owapp: Apparent oxygen isotopic composition of fluid (permill)


% *System parameter*
%P: Porosity of the system dolomite+fluid (from 0 to 1)
%Tsys: Temperature of the system dolomite+fluid (°C)

% *Fractionation factor*
%a: Oxygen isotope fractionation factor between dolomite and fluid
%aC: Carbon isotope fractionation factor between dolomite and fluid

% *Fluid-rock interaction*
%n: Model iteration (i.e pore volume of fluid) ([ ])
%N: Cumulative water to rock weight ratio at a given stage in the
%interaction process (dimensionless)
%nf: Total number of iterations ([ ])
%V0: Initial dolomite volume (from 0 to 1)
%Vrec: Dolomite volume which react with the fluid (from 0 to V0)
%Vnr: Dolomite volume which does not react with the fluid (from 0 to V0)
%inc: Iteration at which dolomite has completely recrystallized ([ ])
%Ftot: total weight fraction of fluid in the system relative to V0 ([ ])
%F: weight fraction of fluid in the system relative to Vrec ([ ])

% *Equilibrium*
%IEQ_O: extent to which dolomite oxygen isotopes have reach equilibrium
%with the oxygen isotopes of the fluid.
%IEQ_C: extent to which dolomite carbon isotopes have reach equilibrium
%with the carbon isotopes of the fluid.
%IEQ_T47: extent to which dolomite clumped isotopes have reach equilibrium
%with the fluid temperature.

%% SECTION 1: INITIAL CONDITIONS OF THE SYSTEM DOLOMITE (D0) AND FLUID (W0)

%*Initial carbon concentration of the system dolomite+fluid*
cCf=cst.cCf; 
cCdol=cst.cCdol; 

%*Initial oxygen concentration of the system dolomite+fluid*
cOf=cst.cOf; 
cOdol=cst.cOdol; 
Rhodol=cst.Rhodol;
Rhof=cst.Rhof; 

%*Initial isotopic composition of dolomite and fluid*
[~,d18Odol0,d18Ow0]=matthewDW(ic.T,'',ic.d18Ow);
D470=daviesT2D(ic.T); 
d13Ctdc=ic.d13Ctdc; 
d13Cdol0=ic.d13Cdol;

%*Initial temperature, porosity and dolomite volume of the system* 
T0=ic.T; 
V0=1-ic.P; 
P=ic.P;

%*Initialization of the case number*
j=0;

%*The user decides the dolomite volume which recrystallize 
%in each pore volume of fluid n*
x=input(['Enter the dolomite volume which recrystallizes at each n'...
    '\n Enter 100 if recrystallization is 100% complete at each n'...
    '\n Enter a number <100 if dolomitization is incomplete'...
    '\n (e.g. enter 1 if 1% of dolomite recrystallizes)\n']);
%%
if x==100
    inc=0;
else
    inc=V0/(x/100);
end

nf=10000; %total number of runs

%% SECTION 2: RUN THE MODEL FOR DIFFERENT CASES

for Tsys=ic.T+1:1:100
    
    j=j+1; %case number
    
    %Initialization of the lists   
    listTsys=[];listd18O=[];
    listN=[]; listF=[];listd18Odol=[];listd18Owapp=[];listd13Cdol=[];
    listT47=[];listVrec=[];listIEQ_O=[]; listIEQ_C=[];listIEQ_T47=[];
   
    %% INITIALIZATION OF THE MODELS PARAMETERS
    n=0;
    Vrec=0;
    d18Ow=d18Ow0;
    d18Odol=d18Odol0;
    d13Cdol=d13Cdol0;
    
    % Compute a at Tsys
    [a,~,~]=matthewDW(Tsys,'',''); 
    aC=1; 
    
    % Compute the dolomite oxygen isotopic composition at equilibrium 
    [~,d18Oeq,~]=matthewDW(Tsys,'',d18Ow0);
    
    %% ITERATIVE CALCULATION OF ISOTOPIC ABUNDANCES 
    
    while n<nf            
            n=n+1; %update iteration number 
            
            if n<inc
                Vrec=(Vrec+(V0/inc));
                Vrec=round(Vrec,10); 
                Vnr=V0-Vrec; 
            
                else % after n=inc, recrystallization is complete
                Vrec=V0; Vnr=0;  
            end

            % Compute the weight fraction of fluid
            F=P*1/(P+(Vrec*Rhodol)); %F=Ftot if complete recrystallization
            Ftot=P*1/(P+((1-P)*Rhodol));
            
            % Compute the cumulative water to rock ratio N
            N=n*Ftot/(1-Ftot); %cumulative water to rock ratio
            
            %Compute the oxygen and carbon concentration in the system
            %dolomite+fluid
            cO=(F*cOf)+(1-F)*cOdol;
            cC=(F*cCf)+(1-F)*cCdol;
            
            %Mass balance oxygen
            d18O=((d18Ow*F*cOf)+(d18Odol*(1-F)*cOdol))/cO;
            d18Odolrec=((d18O*cO*a)-(1000*cOf*F*(1-a)))...
                /((cOdol*(1-F)*a)+cOf*F);
            d18Odol=((Vrec*d18Odolrec)+(Vnr*d18Odol0))/V0;
            
            %Mass balance carbon
            d13C=((d13Ctdc*F*cCf)+(d13Cdol*(1-F)*cCdol))/cC;
            d13Cdolrec=((d13C*cC*aC)-(1000*cCf*F*(1-aC)))...
                /((cCdol*(1-F)*aC)+cCf*F);
            d13Cdol=((Vrec*d13Cdolrec)+(Vnr*d13Cdol0))/V0;
            
            %Dolomite "clumped" temperature
            T47=((Vrec*Tsys)+(Vnr*T0))/V0;

            % Apparent oxygen isotopic composition of the fluid
            [~,~,d18Owapp] = matthewDW(T47,d18Odol,'') ;

            %RESULT
            
            listTsys(n,1)=Tsys;
            listN(n,1)=N;
            listF(n,1)=F;
            listd18Odol(n,1)=d18Odol-30.92;
            listd18Owapp(n,1)=d18Owapp;
            listd13Cdol(n,1)=d13Cdol;
            listT47(n,1)=T47;            
            listVrec(n,1)=Vrec;
            listIEQ_O(n,1)=abs(d18Odol-d18Odol0)/abs(d18Oeq-d18Odol0);
            listIEQ_C(n,1)=abs(d13Cdol-d13Cdol0)/abs(d13Ctdc-d13Cdol0);
            listIEQ_T47(n,1)=abs(T47-T0)/abs(Tsys-T0);
    end
    %% Final Array
    M(:,:,j)=[listTsys,listN,listF,listd18Odol,listd18Owapp,listd13Cdol,...
        listT47,listVrec,listIEQ_O,listIEQ_C,listIEQ_T47];

end

%% SECTION 3: FIND SELECTED N VALUES IN M

%In practice, the user has a dataset of dolomite oxygen and clumped
%isotopes and wants to know if these dolomites have recrystallized, 
% and if so, under which water to rock ratio (N). 
%A good way to find N is to plot the dataset on a T=f(d18Odol) 
%and to use the model to plot N isolines on the same figure.

%It is straigthforward to plot T47=f(N) and d18Odol=f(N) for the 
%different cases using the array M. But it is not as straightforward to
%plot N isolines on a T47=f(d18O) plot. To be able to plot N isolines on a 
%T=f(d18Odol) plot, we built the following two-steps process:

%The first step is to extract few selected values of N from all the
%results stored in M (e.g. N=1,5,100 etc.). 
%The second step is to create a best fit function between T47 and d18O for
%each value of N.

% The user select the values of N in chooseN: 
chooseN={100,50,10,5,1,F/(1-F)}; 

%The following loop is going to read through all the M matrices to find the
%the row of M corresponding to selected values of N and store them in a
%matrix A

A=NaN(nf,size(M,2),size(chooseN,2)); %initialization of A

for j=1:size(chooseN,2) %read through all the selected N values
    
    ChosenRow=[];
    
    for i=1:1:size(M,3) %read through all the cases stored in M 
        
        %Find the row index of M where the selected value of N is located.
        %Careful! The model outputs stored in M may not contain 
        %the selected value of N. It is necessary to look for the closest
        %value between the selected N value and the model N value:
        ab=abs(chooseN{j}-M(:,2,i));
        idx=find(ab==min(ab)); 

        %If the difference between the selected value of N 
        %and the value of N stored in M is bigger than 0.1, break the loop.
        if ab(idx)>0.1 
            break
        end
        
        %If the value of N stored in M is equal or close (+/-0.1) to the
        %selected N value, store the row of M in ChosenRow.
        
        ChosenRow=[ChosenRow;M(idx,:,i)];
    end
    
    %Store all the chosen rows in A
    A(1:size(ChosenRow,1),:,j)=ChosenRow; 
    
end

%% SECTION 4: CREATE BEST FIT FUNCTIONS FOR N ISOLINES

Fits=[];
%Fitfunc: function to fit with constant parameters x(1), x(2) and x(3)
Fitfunc=@(x,xdata)x(1)*power(xdata,2)+x(2)*xdata+x(3); 
    
for j=1:size(chooseN,2)
    
    xdata=[d18Odol0-30.92;A(~isnan(A(:,1,1)),4,j)]; %d18Odol
    ydata=[ic.T;A(~isnan(A(:,1,1)),7,j)];%N
    
    xo=[-1,1,20]; %first guess for x(1), x(2) and x(3)
    options=optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000);
    x=lsqcurvefit(Fitfunc,xo,xdata,ydata,[],[],options);  
    
    GOF=goodnessOfFit(Fitfunc(x,xdata),ydata,'NMSE');
    
    Fits=[Fits;chooseN{j},x,{xdata},{ydata},GOF];  
end

%% SECTION 5: PLOTS
% PLOT T($\Delta_{47}$) AS A FUNCTION OF $\delta^{18}O_{dol}$
col=vectorcolor.col;
 
figure(1)

plot_isolines_d18Ow([-Inf Inf],[-Inf Inf]);
hold on;

lgd=[];
for i=size(M,3)-40:-5:1
    p1=plot(M(:,4,i),M(:,7,i),'k-.','LineWidth',1.5,'Color',col(i,:),...
        'DisplayName',[sprintf('%g', M(1,1,i)) 'ºC']); hold on;
    lgd=[lgd p1];
end

for i=1:size(Fits,1)
    scatter(cell2mat(Fits(i,3)),cell2mat(Fits(i,4)),'ko'); hold on;
    p2=plot(cell2mat(Fits(i,3)),Fitfunc(cell2mat(Fits(i,2)),...
        cell2mat(Fits(i,3))),'','LineWidth',3,'DisplayName',...
        ['N=', num2str(Fits{i},'%.2G')]); hold on;
    lgd=[lgd p2];
end

plot(d18Odol0-30.92,ic.T,'MarkerSize',20,'Marker','d','MarkerFaceColor',...
    [1 0.4 0.6],'Color','k','LineWidth',2);
text(d18Odol0-30.92,ic.T-3,'D0','FontSize',20,'FontWeight','bold',...
    'HorizontalAlignment','center')

ylabel('T(\Delta_{47 dol}) (ºC)','Fontsize',20);
xlabel('\delta^{18}O_{dol}','Fontsize',20);
xlim([-8 3]); 
ylim([15 60]);
set(findobj(gcf,'type','axes'),'FontSize',20,'FontWeight','Normal', ...
    'LineWidth', 2,'XDir','normal','TickLength',[0.01, 0.01]);
l=legend(lgd);
l.Location='West'; l.FontSize=10;

% PLOT EQUILIBRIUM
figure(3)
lgd=[];

for j=1:1:size(M,3)
        p1=plot([0;M(:,2,j)],[0;M(:,9,j)],'LineWidth',2,'Color','r',...
            'DisplayName','\delta^{18}O_{dol}'); hold on;
        p2=plot([0;M(:,2,j)],[0;M(:,10,j)],'LineWidth',2,'Color','g',...
            'DisplayName','\delta^{13}C_{dol}'); hold on;
        p3=plot([0;M(:,2,j)],[0;M(:,11,j)],'LineWidth',2,'Color','b',...
            'DisplayName','T(\Delta_{47dol})'); hold on;
end
l=legend([p1(1) p2(1) p3(1)]);
l.Location='East';
xlim([0 500])
xlabel('N');
ylabel('IEQ');
set(findobj(gcf,'type','axes'),'FontSize',20,'FontWeight','Normal',...
    'LineWidth', 2,'XDir','normal','TickLength',[0.01, 0.01],...
'XScale', 'log');
grid('on')

% EVOLUTION OF MODEL PARAMETERS DURING FLUID/ROCK INTERACTION
figure(4)

subplot(2,2,1)
for k=size(M,3)-40:-5:1
    plot(M(:,2,k),M(:,1,k),'k-.','LineWidth',2,'Color',col(k,:)); hold on;
end
ylim([0 100])
xlabel('N'); ylabel('T_{sys} (ºC)'); title('T_{sys} = f(N)')

subplot(2,2,2)
plot(M(:,2,1),d18Ow0*ones(1,size(M(:,2,1),1))','r-','LineWidth',2);
ylabel('\delta^{18}O_{w (sys)}'); xlabel('N'); 
title('\delta^{18}O_{w (sys)}=f(N)')

subplot(2,2,3)
plot(M(:,2,1),M(:,8,1)./(1-P),'r-','LineWidth',2); hold on;
ylim([0 1])
ylabel('V_{rec}'); xlabel('N');title('V_{rec} = f(N)')

subplot(2,2,4)
plot(M(:,2,1),M(:,3,1),'r-','LineWidth',2); hold on;
ylim([0 1])
ylabel('F'); xlabel('N'); title('F = f(N)')

set(findobj(gcf,'type','axes'),'FontSize',15,'FontWeight','Normal',...
    'LineWidth', 1.5,'XDir','normal','TickLength',[0.01, 0.01]);
