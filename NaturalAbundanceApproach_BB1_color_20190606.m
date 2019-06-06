% Add dismutation into (Peters et al. 2016) model
% Data from ETSP2013 BB1
% Xin Sun (last update: June 6th, 2019)


% 1. INPUTS for the model

delta15Nitrate = delta15Nitrate0./1000; % no unit
delta15Nitrite = delta15Nitrite0./1000; % no unit

% ! Difference-1 add Dismutation, and make assumption for fractionation factors:
alpha_DIS = 0.978; %set to use alpha_nxr, since OMZ NOB have Nxr clustered with 'normal' NOB
alpha_DISN2 = 1.018; %set to use alpha_nir since OMZ NOB have NirK
% ! Difference-2: remove 'canonical' NO2- oxidation in ODZ core from model since there is no O2 available
% ! Difference-3 Anammox stoichiometry change based on Oshiki et al. 2016:
c = 0.16 % it was 0.3 in Peters et al. (2016)'s model


% These numbers were from caption of Fig.8 or Table 1 in Peters et al. 2016
% alpha_NXR = 0.978;
alpha_NIR = 1.018; 
alpha_AMX = 1.016;
alpha_NXRAMX = 0.969;
alpha_NAR = 1.018;

% Nitrate1514 means [15NO3]/[14NO3]
Nitrate1514 = (1 + delta15Nitrate) .* 0.003667./(1-0.003667);
Nitrite1514 = (1 + delta15Nitrite) .* 0.003667./(1-0.003667);
Nitrate15 = Nitrate ./ (1./Nitrate1514+1);
Nitrite15 = Nitrite ./ (1./Nitrite1514+1);
Nitrate14 = Nitrate - Nitrate15;
Nitrite14 = Nitrite - Nitrite15;

% 2. Net production or consumption rates of nitrite and nitrate predicted by vertical advection + diffusion
% use standard case values from Table 1 in Peters et al. 2016
w = -2*10^(-7); % unit: m/s
D = 4*10^(-5); % unit: m2/s
% central differencing for advection
for i = 2:(length(Z)-1);
    R14Nitrate(i) = 24*60*60* (w*(Nitrate14(i+1)-Nitrate14(i-1))/(Z(i+1)-Z(i-1)) - ((D*(Nitrate14(i+1)-Nitrate14(i))/(Z(i+1)-Z(i)))-(D*(Nitrate14(i)-Nitrate14(i-1))/(Z(i)-Z(i-1))))/((Z(i+1)-Z(i-1))/2));
    R14Nitrite(i) = 24*60*60* (w*(Nitrite14(i+1)-Nitrite14(i-1))/(Z(i+1)-Z(i-1)) - ((D*(Nitrite14(i+1)-Nitrite14(i))/(Z(i+1)-Z(i)))-(D*(Nitrite14(i)-Nitrite14(i-1))/(Z(i)-Z(i-1))))/((Z(i+1)-Z(i-1))/2));
    R15Nitrate(i) = 24*60*60* (w*(Nitrate15(i+1)-Nitrate15(i-1))/(Z(i+1)-Z(i-1)) - ((D*(Nitrate15(i+1)-Nitrate15(i))/(Z(i+1)-Z(i)))-(D*(Nitrate15(i)-Nitrate15(i-1))/(Z(i)-Z(i-1))))/((Z(i+1)-Z(i-1))/2));
    R15Nitrite(i) = 24*60*60* (w*(Nitrite15(i+1)-Nitrite15(i-1))/(Z(i+1)-Z(i-1)) - ((D*(Nitrite15(i+1)-Nitrite15(i))/(Z(i+1)-Z(i)))-(D*(Nitrite15(i)-Nitrite15(i-1))/(Z(i)-Z(i-1))))/((Z(i+1)-Z(i-1))/2));
end

% 3. Run model with dismutatition, without canonical nitrite oxi:
% below 5 equations used to solve 4 unknown rates: Fnir14 Fnar14 Famx14 Fdis14
%{
%Equations:
% using equations 7a,b,c, 11a,b in Peter et al. (2016), and same to their model, assuming R14NH4 = 0)

0          =                     0.11 * 14Fnir +                    0.07 * 14Fnar                                                         -14Famx;
R14Nitrite =                          - 14Fnir +                           14Fnar                                                 -(1+c) * 14Famx                     - (5/3)*14Fdis;
R14Nitrate =                                                             - 14Fnar +                                                    c * 14Famx +                           14Fdis;
R15Nitrite = -1/alpha_NIR*Nitrite1514 * 14Fnir + 1/alpha_NAR*Nitrate1514 * 14Fnar - (1/alpha_AMX*Nitrite1514+c/alpha_NXRAMX*Nitrite1514) * 14Famx - (1/alpha_DIS*Nitrite1514 + 1/alpha_DISN2*Nitrite1514*(2/3))*14Fdis;
R15Nitrate =                                   - 1/alpha_NAR*Nitrate1514 * 14Fnar +                           c/alpha_NXRAMX*Nitrite1514 * 14Famx + 1/alpha_DIS*Nitrite1514 * 14Fdis;
%}
% solve above equations:
for i=1:(length(R14Nitrite));
    d = [0
    R14Nitrite(i)
    R14Nitrate(i)
    R15Nitrite(i)
    R15Nitrate(i)];
C = [0.11 0.07 -1 0
    -1 1 -(1+c) -(5/3)
    0 -1 c 1
    (-1/alpha_NIR*Nitrite1514(i)) (1/alpha_NAR*Nitrate1514(i)) (-(1/alpha_AMX*Nitrite1514(i)+c/alpha_NXRAMX*Nitrite1514(i))) (-(1/alpha_DIS*Nitrite1514(i) + 1/alpha_DISN2*Nitrite1514(i)*(2/3)))
    0 (-1/alpha_NAR*Nitrate1514(i)) (c/alpha_NXRAMX*Nitrite1514(i)) (1/alpha_DIS*Nitrite1514(i))];
% solve x = [14Fnir 14Fnar 14Famx 14Fdis] 
% using nonnegative least squares (NNLS) optimization routine in MATLAB.
[x,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d);

% Print results:
Fnir14(i) = x(1);
Fnir(i) = Fnir14(i) * (1 + 1/alpha_NIR*Nitrite1514(i));
NitriteReduction_nMNperDay(i) = Fnir(i) * 1000;

Fnar14(i) = x(2);
Fnar(i) = Fnar14(i) * (1 + 1/alpha_NAR*Nitrate1514(i));
NitrateReduction_nMNperDay(i) = Fnar(i) * 1000;

Famx14(i) = x(3);
Famx(i) = Famx14(i) * (1 + 1/alpha_AMX*Nitrite1514(i));
Anammox_nMNperDay(i) = Famx(i) * 1000;% anammox rates here is nM-NO2- used per day.

Fdis14(i) = x(4);
Fdis(i) = Fdis14(i) * (1 + 1/alpha_DIS*Nitrite1514(i));
NO2_dismutation_nMNO3perDay(i) = Fdis(i) * 1000;
end

% plot data
figure1 = figure;
% Create subplot
subplot1 = subplot(1,4,2,'Parent',figure1,'FontSize', 16, 'Ydir', 'reverse','XTick',(0:5:15));
xlim(subplot1,[0 15]);
ylim(subplot1,[0 450]);
box(subplot1,'on');
hold(subplot1,'on');
% ODZ grey areas
%
ODZ=0:100;
curve1=ODZ-ODZ+130;
curve2=ODZ-ODZ+370;
plot(ODZ,curve1);
hold on
plot(ODZ,curve2);
ooo = [ODZ, fliplr(ODZ)];
inBetween = [curve1, curve2];
fill(ooo, inBetween, [0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
hold on
%
plot(NO2_dismutation_nMNO3perDay, Z(1:(length(NO2_dismutation_nMNO3perDay))), 'blueo-','MarkerSize',12,'Linewidth',2);
xlabel('NO_2^- dismuation (nM-NO_3^-/d)','FontSize',20,'FontWeight','normal');
% ylabel('Depth (m)','FontSize',22,'FontWeight','bold');


% plot nitrate reduction rate
subplot2 = subplot(1,4,1,'Parent',figure1,'FontSize', 16, 'Ydir', 'reverse','XTick',(0:20:60));
xlim(subplot2,[0 60]);
ylim(subplot2,[0 450]);
box(subplot2,'on');
hold(subplot2,'on');
% ODZ grey areas
%
ODZ=0:120;
curve1=ODZ-ODZ+130;
curve2=ODZ-ODZ+370;
plot(ODZ,curve1);
hold on
plot(ODZ,curve2);
ooo = [ODZ, fliplr(ODZ)];
inBetween = [curve1, curve2];
fill(ooo, inBetween, [0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
hold on
%
plot(NitrateReduction_nMNperDay, Z(1:length(NitrateReduction_nMNperDay)), 'blacko-','MarkerFaceColor','black','MarkerSize',12,'Linewidth',2);
xlabel('NO_3^- reduction (nM-N/d)','FontSize',20,'FontWeight','normal');
ylabel('Depth (m)','FontSize',22,'FontWeight','normal');

%plot rates of nitrogen loss processes, including nitrite dismutation, denitrification and anammox
subplot3 = subplot(1,4,3,'Parent',figure1,'FontSize', 16, 'Ydir', 'reverse','XTick',(0:5:20));
xlim(subplot3,[0 10]);
ylim(subplot3,[0 450]);
box(subplot3,'on');
hold(subplot3,'on');
% ODZ grey areas
%
ODZ=0:100;
curve1=ODZ-ODZ+130;
curve2=ODZ-ODZ+370;
plot(ODZ,curve1);
hold on
plot(ODZ,curve2);
ooo = [ODZ, fliplr(ODZ)];
inBetween = [curve1, curve2];
fill(ooo, inBetween, [0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
hold on
%
plot(NitriteReduction_nMNperDay, Z(1:length(NitriteReduction_nMNperDay)), 'blacko-','MarkerFaceColor','black','MarkerSize',12,'Linewidth',2); %'Color',[253/255 174/255 97/255], 
hold on
% Note that previous anammox is nM-nitrite consumed, here times 2 to make it into nM-N of N2 produced per day.
plot(2*Anammox_nMNperDay, Z(1:length(Anammox_nMNperDay)), 'Redo-','MarkerSize',12,'Linewidth',2);
xlabel('N_2 production (nM-N/d)','FontSize',20,'FontWeight','normal');
%ylabel('Depth (m)','FontSize',22,'FontWeight','bold');
hold on
plot (NO2_dismutation_nMNO3perDay./3.*2, Z(1:length(NO2_dismutation_nMNO3perDay)), 'blueo-','MarkerSize',12,'Linewidth',2)

% plot total nitrogen loss rates in models with or without dismutation.
subplot3 = subplot(1,4,4,'Parent',figure1,'FontSize', 16, 'Ydir', 'reverse','XTick',(0:5:20));
xlim(subplot3,[0 20]);
ylim(subplot3,[0 450]);
box(subplot3,'on');
hold(subplot3,'on');
% ODZ grey areas
%
ODZ=0:100;
curve1=ODZ-ODZ+130;
curve2=ODZ-ODZ+370;
plot(ODZ,curve1);
hold on
plot(ODZ,curve2);
ooo = [ODZ, fliplr(ODZ)];
inBetween = [curve1, curve2];
fill(ooo, inBetween, [0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
hold on
%
plot (NitriteReduction_nMNperDay+2*Anammox_nMNperDay+NO2_dismutation_nMNO3perDay./3.*2, Z(1:length(NO2_dismutation_nMNO3perDay)), 'blueo-','Linewidth',2,'MarkerSize',12)
hold on
plot (Nlosswithoutdismuation, Z(1:length(NO2_dismutation_nMNO3perDay)), 'black*--','Linewidth',2,'MarkerSize',12)
xlabel('Total N_2 production (nM-N/d)','FontSize',20,'FontWeight','normal');

