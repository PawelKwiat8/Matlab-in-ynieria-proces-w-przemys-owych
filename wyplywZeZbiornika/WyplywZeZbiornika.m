clc;clear;

Data=readtable('IPP/wyplyw_ze_zbiornika/dysza.txt', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');

Data.Properties.VariableNames = {'over_p', 'q'};

% Konwersja danych z kolumn 'over_p' i 'q' z tekstowych na numeryczne
Data.over_p = str2double(strrep(Data.over_p, ',', '.'));
Data.q = str2double(strrep(Data.q, ',', '.'));

D=Data(25:end-17,:);

q_DN=(D.q) ; % przepływ [Nl/min]
p0_D=100000+D.over_p*100000; %cisnienie w zbiorniku [Pa]

%wconst
ad=1.4;
rho1 = 1.293; % gestość powitrza [kg/m^3]
p1=100000;%[Pa]


%zamiana Nl na l
rho0_D=rho1.*(p0_D./p1).^(1./ad); %gestość w zbiorniku[kg/m^3]
N=sqrt(rho0_D./rho1).*sqrt(p0_D./p1);% wsp do zamiany 
q_D=q_DN./N;% przepływ [l/min]


% Liniowe dopasowanie danych
linFit = fit(p0_D, q_D, 'poly1');% 'poly1' oznacza wielomian pierwszego stopnia, czyli model liniowy

p0_linfit=min(p0_D):500:max(p0_D);
q_linfit = feval(linFit, p0_linfit);

figure;
plot(p0_D, q_D);
hold on;
plot(p0_linfit, q_linfit,MarkerSize=0.5);
title(' przepływ od ciśnienia');
xlabel('Ciśnienie w zbiorniku p0 [Pa]');
ylabel('Przepływ q [l/min]');
%% Dane obliczeniowe
A=0.8*10^(-6); %przekrów [m^2]
rho1 = 1.293; % gestość powitrza [kg/m^3]
p1=100000;%[Pa]
p0= min(p0_D):500:max(p0_D);%[Pa]
rho0=rho1.*(p0./p1).^(1./ad);%gestość w zbiorniku[kg/m^3]
beta = (2/(ad+1))^(ad/(ad-1));


initialC = 1; % Początkowa wartość c, dostosuj zgodnie z potrzebami

optimalC = fminsearch(@(c) fitError(c, p0, q_linfit', A, ad, rho0), initialC);

disp(['Optymalna wartość c: ', num2str(optimalC)]);

q=(A*optimalC*sqrt(ad.*p0.*rho0*(2/(ad+1))^((ad+1)/(ad-1))))./rho0*1000*60;

figure
plot(p0,q);
hold on;
plot(p0_linfit, q_linfit,MarkerSize=0.5);
title(' przepływ od ciśnienia współczynnik kontrakcji: ',optimalC);
xlabel('Ciśnienie w zbiorniku p0 [Pa]');
ylabel('Przepływ q [l/min]');
legend('obliczenia','eksperyment')

%%
function error = fitError(c, p0, q_linfit, A, ad, rho0)
    q_calc = (A * c * sqrt(ad * p0 .* rho0 * (2/(ad+1))^((ad+1)/(ad-1)))) ./ rho0 * 1000 * 60; % Obliczenie przepływu [l/min]
   
    error = sum((q_calc - q_linfit).^2); % Suma kwadratów różnic
end