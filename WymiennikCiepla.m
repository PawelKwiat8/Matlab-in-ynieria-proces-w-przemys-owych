clc;clear;

fp=0.5; % czestotliwość probkowania
Tp=1/fp;
l = 0:(Tp):(100-(Tp));
N = length(l);
k_R = 1.8; %nibywzomcnienie nibyregulatora
K = 5; %współczynnik wynikający z rodzaju wymiennika
r = 0.02; %promień rurki wymiennuka ciepła [m]
G1 = 2; %zużycie cieczy1 grzejącej [kg]
G2 = 2; %zużycie cieczy2 ogrzewanej[kg]
Cw1=4190; %ciepło właściwe cieczy1 grzejącej[J/(kgK)]
Cw2=4190; %ciepło właściwe cieczy2 ogrzewanej[J/(kgK)]
T_in = 273.15+80; %Temperatura cieczy2 grzejącej na wejściu [K]
t_in = 273.15+20; %Temperatura cieczy ogrzewanej na wejściu [K]

T_out = zeros(1,N);%Temperatura cieczy2 grzejącej na wyjściu [K]
T_out(1) = T_in;
t_out = zeros(1,N);%Temperatura cieczy ogrzewanej na wyjściu [K]
t_out(1) = t_in;


H_gwiazdka = zeros(1,N);
A =zeros (1,N); %powierzchnia wymiany ciepła [m^2]
e = zeros(1,N); %nibyodchyłka regulacji
H = zeros(1,N); %wnikanie ciepła [J]
Theta_min = zeros(1,N); %Różnica temperatur na wyjściu [K]
Theta_max = T_in-t_in; %Różnica temperatur na wejściu [K]
Theta_sr = zeros(1,N); %Średnia różnica temperatur [K]
Theta_sr(1) = Theta_max;

for i = 1:N
    if i>1
        e(i)=H_gwiazdka(i-1)-H(i);
        H(i)=H(i-1)+k_R*e(i);
    end
    A(i) = 2*pi*r*l(i); %powierzchnia wymiany ciepła [m^2]
    T_out(i)=T_in-(H(i)/(G1*Cw1));
    t_out(i)=t_in+(H(i)/(G2*Cw2));
    Theta_min(i) = T_out(i)-t_out(i);
    if Theta_min(i) == Theta_max
        Theta_sr(i) = Theta_max;
    else
        Theta_sr(i) = (Theta_max-Theta_min(i))/log(Theta_max/Theta_min(i));
    end
    H_gwiazdka(i) = K*A(i)*Theta_sr(i);
    
end
% Plotting
figure;

% Heat input over time

plot(l, H, 'LineWidth', 2);
title('wymiana ciepła w zależności od długości wymiennika');
xlabel('l - długość [m]');
ylabel('H - ciepło [J]');

% Outlet temperatures over time
figure;
plot(l, T_out, 'LineWidth', 2); 
hold on;
plot(l, t_out, 'LineWidth', 2); 
hold off;
lgd = legend('T_out - Temperatura cieczy2 grzejącej na wyjściu', 't_out - Temperatura cieczy ogrzewanej na wyjściu '); % Tworzenie legendy
title('Temperatura na wyjściu w zależności od długości wymiennika'); % Tytuł wykresu
xlabel('długość [m]'); % Etykieta osi X
ylabel('Temperatura [K]'); % Etykieta osi Y

% Powiększenie legendy
set(lgd, 'FontSize', 12)

% Regulation error over time
figure;
plot(l, e, 'LineWidth', 2);
title('odchyłka regulacji');
xlabel('długość [m]');
ylabel('e [J]');


figure;
plot(l, Theta_sr, 'LineWidth', 2);
title('Średnia różnica temperatura');
xlabel('długość [m]');
ylabel('Theta_sr [K]');
