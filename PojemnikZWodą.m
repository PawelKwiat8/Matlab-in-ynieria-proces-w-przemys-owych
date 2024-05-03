clc;clear;
fp=0.5; % czestotliwość probkowania
Tp=1/fp;
t = 0:(Tp):(500-(Tp));
N = length(t);
Czas=0;

wyp=0.8; %wypełnienie kubka kawą
mom=2000; %moment dolania mleka [s]
Cp_w=4190; %ciepło właściwe wody [J/(kgK)]
Cp_s=900; %ciepło właściwe ścianek [J/(kgK)]
r=0.06; % promień kubka [m]
h=0.1; %wysokość kubka [m]
d=0.003; % szerokość scianek [m]
C_s = 5.667*0.92; %stała promieniowania porcelana (C0*stopień czarności) [W/(m^2*K^4)]
C_w = 5.667*0.98; %stała promieniowania woda  [W/(m^2*K^4)]
C_p = 5.667*0.93;%stała promieniowania tynki  [W/(m^2*K^4)]
rho_w=997; % gestość wody[kg/m^3]
rho_s=2300; % gestość scianek po[kg/m^3]
alfa_woda_ks=400; % Współczynnik wnikania ciepła woda,Konwekcja swobodna  [W/(mK)] (250 ... 600)
alfa_woda_kw=1000; % Współczynnik wnikania ciepła woda,Konwekcja wymuszona  [W/(mK)] (500 ... 10 000)
alfa_powietrze_ks=15; % Współczynnik wnikania ciepła woda,Konwekcja swobodna  [W/(mK)] (3 ... 20)
alfa_powietrze_kw=120; % Współczynnik wnikania ciepła woda,Konwekcja wymuszona  [W/(mK)] (10 ... 150)
Tw=zeros(1,N); %temperatura wody[K]
Tw(1)=95+273.15;
Ts=zeros(1,N); %temperatura ścianek[K]
T_p=ones(1,N)*(4+273.15); %temperatura powietrza[K] , constant
Ts(1)=23+273.15;
Qw = zeros(1,N); %ciepło wody  [J]
Qw_k = zeros(1,N); %ciepło wody konwekcyjne [J]
Qw_r = zeros(1,N); %ciepło wody radiacyjne oddane[J]
Qp = zeros(1,N); % ciepło powietrza [J]
Qp_k = zeros(1,N); % ciepło powietrza konwekcyjne [J]
Qp_r = zeros(1,N); % ciepło powietrza radiacyjne oddane[J]
V_w=pi*r^2*h*wyp; %objętość wody[m^3]
V_s=pi*(r+d)^2*h*wyp - V_w ; %objętość ścianek[m^3]
Ab=2*pi*(r+d)*h; %pole powierzchni bocznej [m^2]
Ab_w=2*pi*(r)*h*wyp; %pole powierzchni bocznej wewnatrz woda [m^2]
Ag=pi*(r)^2; %pole powierzchni górnej [m^2]
mom=mom*fp; %kalkulacja czasu na moment próbkowania

for i = 1:N
    Qw_r(i)= C_w*Ab_w*(Tw(i)/100)^4 - C_s*Ab_w*(Ts(i)/100)^4 + C_w*Ag*(Tw(i)/100)^4 - C_p*Ag*(T_p(i)/100)^4;  % ciepło radiacyjne woda-scianki i woda-otoczenie(góra)
    Qw_k(i) = alfa_woda_ks * Ab_w * (Tw(i) - Ts(i)) + alfa_woda_ks * Ag * (Tw(i) - T_p(i)) ;  % ciepło konwekcyjne woda - scianki i woda otoczenie(góra)
    Qw(i) = Qw_r(i)+Qw_k(i); % Całkowite Ciepło oddane przez wode

    dTw_dt = Qw(i) / (V_w * rho_w * Cp_w);

    if(i~=N)
        Tw(i + 1) = Tw(i) - dTw_dt * Tp; % Nowa temperatura wody
    end

    Qp_r(i)=C_s*Ab*(Ts(i)/100)^4-C_p*Ab*(T_p(i)/100)^4 ; % ciepło radiacyjne otoczenie - scianki
    Qp_k(i)= -alfa_powietrze_kw * Ab * (T_p(i)-Ts(i)); % ciepło konwekcyjne otoczenie - scianki
    Qp(i) = Qp_k(i) + Qp_r(i); % całkowite Ciepło pobrane do powietrza

    dTs_dt = (Qw(i) - Qp(i)) / (V_s * rho_s * Cp_s);
    if(i~=N)
        Ts(i + 1) = Ts(i) + dTs_dt * Tp; % Nowa temperatura ścianek
    end

    if Tw(i)>60+273.15
        Czas=Czas+1;
    end

    if i==mom %moment wlania mleka do kawy
        % Ab_w=Ab;
        % V_w=pi*r^2*h;
        % Qw(i)=Qw(i-1);
        % Tw(i + 1)=(1-wyp)*(4+273.15)+wyp*Tw(i);
        % Qp(i)=Qp(i-1);
        % Ts(i + 1)=Ts(i);
        %T_p(mom:N) = (23+273.15);
    end
end
ond = zeros(1,N) + 60;
Czas1=Czas/5;
tmin=t./60;
Tw_C=Tw-273.15;
Ts_C=Ts-273.15;
T_p_C=T_p-273.15;


w = sum(Qp_r)/sum(Qp); %stosunek Qp_r do Qp
z = sum(Qw_r)/sum(Qw); %stosunek Qw_r do Qw

figure;
hold on;

plot(tmin, Tw_C, 'b', 'LineWidth', 1.5);
plot(tmin, Ts_C, 'r', 'LineWidth', 1.5);
plot(tmin, T_p_C, 'g', 'LineWidth', 1.5);
plot(tmin, ond, '--k')
plot(tmin, xline(Czas/60/fp), '--k')
plot(Czas/60/fp,60,'k*')

xlabel('Czas [min]');
ylabel('Temperatura wody[C]');
legend('Tw','Ts','Tp');
grid on;
hold off;
title(['temperatura wody od czasu (czas po którym osiągnięto temperaturę przyjemności: ', num2str(Czas1), 's)'])

figure;
hold on;
plot(tmin, Qp, 'r', 'LineWidth', 1.5);
plot(tmin, Qp_r, '--', 'Color', [1 0.5 0], 'LineWidth', 0.5); % Pomarańczowy
plot(tmin, Qp_k, '--', 'Color', [0.5 0 1], 'LineWidth', 0.5); % Fioletowy
plot(tmin, Qw, 'b', 'LineWidth', 1.5);
plot(tmin, Qw_r, 'c--', 'LineWidth', 0.5);
plot(tmin, Qw_k, 'g--', 'LineWidth', 0.5);


xlabel('Czas [min]');
ylabel('ciepło [J]');
legend('Qp','Qp_r','Qp_k','Qw','Qw_r','Qw_k');
grid on;
hold off;
title(['ciepło od czasu , Qp_r / Qp = ',num2str(round(w,3)),'  Qw_r / Qw = ',num2str(round(z,3)) ])




