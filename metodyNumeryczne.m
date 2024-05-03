clear;close all;clc;
%Parametry toluen
A_T = 4.07827;
B_T = 1343.943;
C_T = -53.773;
x = 0.5;
%Parametry benzen
A_B = 4.72583;
B_B = 1660.652;
C_B = -1.461;

P_a = 100000; %Ciśnienie atmosferyczne [Pa]
T = 323.15:0.1:383.15; %zakres temperatur [K]

Ps_T=100000.*antoine(A_T,B_T,C_T,T); %ciśnienie par toluen [pa]
Ps_B=100000.*antoine(A_B,B_B,C_B,T); %ciśnienie par benzen [pa]

P=x.*Ps_T+(1-x).*Ps_B;

figure
plot(T,Ps_T)
hold on
plot(T,Ps_B)
plot(T,P)
hold off
legend('Toluen','Benzen','suma','Location','eastoutside')
xlabel('Temperatura[K]');
ylabel('cisnienie par[Pa]')

y=x.*Ps_T./P;

f=@(T_x)x*100000*antoine(A_T,B_T,C_T,T_x)+(1-x)*100000.*antoine(A_B,B_B,C_B,T_x)-P_a; %Funkcja do znalezienia punktu przecięcia

p = arrayfun(f,T);

%metoda bisekcji
[szukane_T_BI,c_history_BI,elapsedTime_BI] = bisectionMethod(f,300,400,0.01,T,p);

%% metoda newtona rapsona
[szukane_T_RN,c_history_RN,elapsedTime_RN] = N_R_Method(f,300,400,0.01,T,p);

%% flasi
[szukane_T_F,c_history_F,elapsedTime_F] = falsiMethod(f,300,400,0.01,T,p);

%% explicitConvergence
[szukane_T_eC,c_history_eC,elapsedTime_eC] = explicitConvergence(f,300,400,0.01,T,p,1.2);

%% zawartość procentowa par
y_toluen=x*100000*antoine(A_T,B_T,C_T,szukane_T_F)/P_a*100;
y_benzen=x*100000*antoine(A_B,B_B,C_B,szukane_T_F)/P_a*100;
disp(['Temperatura wrzenia mieszanki: ' num2str(szukane_T_eC) 'K'])
disp(['zawartość toluenu w parach: ' num2str(y_toluen) '%']);
disp(['zawartość benzenu w parach: ' num2str(y_benzen) '%']);

%%
function P =antoine(A,B,C,T) %ciśnienie z równania Antione'a
    P=10.^(A-(B./(T+C)));
end


%% implementacja bisekcji
function [c, c_history,elapsedTime] = bisectionMethod(f,a,b,error,T,p)

%standardowe
a=min(T);
b=max(T);
c=(a+b)/2;
tic;

while abs(f(c))>error
    if f(c)<0&&f(a)<0
        a=c;
    else
        b=c;
    end
    c=(a+b)/2;
end
elapsedTime=toc;


%ulepszona wersja 
a=min(T);
b=max(T);
c=(a+b)/2;
c_history = c;

while abs(f(c))>error
    if f(c)<0&&f(a)<0
        a=c;
    else
        b=c;
    end
    c=(a+b)/2;
    c_history = [c_history, c];
end

figure
plot(T, p+100000);
title(sprintf('Metoda bisekcji - Czas: %.5f s, liczba iteracji: %d', elapsedTime, length(c_history)));
xlabel('Temperatura[K]');
ylabel('cisnienie par[Pa]')
hold on;
plot(xlim, [100000 100000], '--k')
scatter(c_history, f(c_history)+100000 ,'+','r');
hold off;

end

%% metoda newtona rapsona
function [c, c_history,elapsedTime] = N_R_Method(f,a,b,error,T,p)

h = 1e-5; % Step size for numerical derivative
df = @(x) (f(x + h) - f(x - h)) / (2*h); % Central difference method

a=min(T);
b=max(T);
c=(a+b)/2;

tic;
while abs(f(c))>error
    c=c-f(c)/df(c);
end
elapsedTime=toc;

a=min(T);
b=max(T);
c=(a+b)/2;
c_history = c;
df_history= df(c);

while abs(f(c))>error
    c=c-f(c)/df(c);
    c_history = [c_history, c];
    df_history = [df_history,df(c)];
end

figure
plot(T, p+100000);
title(sprintf('Metoda Newtona Rapsona - Czas: %.5f s, liczba iteracji: %d', elapsedTime, length(c_history)));
xlabel('Temperatura[K]');
ylabel('cisnienie par[Pa]')
hold on;
plot(xlim, [100000 100000], '--k')
scatter(c_history, f(c_history)+100000 ,'+','r');

% Drawing tangent lines
for k = 1:length(df_history)
    
    x_start = c_history(k);
    y_start = f(c_history(k));

    x_intercept = x_start - y_start / df_history(k);
   
    if x_start < x_intercept
        x_range = linspace(x_start, x_intercept, 100);
    else
        x_range = linspace(x_intercept, x_start, 100);
    end
 
    y_line = df_history(k) * (x_range - x_start) + y_start;
    plot(x_range, y_line+100000, 'r--', 'LineWidth', 0.5); 
   
    text(c_history(k), y_start+100000, [' ' num2str(k)], 'VerticalAlignment', 'top', 'Color', 'r', 'FontSize', 7);
end

hold off;

end


%% metoda falsi
function [c, c_history,elapsedTime] = falsiMethod(f,a,b,error,T,p)

a=min(T);
b=max(T);
c = (a * f(b) - b * f(a)) / (f(b) - f(a));

tic;
while abs(f(c))>error
    c = (a * f(b) - b * f(a)) / (f(b) - f(a)); % Calculate the false position
    
    
    if f(a) * f(c) < 0
        b = c; 
    else
        a = c; 
    end
end
elapsedTime=toc;


a=min(T);
b=max(T);
c = (a * f(b) - b * f(a)) / (f(b) - f(a));
c_history = c;
figure
plot(T, p+100000);

xlabel('Temperatura[K]');
ylabel('cisnienie par[Pa]')
hold on;
plot([a, b], [f(a)+100000, f(b)+100000], 'Color', [0.5, 0.5, 0.5] , 'LineWidth', 0.25);
while abs(f(c))>error
    c = (a * f(b) - b * f(a)) / (f(b) - f(a)); % Calculate the false position
    c_history = [c_history, c]; % Store the current approximation
    
    % Decide which interval to keep for the next iteration
    if f(a) * f(c) < 0
        b = c; % Root lies between a and c
    else
        a = c; % Root lies between c and b
    end

    plot([a, b], [f(a)+100000, f(b)+100000],  'Color', [0.5, 0.5, 0.5] , 'LineWidth', 0.25);
end

title(sprintf('Metoda Falsi - Czas: %.5f s, liczba iteracji: %d', elapsedTime, length(c_history)));
plot(xlim, [100000 100000], '--k')
scatter(c_history, f(c_history)+100000 ,'+','r');
hold off;

end
%%
function [c, c_history, elapsedTime] = explicitConvergence(f, a, b, error, T, p, beta)
a=min(T);
b=max(T);
c=(a+b)/2;

tic;
while abs(f(c)) > error
    c_calc = (a * f(b) - b * f(a)) / (f(b) - f(a)); % Przykład z metody falsi

    c = c + beta * (c_calc - c);

    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
end
elapsedTime = toc;

a=min(T);
b=max(T);
c=(a+b)/2;
c_history = [];

figure
plot(T, p+100000);

xlabel('Temperatura[K]');
ylabel('cisnienie par[Pa]')
hold on;

plot([a, b], [f(a)+100000, f(b)+100000], 'Color', [0.5, 0.5, 0.5] , 'LineWidth', 0.25);

while abs(f(c)) > error
     
    c_calc = (a * f(b) - b * f(a)) / (f(b) - f(a)); 
    
    c = c + beta * (c_calc - c);
   
    c_history = [c_history, c];

    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    plot([a, b], [f(a)+100000, f(b)+100000],  'Color', [0.5, 0.5, 0.5] , 'LineWidth', 0.25);
end

title(sprintf('explicitConvergence - Czas: %.5f s, liczba iteracji: %d', elapsedTime, length(c_history)));
plot(xlim, [100000 100000], '--k')
scatter(c_history, f(c_history)+100000 ,'+','r');
hold off;
end