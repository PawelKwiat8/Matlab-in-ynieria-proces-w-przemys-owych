%% tworzenie Qwe
clc;clear;
fp=1000; % czestotliwość probkowania
Tp=1/fp;
t = 0:Tp:(0.2-Tp);
N = length(t);
f=50; %czestotliwośc pracy pompy [Hz]
T=1/f;
A=0.785; % w wolfram wyszlo tyle aby uzyskać 50ml na skok
Qwe = A*sin(2*pi*f*t);
Qwe = max(Qwe, 0); %ograniczenie przepłwyu do wartości dodatnich
%Qwe = abs(Qwe) ; %przy dwwustronnej

figure
plot(t,Qwe);
xlabel('s')
ylabel('l/s')

% sprawdzenie przepływu na jeden okres/skok 
fun = @(t) abs(A * sin(2*pi*f*t));
area_under_curve = integral(fun, 0, T/2);

%%fprintf('Pole powierzchni pod krzywą: %f\n', area_under_curve);

ro = 1000; %gestosc wody [kg/m^3]
g = 9.81;%przyspieszenie ziemskie [m/s^2]
R=8.314; %stała gazowa 
T=293.15; %stopnie [K]
r=1; %promien [m]
d=2; % dłogośc powietrznika [m]
A=pi*r^2; %srednica [m]
V=A*d; % objetośc gazu [m^3]
pt=15000; %cisnienie techniczne [Pa]
m = pt*V/R/T;
cv = 0.1;
z=zeros(1,N);
Qwy=zeros(1,N);
Pg=zeros(1,N);
P=zeros(1,N);

for i = 1:N
    if (i==1) 

        Pg(i)=m*R*T/(V-A*z(i));

        P(i)=ro*g*z(i)+Pg(i);
    else

        z(i)=z(i-1)+Tp*(Qwe(i-1)-Qwy(i-1))/A;

        Pg(i)=m*R*T/(V-A*z(i));

        P(i)=ro*g*z(i)+Pg(i);

        Qwy(i)=cv*sqrt(P(i)-pt);

    end
    
end

figure;
hold on;

plot(t, Qwe, 'b', 'LineWidth', 1.5);
plot(t, Qwy, 'r--', 'LineWidth', 1.5);

xlabel('Czas (s)');
ylabel('Przepływ (l/s)');
legend('Qwe', 'Qwy');
grid on;
hold off;
 