% Bestimme Naeherungsloesung fuer ein Feder-Daempfer-Masse-System:
% m*y''(t) + D*y'(t) + k*y(t) = F_0*cos(omega*t-Delta)
% mit y(t=0)=y_0, y'(t=0)=ydot_0 und m,k,D = const

%% HAUPTFUNKTION
% =========================================================================
function AWP_Loeser()

clc
clear
format shortE
close all

% Modellparameter
m     = 1;% Masse
k     = 1;% Federhaerte
D     = 0.1;% Daempfung
omega = 0.2;% Anregungskreisfrequenz
Delta = 0;% Phasenwinkel
F_0   = 0;% Kraftamplitude
kappa_max = 0.6;
kappa_min = 0.3;

% Simulationsparameter
y_0   = [1,0];% Anfangswerte
t_0   = 0;% Startzeitpunkt
t_end = 10;% Endzeitpunkt

% Diskretisierungsparameter
N       = 200; % Anzahl Intervallschritte
h       = 0.1; % Schrittweite
t       = t_0;% Laufvariable
yh      = y_0;% Loesungsvariable
n_start = 2; % Startindex Zeititeration

% Loesungsvektor anlegen
y_numerisch      = zeros(3,N+1);    % erste Zeile = Stuetzstellen t, zweite Zeile = Loesungswerte y1, dritte Zeile = Loesungswerte y2
y_numerisch(:,1) = [t_0, y_0(1), y_0(2)];
analytisch_auflosung=1000;  %Auflösung für plot
y_analytisch=zeros(2,analytisch_auflosung);
y_anregung=zeros(1,analytisch_auflosung);

n=n_start;
% Zeit-Iterationsschleife .............................................
while t < t_end


    % Berechnung Loesungsinkrement
    [dy,kappa] = RK4(t,yh,m,k,D,F_0,omega,Delta,h);
    if kappa>kappa_max
        h=h/2;
        [dy,kappa] = RK4(t,yh,m,k,D,F_0,omega,Delta,h);
    end

    if (t+h)>t_end %Achte darauf, dass die endzeit eingehalten genau wird
        h=t_end-t;
        [dy,kappa] = RK4(t,yh,m,k,D,F_0,omega,Delta,h);
    end


    % Lauf- und Loesungsvariable inkrementieren
    t=t+h;
    yh= dy+yh;
    % Werte in Loesungsvektor speichern
    y_numerisch(:,n) = [t;yh'];

    if kappa<kappa_min
        h=2*h;
    end
    %Zähler erhöhen
    n=n+1;
end % .................................................................
%letzte Erhöhung entfernen
n=n-1;
% Analytische Werte berechnen
y_analytisch(1,:)=linspace(t_0,t_end,analytisch_auflosung);
y_analytisch(2,:)=y1_erzwungen_exakt(y_analytisch(1,:),k,D,m,omega,F_0,Delta,y_0);
%Anregung Werte
y_anregung=anregung(y_analytisch(1,:),F_0,omega,Delta);
% quadratischen Fehler berechnen und ausgeben
Err=0;
for i=n_start:(n) %n=n+1
    y_ref= y1_erzwungen_exakt(y_numerisch(1,i),k,D,m,omega,F_0,Delta,y_0);
    diff = y_ref-y_numerisch(2,i);
    Err = Err+diff^2;
end
Err=sqrt((1/(n-1))*Err)

% analytische und numerische Loesung plotten
plot(y_numerisch(1,:),y_numerisch(2,:),'b+')
hold on
plot(y_analytisch(1,:),y_analytisch(2,:),'r')
plot(y_analytisch(1,:),y_anregung,'g')
legend({'Numerisch','Analytisch','Anregung'})


end % Funktion AWP_Loeser()


%% UNTERFUNKTIONEN
% =========================================================================

% rechte Seite (RHS) der DGL ----------------------------------------------
function [RHS] = f(t,y,m,k,D,F_0,omega,Delta)
RHS=[y(2),(-1/m)*(D*y(2)+k*y(1)-anregung(t,F_0,omega,Delta))];
end

% Erregungskraft ----------------------------------------------------------
function [F] = anregung(t,F_0,omega,Delta)
F = F_0*cos(omega*t-Delta);
end

% Loesungsinkrement dy des klassischen Runge-Kutta-Verfahrens -------------
function [dy,kappa] = RK4(t,y,m,k,D,F_0,omega,Delta,h)
k1 = f(t,y,m,k,D,F_0,omega,Delta);
k2 = f(t+h*0.5,y+0.5*h*k1,m,k,D,F_0,omega,Delta);
k3 = f(t+h*0.5,y+0.5*h*k2,m,k,D,F_0,omega,Delta);
k4 = f(t+h,y+h*k3,m,k,D,F_0,omega,Delta);
kappa= h*abs((2/h)*(k3-k2)./(k2-k1));
kappa = max(kappa);%Man nehme das Maximum, da kappa besonders nach oben beschränkt ist 
dy = (h/6)*(k1+2*(k2+k3)+k4);
end

% analytische Loesung y(t) fuer erzwungene Schwingung ---------------------
function [output] = y1_erzwungen_exakt(t,k,D,m,omega,F_0,Delta,y_0)
gamma=D/m;
omega_0=sqrt(k/m);
omega_d=sqrt(omega_0^2-(gamma*0.5)^2);
Phi= atan2(-gamma*omega,omega_0^2-omega^2);
A1=1/(sqrt(m^2 *((omega^2-omega_0^2)^2 +(gamma*omega)^2)));
phi_d= atan((1/omega_d)* ( (omega*A1*F_0*sin(-Delta+Phi)+y_0(2))/(A1*F_0*cos(-Delta+Phi)-y_0(1)) - gamma*0.5 ));
A2=(y_0(1)-A1*F_0*cos(-Delta+Phi))/cos(phi_d);
output = y_0(1)*exp(-0.5*gamma*t).*(cos(omega_d*t)+gamma*sin(omega_d*t)/(2*omega_d));
%output= A1*F_0*cos(omega*t-Delta+Phi)+A2*exp(-0.5*gamma*t).*cos(omega_d*t+phi_d);
end