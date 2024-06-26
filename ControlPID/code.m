clear
close all
clc
%Paso1: Cargar los datos reales de la respuesta del proceso y asignarlos a
%una variable
%load("/Users/josevarela/Library/Mobile Documents/com~apple~CloudDocs/TEC/Semestre 5/Análisis de sistemas de control/Reto/IdentPA.txt")
IdentPA=load("datosFinales.txt");

t=IdentPA(:,1);
u=IdentPA(:,2);
y=IdentPA(:,3);

%Paso2: Especificar la magnitud del escalon
opt = stepDataOptions('StepAmplitude',50);

%Paso3: Graficar la respuesta real y la modelacion 
figure
    subplot(2,1,1);
    plot(t,y,'r-','LineWidth',2)
    xlim([0 0.6]);
    ylim([0 40]);
    ax=gca;
    ax.FontSize = 15;
    ax.XTick = [0:0.05:0.6];
    ax.YTick = [0:10:40];

    title('Respuesta del Proceso Real','FontName','Arial','FontSize',18);
    ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
    legend('Variable de Proceso','Location', 'southeast');
    grid on;
        
    subplot(2,1,2);
    plot(t,u,'b--','LineWidth',2)
    xlim([0 0.6]);
    ylim([0 60]);
    ax=gca;
    ax.FontSize = 15;
    ax.XTick = [0:0.05:0.6];
    ax.YTick = [0:10:60];

    xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
    ylabel('Leds encendidos [%]','FontSize',15,'FontWeight','bold');
    legend('Manipulación','Location', 'southeast');
    grid on;    
    
figure
    plot(t,y,'r-','LineWidth',2)
    xlim([0 0.6]);
    ylim([0 30]);
    ax=gca;
    ax.FontSize = 15;
    ax.XTick = [0:0.05:0.6];
    ax.YTick = [0:5:25];

    title('Respuesta del Proceso Real','FontName','Arial','FontSize',18);
    xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
    ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
    legend('Variable de Proceso','Location', 'southeast');
    grid on;
%% Parámetros de Identificación
i=80; %  Intersección de la tangente (Índice Abscisa)
Mg = 50; % Magnitud de la entrada escalón
valFin = 22; % Valor final respuesta proceso
valInt = 0; % Valor incial respuesta proceso

%% Metodo Ziegler-Nichols
% Cálculo de tangente
dy=diff(y)./diff(t);
tang=(t-t(i))*dy(i)+y(i);
lim = ones(size(t))*valFin;
org = zeros(size(t));

% Cálculo de valores de los parámetros
K_zn = (valFin-valInt)/Mg;
Teta_zn = -((t(1)-t(i))*dy(i)+y(i))/dy(i);
Tao_zn = (valFin -((t(1)-t(i))*dy(i)+y(i)))/dy(i) - Teta_zn;

% Visualización del método
figure
hold on

plot(t,y, 'k', 'LineWidth', 2.5); 
plot(t,tang, 'r', 'LineWidth', 1.5);
plot(t,lim, 'b--', 'LineWidth',1);
plot(t,org, 'b--', 'LineWidth',1);

scatter(Tao_zn + Teta_zn, valFin,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)',  Tao_zn + Teta_zn, valFin);
text(Tao_zn + Teta_zn + 0.01,valFin+1,textString,'FontSize', 15)
scatter(Teta_zn, valInt,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)', Teta_zn, valInt);
text(Teta_zn+ 0.01,valInt+1,textString,'FontSize', 15)

xlim([0 0.6]);
ylim([-1 30]);
title('Método de Ziegler-Nichols','FontName','Arial','FontSize',18);
xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
hold off


fprintf('Método de Ziegler-Nichols: k = %g,  τ= %g, θ = %g\n', K_zn, Tao_zn, Teta_zn);
G_zn = tf(K_zn, [Tao_zn, 1], "InputDelay", Teta_zn)

 %% Metodo Miller

% Cálculo de tangente
dy=diff(y)./diff(t);
tang=(t-t(i))*dy(i)+y(i);
org = zeros(size(t));


%Encontrar cordenada t63
val63 = valFin*(1-exp(-1));
[~, i63] = min(abs(val63-y));
t63 = t(i63);

% Cálculo de modelo
K_mi = (valFin-valInt)/Mg;
Teta_mi = -((t(1)-t(i))*dy(i)+y(i))/dy(i);
Tao_mi = t63 - Teta_mi;


figure
hold on

plot(t,y, 'k', 'LineWidth', 2.5); 
plot(t,tang, 'r', 'LineWidth', 1.5);
plot(t,org, 'b--', 'LineWidth',1);
plot([0 t63], [val63 val63], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');
plot([t63 t63], [val63 0], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');


scatter(Teta_mi, valInt,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)', Teta_mi, valInt);
text(Teta_mi+ 0.01,valInt+1,textString,'FontSize', 15)

scatter(t63, val63,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)',  t63, val63);
text(t63 + 0.01,val63+1,textString,'FontSize', 15)

xlim([0 0.6]);
ylim([-1 30]);
title('Método de Miller','FontName','Arial','FontSize',18);
xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
hold off

fprintf('Método de Millers: k = %g,  τ= %g, θ = %g\n', K_mi, Tao_mi, Teta_mi);
G_mi = tf(K_mi, [Tao_mi, 1], "InputDelay", Teta_mi)

%% Metodo analítico


%Encontrar cordenada t63 y t28
val63 = valFin*(1-exp(-1));
[~, i63] = min(abs(val63-y));
t63 = t(i63);
val28 = valFin*(1-exp(-1/3));
[~, i28] = min(abs(val28-y));
t28 = t(i28);

% Cálculo de modelo
K_an = (valFin-valInt)/Mg;
Tao_an= 3/2 * (t63-t28);
Teta_an= t63 - Tao_an;


figure
hold on

plot(t,y, 'k', 'LineWidth', 2.5); 
plot([0 t63], [val63 val63], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');
plot([t63 t63], [val63 0], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');
plot([0 t28], [val28 val28], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');
plot([t28 t28], [val28 0], 'Color','blue', 'LineWidth',1 ,'LineStyle','--');

scatter(t28, val28,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)',  t28, val28);
text(t28 + 0.01,val28+1,textString,'FontSize', 15)


scatter(t63, val63,'MarkerEdgeColor',[0 .5 .5],...
                        'MarkerFaceColor',[0 .7 .7],...
                        'LineWidth',1.5)
textString = sprintf('(%g, %g)',  t63, val63);
text(t63 + 0.01,val63+1,textString,'FontSize', 15)

xlim([0 0.6]);
ylim([0 30]);
title('Método de Analítico','FontName','Arial','FontSize',18);
xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
hold off

fprintf('Método de Analítico: k = %g,  τ= %g, θ = %g\n', K_an, Tao_an, Teta_an);
G_an = tf(K_an, [Tao_an, 1], "InputDelay", Teta_an)

%% Comparación de métodos
opt = stepDataOptions('StepAmplitude',Mg);

% Respuestas a entrada escalón:
Y_zn = step(G_zn,t,opt);
Y_mi=step(G_mi,t,opt);
Y_an=step(G_an,t,opt);

% Graficar respuestas
figure
hold on
plot(t,y,'r-','LineWidth',2)
plot(t,Y_zn,'g-','LineWidth',2)
plot(t,Y_mi,'b-','LineWidth',2)
plot(t,Y_an,'y-','LineWidth',2)
xlim([0 0.6]);
ylim([0 30]);
ax=gca;
ax.FontSize = 15;
ax.XTick = 0:0.05:0.6;
ax.YTick = 0:5:25;
title('Respuesta del Proceso Real','FontName','Arial','FontSize',18);
xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
legend('Real','Ziegler&Nichols','Miller','Analitico','Location', 'southeast');
grid on;
hold off

% Cálculo de suma de errores al cuadrado
SSE_zn=sum((y-Y_zn).^2);
SSE_mill=sum((y-Y_mi).^2);
SSE_an=sum((y-Y_an).^2);
fprintf(['Cálculo de suma de errores al cuadrado:\n' ...
    'Método de Ziegler-Nichols: SSE = %g\n' ...
    'Método de Millers:  SSE= %g\n' ...
    'Método de Analítico: SSE = %g\n'], SSE_zn, SSE_mill, SSE_an);
%% Selección de mejor identificación 
minSSE = min([SSE_zn, SSE_mill, SSE_an]);
if minSSE == SSE_zn
    fprintf("El Método de Ziegler-Nichols tiene el menor SSE.\n")
    K_bt = K_zn;
    Tao_bt = Tao_zn;
    Teta_bt = Teta_zn;
    met = "Ziegler-Nichols";
elseif minSSE == SSE_mill
    fprintf("El Método de Millers tiene el menor SSE.\n")
    K_bt = K_mi;
    Tao_bt = Tao_mi;
    Teta_bt = Teta_mi;
    met = "Millers";
else 
    fprintf("El Método de Analítico tiene el menor SSE.\n")
    K_bt = K_an;
    Tao_bt = Tao_an;
    Teta_bt = Teta_an;
    met = "Analítico";
end 
G_1 =  tf(K_bt, [Tao_bt, 1], "InputDelay", Teta_bt);

%% Aproximación de Padé
delayPade = tf([-Teta_bt/2 1], [Teta_bt/2 1]);
Gp_nd = tf([0 K_bt], [Tao_bt 1]);
Gp_2 = Gp_nd*delayPade;
% Se muestra comparación de mejor identificación y su aproximación de Padé
figure 
hold on
Y_pd=step(Gp_2,t,opt);
Y_bt=step(G_1,t,opt);
plot(t,Y_pd,'k-','LineWidth',2)
plot(t,Y_bt,'b-','LineWidth',2)
xlim([0 0.6]);
ylim([-3 27]);
ax=gca;
ax.FontSize = 15;
ax.XTick = 0:0.05:0.6;
ax.YTick = 0:5:25;
title('Respuesta del Proceso Real','FontName','Arial','FontSize',18);
xlabel('Tiempo [seg]','FontSize',15,'FontWeight','bold');
ylabel('Luminosidad [%]','FontSize',15,'FontWeight','bold');
legend('Aproximacón de Padé', met,'Location', 'southeast');
grid on;
hold off
%% Root Locus de aproximación de Padé
figure
rlocus(Gp_2);
title(' Root Locus de aproximación de Padé','FontName','Arial','FontSize',15);
%% Ejemplos de polos de distintos comportamientos 
figure
hold on
rlocus(Gp_2);
scatter(-34, 0, 100,'b','*','LineWidth',1.5);
scatter([-50 -20.4], [0 0], 100,'g','x','LineWidth',1.5);
scatter([-27 -27], [34.4 -34.4], 100,'r','x','LineWidth',1.5);
legend('Root Locus Aproximacón de Padé', 'Ejemplo de polos de sistema Críticamente amortiguado', 'Ejemplo de polos de sistema Sobreamortiguado', 'Ejemplo de polos de sistema Subamortiguado'  ,'Location', 'southeast');
title('Ejemplos de polos de distintos comportamientos ','FontName','Arial','FontSize',15);
hold off
%% Respuesta al escalon del sistema G2 retroalimentado 
G2_fb = feedback(Gp_2,1);
figure
step(G2_fb);
title('Respuesta al escalon del sistema G2 retroalimentado ','FontName','Arial','FontSize',15);
%% Respuesta al escalon del sistema G2 retroalimentado con ganacia k_U
figure
k_u = 0.07675/0.008030;
G2_fb_ku = feedback(k_u*Gp_2,1);
step(G2_fb_ku);
xlim([0 0.6]);
ylim([-1 3]);
title('Respuesta al escalon del sistema G2 retroalimentado con ganacia k_U','FontName','Arial','FontSize',15);
%% Coeficientes de controladores P, PI y PID con Primer método de Ziegler-Nichols
% Primer método de Ziegler-Nichols
C_ZN1_P = Tao_bt/(K_bt*Teta_bt);
C_ZN1_PI = [0.9*Tao_bt/(K_bt*Teta_bt) 3.3*Teta_bt];
C_ZN1_PID = [1.2*Tao_bt/(K_bt*Teta_bt) 2*Teta_bt 0.5*Teta_bt];
%% Coeficientes de controladores PID con criterios integrales
% PID ISE
C_ISE_KC = 1.474/K_bt*(Teta_bt/Tao_bt)^-0.97;
C_ISE_TI = Tao_bt/1.115*(Teta_bt/Tao_bt)^-0.753;
C_ISE_TD = 0.550*Tao_bt*(Teta_bt/Tao_bt)^0.948;
C_ISE = [C_ISE_KC C_ISE_TI C_ISE_TD];

% PID IAE cambio de referencia
C_IAE_REFE_KC = 1.086/K_bt*(Teta_bt/Tao_bt)^-0.869;
C_IAE_REFE_TI = Tao_bt/(0.74-0.130*(Teta_bt/Tao_bt));
C_IAE_REFE_TD = 0.348*(Tao_bt)*(Teta_bt/Tao_bt)^0.914;
C_IAE_REFE = [C_IAE_REFE_KC C_IAE_REFE_TI C_IAE_REFE_TD];

% PID IAE Rechazo de perturbaciones
C_IAE_RECHA_KC = 1.435/K_bt*(Teta_bt/Tao_bt)^-0.921;
C_IAE_RECHA_TI = Tao_bt/0.878*(Teta_bt/Tao_bt)^0.749;
C_IAE_RECHA_TD = 0.482*Tao_bt*(Teta_bt/Tao_bt)^1.137;
C_IAE_RECHA = [C_IAE_RECHA_KC C_IAE_RECHA_TI C_IAE_RECHA_TD];

% PID ITAE cambio de referencia
C_ITAE_REFE_KC = 0.965/K_bt*(Teta_bt/Tao_bt)^-0.855;
C_ITAE_REFE_TI = Tao_bt/(0.796-0.147*(Teta_bt/Tao_bt));
C_ITAE_REFE_TD = 0.308*Tao_bt*(Teta_bt/Tao_bt)^0.992;
C_ITAE_REFE = [C_ITAE_REFE_KC C_ITAE_REFE_TI C_IAE_REFE_TD];

% PID ITAE Rechazo de perturbaciones
C_ITAE_RECHA_KC = 1.357/K_bt*(Teta_bt/Tao_bt)^-0.947;
C_ITAE_RECHA_TI = Tao_bt/0.842*(Teta_bt/Tao_bt)^0.738;
C_ITAE_RECHA_TD = 0.381*Tao_bt*(Teta_bt/Tao_bt)^0.995;
C_ITAE_RECHA = [C_ITAE_RECHA_KC C_ITAE_RECHA_TI C_ITAE_RECHA_TD];
%% Controlador con PID con el segundo método de Ziegler Nichols
P_u = 0.148-0.0585;
C_ZN2_PID = [k_u/1.7 P_u/2 P_u/8];
%% Respuesta al escalon del sistema Gp_pd con compensador y retroalimentado 
Gc = tf([1 -1/0.044*0.44/0.008030], [1 -0.44/0.008030]);
figure
rlocus(Gc*Gp_2);
title(' Root Locus de aproximación de Padé con compensador','FontName','Arial','FontSize',15);
k_rl = 6;
k_an = 5.104;
G2_cp_rl = feedback(k_rl*Gc*Gp_2,1);
G2_cp_an = feedback(k_an*Gc*Gp_2,1);
figure
step(G2_cp_rl)
title('Respuesta al escalon del sistema G2 con compensador y ganancia','FontName','Arial','FontSize',15);
legend('Respuesta del sistema con k=6' ,'Location', 'southeast');
figure
step(G2_cp_an)
legend('Respuesta del sistema con k=5.104' ,'Location', 'southeast');
title('Respuesta al escalon del sistema G2 con compensador y ganancia','FontName','Arial','FontSize',15);

%% Espacio de estados de función G2 (aprx. de padé)
[Gp_2_num, Gp_2_den] = tfdata(Gp_2);
[A,B,C, D] = tf2ss(Gp_2_num{1}, Gp_2_den{1});

%Revisamos controlabilidad y observabilidad
rank(ctrb(A, B))
det(ctrb(A, B))
rank(obsv(A, C))
det(obsv(A, C))

%% Parámetros del control;
ts = 0.3; %Tiempo de estabilización 
os = 0.1; % Over shoot
ro = 10; % velocidad del observador
%% Cálculo de polos
sigma_d = 4/ts;
z = -log(os)/sqrt(pi^2+log(os)^2);
theta = acos(z);
wd = sigma_d*tan(theta);
%Se proponen polos del sistema
cpoles = [-sigma_d+wd*1i, -sigma_d-wd*1i];
%Se proponen polos del observador
opoles = cpoles*ro;
%% Diseño de controlador y observador 
%Se calculan las ganancias de retroalimentación
K = place(A, B, cpoles);
%Se calculan las ganancias del observador
L = place(A', C', opoles)';
%Se calcula la ganancia Gpf
At = [A-B*K B*K; zeros(size(A)) A-L*C];
Bt= [B; zeros(size(B))];
Ct=[C zeros(size(C))];
sys = ss(At, Bt, Ct, D);
% Se encuentra el valor del estado final del sistema para eliminar el ess
valFin = evalfr(sys,0);
%Se multiplica B por la ganancia Gpf para eliminar ess
B2=B/valFin;
%Se calcula el sistema completo
At = [A-B*K B*K; zeros(size(A)) A-L*C];
Bt= [B2; zeros(size(B))];
Ct=[C zeros(size(C))];
sys= ss(At, Bt, Ct, D);
%Se Revisa el sistema como regulador con condiciones iniciales
x0 = [-0.1 0.3];
%Tiempo para simular el sistema (0 a 3 s)
t=0:.0001:3;
%Se obtiene la respuesta del sistema con condiciones iniciales
[y,t,x] = lsim(sys, ones(size(t)), t, [x0 x0]);
%Gráfica de salida
figure
plot(t, y);
title('Respuesta al escalon del sistema como regulador de condiciones iniciales','FontName','Arial','FontSize',15);
%Se obtiene la respuesta del sistema con entrada de escalon
figure
step(sys)
title('Respuesta al escalon del sistema con controlador y observador','FontName','Arial','FontSize',15);
%% Referencias 
% Para la algunas secciones se utilizaron codigos creados por 
%       Dr. Carlos Sotelo - Dr. David Sotelo
%       Mto. Salvador Leal
%       Análisis de Sistemas de Control
%       Aug-Dec 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%