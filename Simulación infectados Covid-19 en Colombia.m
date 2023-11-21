% Parámetros
P = 20000000;%Población
Beta = 0.1479;%Tasa de transmisión
Gamma = 0.0125;%Tasa de Recuperación
b = 0.086;%Tasa de mortalidad
% Condiciones iniciales
I0 = 1046; %Numero inicial de Infectados (5 de abril del 2020)
T = 60; %Tiempo total de simulación en días
N = 1000; %Número de pasos de tiempo (Si N tiende a infinito las diferentes aproximaciones se aproximaran a la solucioj unica
deltat = T / N;% Tamaño del paso de tiempo
% Número de simulaciones estocásticas
num_simulations = 5; % Este valor se puede cambiar dependiendo de cuantas simulaciones estocasticas se quieran
% Inicialización de matriz para almacenar resultados de simulaciones estocásticas
I_stochastic = zeros(num_simulations, N + 1);
% Método de Euler para resolver la EDE varias veces
for k = 1:num_simulations
 I = zeros(1, N + 1);
 I(1) = I0;
 for j = 1:N
 dW = sqrt(deltat) * randn();
 I(j + 1) = I(j) + ((Beta / P) * (P - I(j)) * I(j) - (b + Gamma) * I(j)) * deltat + sqrt(((Beta / P) * (P - I(j)) * I(j) + (b + Gamma) * I(j))) * dW;
 % Ecuación SIS estocástica utilizando el método de Euler
 end
 I_stochastic(k, :) = I;
end
% Versión determinista
I_deterministic = zeros(1, N + 1);
I_deterministic(1) = I0;
for i = 1:N
 I_deterministic(i + 1) = I_deterministic(i) + ((Beta / P) * (P - I_deterministic(i)) * I_deterministic(i) - (b + Gamma) * I_deterministic(i)) * deltat;
end
% Gráfico de las soluciones
figure;
plot([0:deltat:T], I_stochastic', 'b-', 'LineWidth', 1); hold on
plot([0:deltat:T], I_deterministic, 'r-', 'LineWidth', 2); hold off
xlabel('Tiempo', 'FontSize', 12)
ylabel('Infectados', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
legend('Estocástico 1', 'Estocástico 2', 'Estocástico 3', 'Estocástico 4', 'Estocástico 5', 'Determinista')






