function Cylinder_Model_Code2
%% =========================================================================
%  NAME : LAIMY HAINANA
%  DATE : 23/08/2026
%
%  CYLINDER GRAVITY - DENSITY-CONTRAST INVERSION
%  -------------------------------------------
%   1) Build a synthetic ("observed") gravity anomaly of a buried horizontal cylinder.
%   2) Build an initial-guess model with a wrong density contrast.
%   3) Invert for the density contrast with an overdetermined linear
%      least-squares system  (A'A) dRho = A'd  written explicitly as matrices.
%      Depth and radius are assumed known. The problem is LINEAR in DeltaRho,
%      so one step is exact for noise-free data.
%   4) Repeat with 0.15 % and 3 % Gaussian noise (noise level is a
%      percentage of the peak anomaly).
%% =========================================================================

clc; close all;
rng(1);                                % fixed seed -> repeatable noise

%% ---------------- CONSTANTS AND TRUE MODEL ---------------------------------
G            = 6.674e-11;               % Gravitational constant [m^3 kg^-1 s^-2]
deltaRhoTrue = 1000;                    % True density contrast [kg/m^3]
z_true       = 200;                     % True depth to cylinder centre [m]
R_true       = 50;                      % True cylinder radius [m]
x            = linspace(-500,500,100)'; % Observation points (column vector) [m]

M0 = cylinderMass(R_true, deltaRhoTrue);  % True excess mass [kg]

%% ---------------- FORWARD MODEL = NOISE FREE OBSERVED DATA --------------------------------------------
g_Observed_cylinder = cylinderGravity(G, M0, z_true, x);   % Noise-free anomaly [mGal]

noiseLevel       = 0;                                  % 0 % noise at this stage
gObserved_cylinder = addNoise(g_Observed_cylinder, noiseLevel);

%figure 1
figure;
plot(x, gObserved_cylinder,  'k*', 'LineWidth', 2); axis tight; hold on;
plot(x, g_Observed_cylinder, 'b-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed data (0% noise) vs Calculated Cylinder Gravity Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed (no noise)','Calculated');

%% ---------------- INITIAL-GUESS MODEL --------------------------------------
deltaRhoInitial = 900;                                 % Wrong starting density contrast [kg/m^3]
M1 = cylinderMass(R_true, deltaRhoInitial);              % Initial-guess mass [kg]
g_Calculated_cylinder = cylinderGravity(G, M1, z_true, x); % Response of initial model [mGal]

misfit_cylinder = g_Observed_cylinder - g_Calculated_cylinder;
MSE_cylinder    = calculateMSE(misfit_cylinder);

%figure 2
figure;
plot(x, g_Observed_cylinder,   'k-', 'LineWidth', 3); axis tight; hold on;
plot(x, g_Calculated_cylinder, 'g-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed vs Initial-Guess Cylinder Gravity Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed','Initial-guess model');

%% ---------------- DENSITY INVERSION (LINEAR LEAST SQUARES) -----------------
%  g(x) = G*M*z/(x^2+z^2) * 1e5 ,   M = (4/3)*pi*R^2*DeltaRho
%  A(x) = dg/dDeltaRho = dg/dDeltaRho = 2*G*pi*R^2*z / (x^2+z^2) * 1e5  ,   V = (4/3)*pi*R^2

V = (4/3)*pi*R_true.^2;

A   = sensitivityMatrixCylinder(G, V, z_true, x);        % 100 x 1
d   = g_Observed_cylinder - g_Calculated_cylinder;         % Residual vector, 100 x 1
AT  = A';                                              % 1 x 100
ATA = AT*A;                                            % 1 x 1
ATd = AT*d;                                            % 1 x 1
dRho = inv(ATA) * ATd;                                    % Density correction [kg/m^3]

deltaRho_estimated = deltaRhoInitial + dRho;
M_estimated        = cylinderMass(R_true, deltaRho_estimated);
g_estimated_cylinder = cylinderGravity(G, M_estimated, z_true, x);

Estimated_misfit = g_Observed_cylinder - g_estimated_cylinder;
MSE_estimated    = calculateMSE(Estimated_misfit);

%figure 3
figure;
plot(x, g_Observed_cylinder,  'k-',  'LineWidth', 3); axis tight; hold on;
plot(x, g_estimated_cylinder, 'bo',  'LineWidth', 2);
grid on; box on; hold off;
title('Cylinder Density Gravity Anomaly: Observed vs Inverted');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed','Inverted model');

%figure 4
figure;
plot(x, g_Observed_cylinder,   'ko', 'MarkerSize', 5); axis tight; hold on;
plot(x, g_Calculated_cylinder, 'b--','LineWidth', 2);
plot(x, g_estimated_cylinder,  'r-', 'LineWidth', 2);
hold off; grid on;
xlabel('Horizontal Distance (m)'); ylabel('Gravity Anomaly (mGal)');
title('Cylinder: Observed, Initial-Guess and Inverted Models');
legend('Observed','Initial guess','Inverted');

fprintf('\n');
fprintf('==============================================\n');
fprintf('      CYLINDER GRAVITY INVERSION RESULTS\n');
fprintf('==============================================\n');
fprintf('Number of observations        = %d\n', length(x));
fprintf('Size of A                     = %dx%d\n', size(A,1),   size(A,2));
fprintf('Size of AT                    = %dx%d\n', size(AT,1),  size(AT,2));
fprintf('Size of ATA                   = %dx%d\n', size(ATA,1), size(ATA,2));
fprintf('Size of ATd                   = %dx%d\n', size(ATd,1), size(ATd,2));
fprintf('True density contrast         = %.2f kg/m^3\n', deltaRhoTrue);
fprintf('Initial density contrast      = %.2f kg/m^3\n', deltaRhoInitial);
fprintf('Estimated density contrast    = %.2f kg/m^3\n', deltaRho_estimated);
fprintf('True mass                     = %.3e kg\n', M0);
fprintf('Initial guess mass            = %.3e kg\n', M1);
fprintf('Estimated mass                = %.3e kg\n', M_estimated);
fprintf('MSE of initial model          = %.3e mGal^2\n', MSE_cylinder);
fprintf('MSE of estimated model        = %.3e mGal^2\n', MSE_estimated);
fprintf('Maximum residual magnitude    = %.3e mGal\n', max(abs(Estimated_misfit)));
fprintf('==============================================\n');

%% ---------------- NOISY DATA: 0.15 % AND 3 % -------------------------------
noiseLevels = [0.0015, 0.03];                          % Fractions of peak anomaly
figNo = 5;
for N = 1:numel(noiseLevels)
    pct  = 100*noiseLevels(N);
    gObs = addNoise(g_Observed_cylinder, noiseLevels(N));% Noisy observed data

    d_noisy    = gObs - g_Calculated_cylinder;           % Residual for noisy data
    dRho_n     = inv(ATA) * (AT*d_noisy);                   % Density correction
    deltaRho_n = deltaRhoInitial + dRho_n;
    M_n        = cylinderMass(R_true, deltaRho_n);
    g_est_n    = cylinderGravity(G, M_n, z_true, x);
    MSE_n      = calculateMSE(gObs - g_est_n);

    % figure 5 / 7
    figure;
    plot(x, gObs, 'mo', 'LineWidth', 2); axis tight; hold on;
    plot(x, g_Observed_cylinder, 'k-', 'LineWidth', 2);
    grid on; box on; hold off;
    title(sprintf('Observed data (%.2f%% noise) vs Noise-free Cylinder Anomaly', pct));
    xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
    legend('Observed noisy','Noise-free');

    % figure 6 / 8
    figure;
    plot(x, gObs, 'ko', 'MarkerSize', 5); axis tight; hold on;
    plot(x, g_Calculated_cylinder, 'b--', 'LineWidth', 2);
    plot(x, g_est_n,             'r-',  'LineWidth', 2);
    hold off; grid on;
    xlabel('Horizontal Distance (m)'); ylabel('Gravity Anomaly (mGal)');
    title(sprintf('Cylinder Inversion with %.2f%% added noise', pct));
    legend('Observed noisy','Initial guess','Inverted','Location','best');
    figNo = figNo + 2;

    fprintf('\n--- Noise level = %.2f%% ---\n', pct);
    fprintf('Estimated density contrast = %.2f kg/m^3  (error %.2f %%)\n', ...
            deltaRho_n, 100*abs(deltaRho_n-deltaRhoTrue)/deltaRhoTrue);
    fprintf('MSE after inversion        = %.3e mGal^2\n', MSE_n);
end
end

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function M = cylinderMass(R, deltaRho)                   % Excess mass of a cylinder [kg]
    M = (4/3)*pi*R.^2 .* deltaRho;
end

function g = cylinderGravity(G, M, z, x)                 % Vertical attraction, centre at x = 0 [mGal]
    g = ((2 * G * M * z) ./(x.^2 + z.^2)) * 1e5;       % 1e5 converts m/s^2 -> mGal
end

function MSE = calculateMSE(misfit)                    % Mean-squared error
    MSE = mean(misfit(:).^2);
end

function A = sensitivityMatrixCylinder(G, V, z, x)       % dg/dDeltaRho  (V = cylinder volume)
    A = ((2*G*V*z) ./ (x.^2 + z.^2))* 1e5;  
end

function gN = addNoise(g, level)                       % Gaussian noise = level * peak |g|
    gN = g + level*max(abs(g))*randn(size(g));
end