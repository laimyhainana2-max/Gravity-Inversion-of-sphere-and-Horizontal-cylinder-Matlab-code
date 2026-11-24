function Sphere_Model_Code3_depth2
%% =========================================================================
%  NAME : LAIMY HAINANA
%  DATE : 23/08/2026
%
%  SPHERE GRAVITY - DEPTH INVERSION (NONLINEAR, GAUSS-NEWTON)
%  ----------------------------------------------------------
%   1) Build a synthetic ("observed") gravity anomaly of a buried sphere.
%   2) Start from a wrong depth (z_initial). Density contrast and radius are
%      assumed known, so the sphere mass M0 is known.
%   3) Gravity is NONLINEAR in depth, so the linearised system
%          (A'A) dz = A'd ,   A = dg/dz
%      is solved repeatedly (iterations), re-computing g and A each time
%      with the updated depth, until the depth stops changing.
%   4) Repeat with 0.15 % and 3 % Gaussian noise.
%% =========================================================================

clc; close all;
rng(1);                                % fixed seed -> repeatable noise

%% ---------------- CONSTANTS AND TRUE MODEL ---------------------------------
G            = 6.674e-11;
deltaRhoTrue = 1000;                    % [kg/m^3]
z_true       = 200;                     % True depth to sphere centre [m]
R_true       = 50;                      % True radius [m]
x            = linspace(-500,500,100)'; % Observation points (column) [m]

M0 = sphereMass(R_true, deltaRhoTrue);  % Known mass [kg]

%% ---------------- FORWARD MODEL --------------------------------------------
g_Observed_sphere = sphereGravity(G, M0, z_true, x);
gObserved_sphere  = addNoise(g_Observed_sphere, 0);    % 0 % noise

%figure 1
figure;
plot(x, gObserved_sphere,  'k*', 'LineWidth', 2); axis tight; hold on;
plot(x, g_Observed_sphere, 'b-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed data (0% noise) vs Calculated Sphere Gravity Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed (no noise)','Calculated');

%% ---------------- INITIAL-GUESS MODEL (DEPTH) ------------------------------
z_initial = 180;                                       % Wrong starting depth [m]
g_Calculated_sphere = sphereGravity(G, M0, z_initial, x);

misfit_sphere = g_Observed_sphere - g_Calculated_sphere;
MSE_sphere    = calculateMSE(misfit_sphere);

%figure 2
figure;
plot(x, g_Observed_sphere,   'k-', 'LineWidth', 3); axis tight; hold on;
plot(x, g_Calculated_sphere, 'g-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed vs Initial-Guess Sphere Gravity Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed','Initial-guess model (z = 180 m)');

%% ---------------- DEPTH INVERSION (GAUSS-NEWTON) ---------------------------
%  g(z) = G*M*z/(x^2+z^2)^(3/2)*1e5
%  A(z) = dg/dz = G*M*(x^2 - 2 z^2)/(x^2+z^2)^(5/2)*1e5
maxIter = 20;  tol = 1e-6;                             % tol on |dz| [m]
[z_estimated, zHist, A, AT, ATA, ATd] = gaussNewtonDepth(G, M0, z_initial, x, gObserved_sphere, maxIter, tol);

g_estimated_sphere = sphereGravity(G, M0, z_estimated, x);
Estimated_misfit   = g_Observed_sphere - g_estimated_sphere;
MSE_estimated      = calculateMSE(Estimated_misfit);

%figure 3
figure;
plot(x, g_Observed_sphere,  'k-', 'LineWidth', 3); axis tight; hold on;
plot(x, g_estimated_sphere, 'bo', 'LineWidth', 2);
grid on; box on; hold off;
title('Sphere Gravity Anomaly: Observed vs Inverted');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed','Inverted model');

%figure 4
figure;
plot(x, g_Observed_sphere,   'ko', 'MarkerSize', 5); axis tight; hold on;
plot(x, g_Calculated_sphere, 'b--','LineWidth', 2);
plot(x, g_estimated_sphere,  'r-', 'LineWidth', 2);
hold off; grid on;
xlabel('Horizontal Distance (m)'); ylabel('Gravity Anomaly (mGal)');
title('Sphere: Observed, Initial-Guess and Inverted Models');
legend('Observed','Initial guess','Inverted');

fprintf('\n');
fprintf('==============================================\n');
fprintf('    SPHERE GRAVITY DEPTH INVERSION RESULTS\n');
fprintf('==============================================\n');
fprintf('Number of observations        = %d\n', length(x));
fprintf('Size of A                     = %dx%d\n', size(A,1),   size(A,2));
fprintf('Size of AT                    = %dx%d\n', size(AT,1),  size(AT,2));
fprintf('Size of ATA                   = %dx%d\n', size(ATA,1), size(ATA,2));
fprintf('Size of ATd                   = %dx%d\n', size(ATd,1), size(ATd,2));
fprintf('True depth                    = %.2f m\n', z_true);
fprintf('Initial guess depth           = %.2f m\n', z_initial);
fprintf('Estimated depth               = %.2f m\n', z_estimated);
fprintf('Iterations used               = %d\n', numel(zHist)-1);
fprintf('Known mass used               = %.3e kg\n', M0);
fprintf('MSE of initial model          = %.3e mGal^2\n', MSE_sphere);
fprintf('MSE of estimated model        = %.3e mGal^2\n', MSE_estimated);
fprintf('Maximum residual magnitude    = %.3e mGal\n', max(abs(Estimated_misfit)));
fprintf('==============================================\n');

%% ---------------- NOISY DATA: 0.15 % AND 3 % -------------------------------
noiseLevels = [0.0015, 0.03];
for N = 1:numel(noiseLevels)
    pct  = 100*noiseLevels(N); %converts the fractional noise level into a percentage value 
    gObs = addNoise(g_Observed_sphere, noiseLevels(N));

    z_n     = gaussNewtonDepth(G, M0, z_initial, x, gObs, maxIter, tol);
    g_est_n = sphereGravity(G, M0, z_n, x);
    MSE_n   = calculateMSE(gObs - g_est_n);

    % figure 5 / 7
    figure;
    plot(x, gObs, 'mo', 'LineWidth', 2); axis tight; hold on;
    plot(x, g_Observed_sphere, 'k-', 'LineWidth', 2);
    grid on; box on; hold off;
    title(sprintf('Observed data (%.2f%% noise) vs Noise-free Sphere Anomaly', pct));
    xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
    legend('Observed noisy','Noise-free');

    % figure 6 / 8
    figure;
    plot(x, gObs, 'ko', 'MarkerSize', 5); axis tight; hold on;
    plot(x, g_Calculated_sphere, 'b--', 'LineWidth', 2);
    plot(x, g_est_n,             'r-',  'LineWidth', 2);
    hold off; grid on;
    xlabel('Horizontal Distance (m)'); ylabel('Gravity Anomaly (mGal)');
    title(sprintf('Sphere Depth Inversion with %.2f%% added noise', pct));
    legend('Observed noisy','Initial guess','Inverted','Location','best');

    fprintf('\n--- Noise level = %.2f%% ---\n', pct);
    fprintf('Estimated depth = %.2f m  (error %.2f %%)\n', z_n, 100*abs(z_n-z_true)/z_true);
    fprintf('MSE after inversion = %.3e mGal^2\n', MSE_n);
    
end
end

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function M = sphereMass(R, deltaRho)
    M = (4/3)*pi*R.^3 .* deltaRho;
end

function g = sphereGravity(G, M, z, x)
    g = (G*M*z) ./ ((x.^2 + z.^2).^(3/2)) * 1e5;       % [mGal]
end

function MSE = calculateMSE(misfit)
    MSE = mean(misfit(:).^2);
end

function A = sensitivityMatrixSphere(G, M, z, x)       % dg/dz  (depth Jacobian)
    A = G*M*(x.^2 - 2*z^2) ./ ((x.^2 + z^2).^(5/2)) * 1e5;
end

function gN = addNoise(g, level)
    gN = g + level*max(abs(g))*randn(size(g));
end

function [z, zHist, A, AT, ATA, ATd] = gaussNewtonDepth(G, M, z0, x, gObs, maxIter, tol)
    z = z0;  zHist = z0; %zHist==>store the estimated depth at every single iteration of the Gauss-Newton loop.
    for it = 1:maxIter
        d   = gObs(:) - sphereGravity(G, M, z, x(:));  % Residual with current depth
        A   = sensitivityMatrixSphere(G, M, z, x(:));  % Jacobian at current depth
        AT  = A';  ATA = AT*A;  ATd = AT*d;
        dz  = inv(ATA) * ATd;                               
        zNew = z + dz; % Depth correction [m]
        while zNew <= 0, dz = dz/2; zNew = z + dz; end % keep depth positive
        z = zNew;  zHist(end+1) = z;                   %#ok<AGROW>    => used to suppress the array growth inside a loop
        if abs(dz) < tol, break; end
    end
end