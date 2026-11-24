function Cylinder_model_code_depth
%% =========================================================================
%  NAME : LAIMY HAINANA
%  DATE : 23/08/2026
%
%  HORIZONTAL CYLINDER (2D) GRAVITY - DEPTH INVERSION (GAUSS-NEWTON)
%  -----------------------------------------------------------------
%   Density contrast and radius are assumed known (mass per unit length
%   lambda is known). Gravity is nonlinear in depth, so the linearised
%   system (A'A) dz = A'd is solved iteratively, with A re-computed at each
%   updated depth.
%% =========================================================================

clc; close all;
rng(1);                                 % fixed seed -> repeatable noise

%% ---------------- CONSTANTS AND TRUE MODEL ---------------------------------
G            = 6.674e-11;
deltaRhoTrue = 1000;                    % [kg/m^3]
z_true       = 200;                     % True depth to cylinder axis [m]
R_true       = 50;                      % True radius [m]
x            = linspace(-500,500,100)'; % Observation points (column) [m]

M0 = cylinderMass(R_true, deltaRhoTrue);% Known mass per unit length [kg/m]

%% ---------------- FORWARD MODEL --------------------------------------------

noiseLevel       = 0;                                  % 0 % noise 

g_Observed_cylinder = cylinderGravity(G, M0, z_true, x);
gObserved_cylinder  = addNoise(g_Observed_cylinder, noiseLevel);   % 0 % noise level

%figure 1
figure;
plot(x, gObserved_cylinder,  'k*', 'LineWidth', 2); axis tight; hold on;
plot(x, g_Observed_cylinder, 'b-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed data (0% noise) vs Calculated 2D Horizontal Cylinder Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed (no noise)','Calculated');

%% ---------------- INITIAL-GUESS MODEL (DEPTH) ------------------------------
z_initial = 180;                                          % Wrong starting depth [m]
g_Calculated_cylinder = cylinderGravity(G, M0, z_initial, x);

misfit_cylinder = g_Observed_cylinder - g_Calculated_cylinder;
MSE_cylinder    = calculateMSE(misfit_cylinder);

%figure 2
figure;
plot(x, g_Observed_cylinder,   'k-', 'LineWidth', 3); axis tight; hold on;
plot(x, g_Calculated_cylinder, 'g-', 'LineWidth', 2);
grid on; box on; hold off;
title('Observed vs Initial-Guess Cylinder Gravity Anomaly');
xlabel('Distance (m)'); ylabel('Gravity Anomaly (mGal)');
legend('Observed','Initial-guess model (z = 180 m)');

%% ---------------- DEPTH INVERSION (GAUSS-NEWTON) ---------------------------
%  g(z) = 2*G*lambda*z/(x^2+z^2)*1e5 ,   lambda = pi*R^2*DeltaRho
%  A(z) = dg/dz = 2*G*lambda*(x^2 - z^2)/(x^2+z^2)^2*1e5
maxIter = 20;  tol = 1e-6;                                % tol on |dz| [m]
[z_estimated, zHist, A, AT, ATA, ATd] = ...
        gaussNewtonDepth(G, M0, z_initial, x, gObserved_cylinder, maxIter, tol);

g_estimated_cylinder = cylinderGravity(G, M0, z_estimated, x);
Estimated_misfit     = g_Observed_cylinder - g_estimated_cylinder;
MSE_estimated        = calculateMSE(Estimated_misfit);

%figure 3
figure;
plot(x, g_Observed_cylinder,  'k-', 'LineWidth', 3); axis tight; hold on;
plot(x, g_estimated_cylinder, 'bo', 'LineWidth', 2);
grid on; box on; hold off;
title('2D Cylinder Gravity Anomaly: Observed vs Inverted');
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
fprintf('    CYLINDER GRAVITY DEPTH INVERSION RESULTS\n');
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
fprintf('Known mass per unit length    = %.3e kg/m\n', M0);
fprintf('MSE of initial model          = %.3e mGal^2\n', MSE_cylinder);
fprintf('MSE of estimated model        = %.3e mGal^2\n', MSE_estimated);
fprintf('Maximum residual magnitude    = %.3e mGal\n', max(abs(Estimated_misfit)));
fprintf('==============================================\n');

%% ---------------- NOISY DATA: 0.15 % AND 3 % -------------------------------
noiseLevels = [0.0015, 0.03];
for n = 1:numel(noiseLevels)
    pct  = 100*noiseLevels(n);
    gObs = addNoise(g_Observed_cylinder, noiseLevels(n));

    z_n     = gaussNewtonDepth(G, M0, z_initial, x, gObs, maxIter, tol);
    g_est_n = cylinderGravity(G, M0, z_n, x);
    MSE_n   = calculateMSE(gObs - g_est_n);

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
    plot(x, g_est_n,               'r-',  'LineWidth', 2);
    hold off; grid on;
    xlabel('Horizontal Distance (m)'); ylabel('Gravity Anomaly (mGal)');
    title(sprintf('Cylinder Depth Inversion with %.2f%% added noise', pct));
    legend('Observed noisy','Initial guess','Inverted','Location','best');

    fprintf('\n--- Noise level = %.2f%% ---\n', pct);
    fprintf('Estimated depth = %.2f m  (error %.2f %%)\n', z_n, 100*abs(z_n-z_true)/z_true);
    fprintf('MSE after inversion = %.3e mGal^2\n', MSE_n);
end

%% ---------------- MONTE CARLO: 200 NOISE REALISATIONS ----------------------
nMC = 200;
fprintf('\nMonte Carlo test (%d noise realisations), depth:\n', nMC);
for n = 1:numel(noiseLevels)
    est = zeros(nMC,1);
    for t = 1:nMC
        gObs   = addNoise(g_Observed_cylinder, noiseLevels(n));
        est(t) = gaussNewtonDepth(G, M0, z_initial, x, gObs, maxIter, tol);
    end
    fprintf('  %.2f%% noise: mean = %.2f m, std = %.2f m, mean |error| = %.3f %%\n', ...
            100*noiseLevels(n), mean(est), std(est), 100*mean(abs(est-z_true))/z_true);
end
end

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function M = cylinderMass(R, deltaRho)                 % Mass per unit length [kg/m]
    M = pi*R.^2 .* deltaRho;                           % pi*R^2 (NOT 4/3*pi*R^2)
end

function g = cylinderGravity(G, M, z, x)               % [mGal]
    g = (2*G*M*z) ./ (x.^2 + z.^2) * 1e5;
end

function MSE = calculateMSE(misfit)
    MSE = mean(misfit(:).^2);
end

function A = sensitivityMatrixCylinder(G, M, z, x)     % dg/dz  (depth Jacobian)
    A = 2*G*M*(x.^2 - z.^2) ./ (x.^2 + z.^2).^2 * 1e5;
end

function gN = addNoise(g, level)
    gN = g + level*max(abs(g))*randn(size(g));
end

function [z, zHist, A, AT, ATA, ATd] = gaussNewtonDepth(G, M, z0, x, gObs, maxIter, tol)
    z = z0;  zHist = z0;
    for it = 1:maxIter
        d   = gObs(:) - cylinderGravity(G, M, z, x(:));
        A   = sensitivityMatrixCylinder(G, M, z, x(:));
        AT  = A';  ATA = AT*A;  ATd = AT*d;
        dz  = inv(ATA) * ATd;
        zNew = z + dz;
        while zNew <= 0, dz = dz/2; zNew = z + dz; end % keep depth positive
        z = zNew;  zHist(end+1) = z;                   %#ok<AGROW>
        if abs(dz) < tol, break; end
    end
end
