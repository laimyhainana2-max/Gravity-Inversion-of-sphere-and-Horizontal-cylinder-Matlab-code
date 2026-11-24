[README Gravity Inversion Using MATLAB code.md](https://github.com/user-attachments/files/35380147/README.Gravity.Inversion.Using.MATLAB.code.md)
**Gravity Inversion Using MATLAB code for buried sphere and horizontal cylinder**



Real gravity data contain noise from instrument drift, terrain, regional trends and data processing, so a method must be tested on data whose true answer is known before it is trusted on field data.  Synthetic profiles were therefore generated from known models and contaminated with Gaussian noise.  Four MATLAB scripts were written, one for each combination of body and unknown.  Sphere\_Model\_Code and Cylinder\_Model\_Code2 invert for the density contrast of a sphere and of a horizontal cylinder, with the depth and radius known; this is a linear problem with one parameter and is solved in a single step.  Sphere\_Model\_Code3\_depth2 and Cylinder\_model\_code\_depth invert for the depth of a sphere and of a cylinder, with the density contrast and radius known; because gravity is not linear in depth, the correction is repeated until the depth stops changing.  The sphere is a three-dimensional body and the infinite horizontal cylinder a two-dimensional body, so the two geometries represent the three- and two-dimensional cases of the objectives; all anomalies are calculated along one profile.



**Horizontal Cylinder Model and Sphere Code**

1\. Generating synthetic data (Forward Model) – For density contrast and depth respectively

Each script is a MATLAB function with no inputs or outputs, so a single file holds the whole experiment.  The commands ‘clc’ and ‘close all’ clear the command window and close open figure windows so that earlier runs do not interfere with the current run.  The constants and true model are then defined: the gravitational constant, a density contrast of 1000 kg/m³, a depth (z) of 200 m, a radius (R) of 50 m and 100 evenly spaced observation points from −500 m to 500 m, with the body centred at x = 0.  These values define the noise-free target data set that a gravimeter would record at the surface.  The excess mass (M0) is calculated with the local function sphereMass, equation (6), or cylinderMass, equation (8) (mass per unit length, πR²Δρ), and the observed anomaly is calculated with sphereGravity, equation (5), or cylinderGravity, equation (7).  The factor 10⁵ converts m/s² to mGal.



2\. Initial guess and error calculation

The initial guess is deliberately wrong.  In the density scripts the density contrast starts at 900 kg/m³ and the depth is kept at its true value; in the depth scripts the depth starts at 180 m, with the density contrast and radius (and therefore the mass M0) known.  The forward model is run with the initial guess to give the calculated anomaly (for the density scripts the initial mass M1 is calculated with sphereMass or cylinderMass), and the misfit is calculated with equation (10).  The mean-squared error (equation 9, function calculateMSE) summarises how far the initial guess is from the data.



3\. Linear least squares inversion

The sensitivity (Jacobian) matrix A holds the partial derivatives of the anomaly with respect to the unknown at each of the 100 observation points: equations (11b) and (11c) for the density contrast and equations (11) and (11a) for the depth.  The data residual vector (d) holds the difference between the observed and calculated data.  The matrix A is 100 × 1, so AT is 1 × 100, ATA is 1 × 1 and ATd is 1 × 1, and the over-determined equation (2) reduces to dividing one number by another.  For the density contrast the problem is linear, so one step gives the exact correction.  Depth enters the forward model non-linearly, so A is recalculated for the current depth and equation (2) is solved repeatedly (Gauss-Newton iteration) until the depth correction is smaller than 10⁻⁶ m or 20 iterations have been done.  If a correction would make the depth negative it is halved until the depth is positive, and the depth after every iteration is stored in the array zHist.  In Matlab, the backslash operator is usually used for efficiency but for this research paper the ‘inv’ operator is used so that the algebra is identical to equation (2); it first calculates the inverse matrix and then multiplies it by the data vector, which is slower than the one-step backslash operator, although for a 1 × 1 matrix the difference is negligible.



4\. Model update and verification

The correction is added to the current model to give the estimated density contrast or depth.  The forward model is run with the estimated parameter to give the inverted anomaly, and the new misfit and MSE are calculated.  The initial and final MSE, the maximum residual and the number of iterations are printed to show whether the error has decreased towards zero (noise-free data) or towards the noise variance (noisy data).



5\. Noise addition

This part of the code evaluates how well the code performs in a noisy, realistic environment.  Gaussian noise generated with ‘randn’ is scaled to 0.15 % and 3 % of the peak absolute anomaly and added to the noise-free data.  The random number generator seed is locked using rng(1) so that the noise is the same each time the code is run, which makes debugging and checking easier.  In the density scripts the noisy data are inverted with the same matrices A and ATA, because A does not depend on the data; in the depth scripts the whole iteration is repeated.  Each noise level is run once, so the noisy results are for a single noise realisation.



6\. Visual output

Figure 1 of each script plots the noise-free observed anomaly against the calculated anomaly (gravity anomaly in mGal on the y-axis and distance in m on the x-axis); figure 2 compares the observed anomaly with the initial-guess model; figures 3 and 4 compare the observed anomaly with the inverted model; figures 5 and 7 plot the noisy data (0.15 % and 3 %) against the noise-free anomaly; and figures 6 and 8 compare the noisy data with the initial guess and the inverted model.  Observed data are plotted as symbols and model curves as lines so that the quality of fit can be judged visually as well as numerically.  See the attached appendices 5 – 8 for the visual outputs.

&#x20;

