clear; clc;
%%%%%%%%%%%%%%%%%%%% Material and loading properties %%%%%%%%%%%%%%%%%%%%
t_p_values = [1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-09, 1e-08, 1e-07]; % Sampling pulse durations in seconds unit
Fluence = [1e01, 1e02, 1e03, 1e04, 1e05, 1e06, 1e07, 1e08]; % Sampling laser fluence in J.m^{-3} unit
l = 2500e-9; % Sample length in m unit
t = 10000e-12; % Simulation time in seconds unit
sample_x = 100; % Total number of sampling points w.r.t distance
sample_t = 100; % Total number of sampling points w.r.t time
rho = 8890; % Density in kg/m^{3} unit
kappa = 91; % Lattice thermal conductivity in W.m^{-1}.K^{-1} unit
C_e = 3.6886*10^5; % Electronic specific heat capacity in J.m^{-3}.K^{-1} unit
C_l = 400*rho; % Lattice specific heat capacity in J.m^{-3}.K^{-1} unit
beta = 1/(13.5*1e-9); % Optical absorption coefficient in m^{-1} unit
G = 7.6e+10; % Shear Modulus in Pa unit
B = 1.8e+11; % Bulk Modulus in Pa unit
beta_T = 1.3*10^-5; % Volumetric thermal expansion coefficient in K^{-1} unit
G_el = 3.6*10^17; % Electron-phonon coupling coefficient in W.m^{-3}.K^{-1} unit
T0 = 298; % Reference temperature in K unit
v = -B*beta_T/(B+4*G/3); % Temperature coefficient of wave speed in K^{-1} unit
ce = sqrt((B+4*G/3)/rho); % Bulk sound velocity in m.s^{-1} unit
beta2 = beta^2;
gamma = beta/kappa;
%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%
x = linspace(0,l,sample_x);
t = linspace(0,t, sample_t);
max_sigma_field = zeros(length(t_p_values), length(Fluence)); % Maximum stress w.r.t pulse duration and fluence
min_sigma_field = zeros(length(t_p_values), length(Fluence)); % Minimum stress w.r.t pulse duration and fluence
for time = 1:length(t_p_values)
    t_p = t_p_values(time);
    [a_i, b_i,k,w0] = compute_coefficients(t_p); % Fourier coefficient calculation
    i = linspace(1,10000,10000)'; iw0 = (i.*w0);
    iw02 = iw0.^2; iw02ce2 = iw02./ce^2;
    Int = Fluence./t_p; % Laser intensity in W.m^{-3} unit
    a_i = Int.*a_i'; b_i = Int.*b_i';
    theta_field = zeros(length(x),length(t),length(Fluence)); % Lattice temperature w.r.t length, pulse duration and fluence
    theta_e_field = theta_field; % Electronic temperature w.r.t length, pulse duration and fluence
    u_field = theta_field; % Displacement w.r.t length, pulse duration and fluence
    sigma_field = theta_field; % Stress w.r.t length, pulse duration and fluence
    %%%%%%%%%%%%%%%% Fourier parameters %%%%%%%%%%%%%%%%
    omega_1 = (1i.*iw0.*(C_l + C_e)./kappa) - (C_l.*C_e.*iw02./(kappa.*G_el));
    omega_2 = 1 + 1i.*iw0.*C_l./G_el;
    omega_3 = (B.*beta_T.*T0./G_el) .* (C_e.*iw02 - 1i.*iw0.*G_el)./kappa;
    omega_4 = (B.*beta_T.*T0.*1i.*iw0./G_el);
    Ap = (beta2 + iw02ce2).*gamma./((omega_1 - beta2.*omega_2).*(beta2...
        + iw02ce2)+v.*beta.*(omega_3.*beta + omega_4.*beta.^3));
    Bp = gamma./((omega_1 - beta2.*omega_2).*(beta2 + iw02ce2)./(v...
        .*beta) + omega_3.*beta + omega_4.*beta.^3);
    k1 = -((omega_1 + (iw02ce2.^2.*omega_2.^2 + 2.*iw02ce2.*omega_1.*omega_2 - 4.*omega_4.*iw02ce2.*omega_1.*v - 2.*iw02ce2.*omega_2.*omega_3.*v + omega_1.^2 + 2.*omega_1.*omega_3.*v + omega_3.^2.*v.^2).^(1/2) - iw02ce2.*omega_2 + omega_3.*v)./(2.*(omega_2 - omega_4.*v))).^(1/2);
    k2 = -((omega_1 - (iw02ce2.^2.*omega_2.^2 + 2.*iw02ce2.*omega_1.*omega_2 - 4.*omega_4.*iw02ce2.*omega_1.*v - 2.*iw02ce2.*omega_2.*omega_3.*v + omega_1.^2 + 2.*omega_1.*omega_3.*v + omega_3.^2.*v.^2).^(1/2) - iw02ce2.*omega_2 + omega_3.*v)./(2.*(omega_2 - omega_4.*v))).^(1/2);
    A1 = ((k1.^2 + iw02ce2).*(Ap.*k2.^3.*v - Bp.*beta.*k2.^3 - Bp.*beta.*iw02ce2.*k2 ...
        + Ap.*beta.*iw02ce2.*v + Ap.*iw02ce2.*k2.*v))./(iw02ce2.*v.*(k1 - k2).*(k1.^2 + k1.*k2 + k2.^2 + iw02ce2));
    A2 = -((k2.^2 + iw02ce2).*(Ap.*k1.^3.*v - Bp.*beta.*k1.^3 - Bp.*beta.*iw02ce2.*k1 ...
        + Ap.*beta.*iw02ce2.*v + Ap.*iw02ce2.*k1.*v))./(iw02ce2.*v.*(k1 - k2).*(k1.^2 + k1.*k2 + k2.^2 + iw02ce2));
    B1 = (-v.*k1./(iw02ce2+k1.^2)).*A1;
    B2 = (-v.*k2./(iw02ce2+k2.^2)).*A2;
    Ae1 = omega_2 .* A1 + k1.*omega_4.*B1;
    Ae2 = omega_2 .* A2 + k2.*omega_4.*B2;
    Aep = omega_2 .* Ap - beta.*omega_4.*Bp;
    %%%%%%%%%%%%%%%% Final solution %%%%%%%%%%%%%%%%
    for i = 1:length(x)
        for j = 1:length(Fluence)
            theta = A1.*exp(k1.*x(i)) + A2.*exp(k2.*x(i)) + Ap.*exp(-beta.*x(i));
            theta_e = Ae1.*exp(k1.*x(i)) + Ae2.*exp(k2.*x(i)) + Aep.*exp(-beta.*x(i));
            u = B1.*exp(k1.*x(i)) + B2.*exp(k2.*x(i)) + Bp.*exp(-beta.*x(i));
            sigma = (B + 4.*G/3).*(v.*(Ap.*exp(-beta.*x(i)) + A1.*exp(k1.*x(i)) + A2.*exp(k2.*x(i))) ...
                - Bp.*beta.*exp(-beta.*x(i)) + B1.*k1.*exp(k1.*x(i)) + B2.*k2.*exp(k2.*x(i)));
            theta_field(i,:,j) = sum(abs(theta).*(a_i(:,j).*cos(iw0.*t + angle(theta)) + ...
                b_i(:,j).*sin(iw0.*t + angle(theta)))); % Lattice temperature field 
            theta_e_field(i,:,j) = sum(abs(theta_e).*(a_i(:,j).*cos(iw0.*t + angle(theta_e)) + ...
                b_i(:,j).*sin(iw0.*t + angle(theta_e)))); % Electronic temperature field 
            u_field(i,:,j) = sum(abs(u).*(a_i(:,j).*cos(iw0.*t + angle(u)) + ...
                b_i(:,j).*sin(iw0.*t + angle(u)))); % Displacement field 
            sigma_field(i,:,j) = sum(abs(sigma).*(a_i(:,j).*cos(iw0.*t + angle(sigma)) + ...
                b_i(:,j).*sin(iw0.*t + angle(sigma)))); % Stress field 
            max_sigma_field(time, j) = max(max(sigma_field(:,:,j))); % Extracting maximum temperature w.r.t length and time for each laser fluence and pulse duration
            min_sigma_field(time, j) = min(min(sigma_field(:,:,j))); % Extracting minimum temperature w.r.t length and time for each laser fluence and pulse duration
        end
    end
end
%%%%%%%%%%%%%%%%%%%% Saving data %%%%%%%%%%%%%%%%%%%%
fname = 'max_min_stress_wrt_fluence_and_pulse_duration.mat';
save(fname, 'max_sigma_field', 'min_sigma_field', 't_p_values', 'Fluence');