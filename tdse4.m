function tdse4(I,run_index)
%now with variable intensity!
%  Program to solve the Schrodinger equation using sparse matrix
%  Crank-Nicolson scheme (Particle-in-a-box version)

%adding a potential and rationalising the grids etc - wsb
% now putting in real units and a time-varying field
%have to set up Hamiltonian at each time step take that inside the loop...

%try to implement using a sparse matrix formalism
% now sorted out the field matrix -


%clear all;  
%check for arguments
if nargin <1
    I = 3e14; %setting the intensity in W/cm^2 if not set as an argument.
    run_index = 1; % just a marker for incrementing filenames
end

format compact
h_bar = 1.05e-34;  
mass = 9.11e-31; 
e =1.6e-19; eps0 = 8.85e-12; 
c=3e8; % SI units

%define freq
lambda = 800e-9;            %pump wavelength
w = 2 * pi * 3e8 / lambda;
delta = 7e-15;              %pulse full width half maximum in seconds
phi = pi/2;                 %the CEP


%% setting the timescale
% two important timescales
% dt is the spacing for calculations of the wavefunction psi
% Ntotal is the number of wavefunction calculations
% Ncalc is the number of acceleration calcs, which determines highest harmonic I
% can look at in the spectrum
% Nplot is the number of recorded psi, so that I can plot the image (not relevant to calcs) 
% (best if the total calc times are 2^n)
Ncalc = 2^11;               %number of acceleration calcs; 2048 pts usually for 3fs, 8192 for 10fs
time_window = 4 * delta;    %calculation window
t_zero = 2* delta;          %set the zero time point in the middle of the window
calc_iter = 2^3;            %number of propagation calculations per acceleration calc
plot_iter = 2^2 * calc_iter;%number of propagation calculations per plot points
Ntotal = Ncalc * calc_iter; %Ntotal is bigger than Ncalc
Nplot = Ntotal/ plot_iter;  %Nplot is smaller than Ncalc, or the images are too big
dt = time_window / (Ntotal-1);
t = 0:dt:time_window;


%show the numbers
fprintf('dt step is %g attoseconds, number of Calc steps is %g\n', dt*1e18, Ncalc);
%% setting up the pulse
env = exp( -2 * log(2) *  ( (t-t_zero)/delta) .^2);  % this is the pulse FIELD envelope
carrier = cos(w * (t-t_zero) + phi); % set the CEP properly!!
%now make a pulse
E0 = sqrt( (2*I*1e4)/ (c * eps0) );
E = E0 .* carrier .* env; %the field
%calculate cutoff 
Up= e*E0^2 / (4*mass*w^2);    %in eV
cutoff_eV=15.7596+3.17*Up; %for reference  
%From the time separation of the acceleration points, work out what the max freq you can sample is...
max_energy = 2 * pi * h_bar/(2*dt*calc_iter * e);  %in eV
fprintf('Up = %g eV, cutoff is %g eV, maximum calculated energy is %g eV\n', Up, cutoff_eV, max_energy);

%% setting the length scale:
%number of length points is set using L, the length size, and h, the length
%step. Both are in terms of the Bohr radius right now.

a0 = 5.23e-11;  %Bohr rad - atomic unit length scale
L = 300 * a0; %DONT CHANGE!!   % System extends from -L/2 to L/2
h = 0.1*a0;%DONT CHANGE!! 
% load up the ground state appropriate to this L/h combo
load('300x0x1.mat'); %initialise wavefn appropriately - no check
psi = gs;

x = -L:h:L;   % Coordinates of grid points
Nx = length(x); % number of spatial grid points


%%  Set up the Hamiltonian operator matrix
coeff = -h_bar^2/(2*mass*h^2);
%this next bit sets up the static Coulomb potential as a single vector -
%the leading diagonal of the matrix
Vstatic =   -(e^2 / (4 * pi * eps0))  *  (1./sqrt(a0^2 + x.^2))' ;  %static potential - Soft Coulomb

H1 = ones(Nx,1) * coeff;
H2 = -2 * H1;
H3 = H1;

%% Set up the matrix Q

psi(1) = 0;  psi(Nx)=0;   %needed?
%set up the gobbler funciton to tidy up at boundaries;
gobbler = 1-((sin(0.5 * pi * x/L)).^18); % absorbing boundaries
%abs_bounds=1;  %turn on gobbling

%allocate some memory
psivst = zeros(Nplot, Nx);
accel = zeros(1,Ncalc);
ta = zeros(1,Ncalc);
tplot = zeros(1,Nplot);

psivst(1,:) = psi;      % Record initial condition

plotting = 0; %this decides whether some extra plots are shown

if plotting==1
    figure(2); clf;
end

%% start the loop
a_index = 0;
plot_index=0;

hh = waitbar(plot_index/Nplot, 'Iterations', 'CreateCancelBtn',...
              'setappdata(gcbf,''cancelling'',1)');  % fancy waitbar
setappdata(hh, 'cancelling', 0);
tic;
for ii=1:length(t)
    if getappdata(hh, 'cancelling')
        break
    end

%set up time-dependent Hamiltonian
   Vlaser = e * E(ii) * x' ; %column vector length x, varies with the varying laser field
   lead_diag = 0.5* (ones(Nx, 1) + (0.5*1i*dt/h_bar)*(H2 + Vstatic + Vlaser));
   off_diag = 0.5*(0.5*1i*dt/h_bar)*H1;
  
  Q = [off_diag lead_diag off_diag];  %make a 3 column version of the matrix
  chi = tri_ge(Q,psi);                %use the custom solver rather than internal function
  %chi = Q\psi; %this would use the internal solver, but needs a full (or sparse) matrix
  psi = chi - psi;     
  %if abs_bounds ==1
  psi = gobbler' .* psi;            %damp down reflections from boundaries
  %end
 

% Periodically record values for calculations
% acceleration calculation every calc_iter points
  if( rem(ii,calc_iter) < 1 )  % calculate the acceleration if appropriate
    a_index = a_index+1; %index for calcualtion
    %need the full Hamiltonian at this time
    lead_diag = H2 + Vstatic + Vlaser;
    full_ham = spdiags([H1 lead_diag H3],-1:1,Nx,Nx); %use sparse matrix for accel calc
    psistar = conj(psi);
    ta(a_index) = t(ii);           % Record current time
    % record points for plotting every plot_iter points
      if( rem(ii,plot_iter) < 1 )  % Every plot_iter steps plotrecord 
          plot_index=plot_index+1;  %index for p_plot
          waitbar(plot_index/Nplot, hh);
          %fprintf('Iteration %d of %d\n', plot_index, Nplot);
          tplot(plot_index) = t(ii);  %record times of plotting points
          psivst(plot_index, :) = psi;  % and P(x,t) for plots
      end
    %electron_fraction(plot_index) = sum(psivst(plot_index, :)) * h; % how much electron have we left?

%now calc some expectation values
    term1 = psistar .* (full_ham * (full_ham * (x'.*psi))) ;
    term2 = psistar .* (full_ham * (x' .* (full_ham * psi)));
    term3 = psistar .* x' .* (full_ham * (full_ham * psi)) ;
    accel(a_index) = (h/h_bar^2) * ( sum(term1) - 2* sum(term2) + sum(term3) );
  end
end
delete(hh); %remove waitbar
elapsed_time = toc;
fprintf('Elapsed time was %g seconds\n', elapsed_time);
psivst = psivst';  %switcheroo to conform to easier viewing


%% Plot probability versus position at various times

figure(3); clf;
imagesc(tplot*1e15, x*1e9, log10(psivst.*conj(psivst))); shading flat; colorbar;%colormap jet%
%imagesc(tplot*1e15, x*1e9, (psivst.*conj(psivst))); shading flat; colorbar;%colormap jet%
ylim([-8 8])
%caxis([4 10]);
caxis([-11 -1]);
ylabel('x /nm');  xlabel('t /fs'); zlabel('P(x,t)');
colormap jet
hold on
fudge_factor = 1e-10;
plot(t*1e15, E * fudge_factor , 'k' );
hold off

%% check the acceleration
figure(5);
plot(ta*1e15, real(accel));
xlabel('time /fs');
ylabel('acceleration /ms-2');
title('The acceleration, or dipole moment');

%% announce how long the calculation takes per acceleration calc cycle
fprintf('Time per cycle = %g ms\n', 1000 * elapsed_time / Ncalc);

%% Frequency-domain calculations
accel_sample = accel(1:end);  % ability to look at smaller pieces - not used.
Ew = fft(accel_sample); %take the FT of the acceleration - this is the spectrum
acc_length = length(accel_sample);
f = (1/ (dt * calc_iter) * (1:acc_length)/acc_length);  %calculate frequency scale

odd_harm = (2 * pi * h_bar * c  / (e*lambda)) * (1:2:151); %calculate where the odd harmonics will be
figure(11);
plot(2*pi*h_bar*f/e, log10(Ew.*conj(Ew)),'.-');
hold on;
plot(odd_harm,47,'.r');
hold off;


title('The Harmonic Spectrum');
ylabel('log10(harmonic intensity)')
xlabel('Energy /eV');
xlim([50 100]);
ylim([36 56]);

%% ring that bell
beep; pause(0.2);beep; pause(0.5); beep
%make a datestamped filename
%foo = clock;
%fname = sprintf('run%d%d', hour(now), minute(now));
%fname = sprintf('run%d', run_index);
% now save the data into a sensible file for processing.
%save(fname, 'E', 'accel', 'dt','delta', 't','t_zero', 'x','L','ta', 'tplot', ...
   % 'E0', 'lambda', 'I', 'calc_iter', 'plot_iter','phi');



