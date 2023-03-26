%%
% Diffraction on single slit. Should result in sinc^2 in Fraunhofer limit.

clf; hold on;

% Slit width
a = 25 * 10^(-6);

% Distance between sources and image
L = 1;
% Simulation height
H = 0.1;

% We want to have enough ray to approximate integral, but we do not want to draw
% all of them. Let define how many rays we are needed (lower this number for
% speed up, increase for more smooth and precise picture)
N = 10^2;

sim = SequentialOpticalModel;
sim.setBorders([0, 1.2 * L, -H/2, H/2]);
N_initial = 60;
sim.createRaysFromTemplate('flat', N_initial);

% It is very consuming to draw all rays. Let draw only 50 of them.
sim.n_visible_rays = 50;
% Let place rays next to slits. And first 50 to be more-less representable
sim.rays(1, 1 : sim.n_visible_rays) = linspace(-a/2, a/2, sim.n_visible_rays);
sim.rays(1, sim.n_visible_rays + 1 : N_initial) = linspace(-a/2, a/2, N_initial - sim.n_visible_rays);

% Also we want to switch off autoscale, because we are interested only in
%crossection
sim.autoscale_enabled = 0;
sim.start;

% Let collimated beam propagate a little bitand
sim.freeSpace_new(0.2 * L);

% Now we want to create double slit
sim.obstacle(-H, - a/2);
sim.obstacle(a/2, H);

% Let add some rays, that would not be showed, but will be
%calculated as normal rays. We will use all these rays to compute intensity
sim.createHiddenRays(N);

% This commented out line is lens. Lens with high positive f effectively moves
% away the sources (increases L_eff) and "magnify" slit width.
%
% sim.freeSpace_new(L / 2);
% sim.thinLens_new(L * 2 / 3);

sim.freeSpace_new;
% Let calculate intensity and plot it! First argument is resolution, feel free to change

sim.calcIntensity;

