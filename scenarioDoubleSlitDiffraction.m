%%
% This demonstration is here to show how simple things may imply complex
% phenomena (all physics is about it). Here we are showing, that we able to
% simulate some diffraction phenomenon: double-slit diffraction. We can do it due
% to three things: considering ray phase and changes of phase over the
% propagating, creating new sub rays according to the Huygens–Fresnel principle,
% energy and momentum conservation laws, and calculating intensities at the end
% with respect to phase.
%
% Since for someone it may seem to be controversal, here the sketch of proof:
% there is plane wave decomposition in Fresnel–Kirchhoff diffraction formula, we
% can approximate this integral as summ and consider rays as plane waves. We
% are also using Huygens–Fresnel principle which is consistent with
% Fresnel–Kirchhoff diffraction formulsim. There are no principal restrictions.
% Although, sometimes we still have to use outstanding number of rays, due to
% radial nature of Huygens–Fresnel principle.

clf; hold on;

% Distance between two sources. Default wavelength is 650nm, we like handy
% numbers
D = 650 * 10^(-6);
% Slit width
a = 80 * 10^(-6);

% Distance between sources and image. Using formula for two-slit diffraction
% (d = L * lambda / D) we know what are we looking for: 1mm between maximas of
% intensity (if slit width -> 0, otherwise it was slightly less due to sinc^2
L = 10;
% Simulation height
H = 0.1;

% We want to have enough ray to approximate integral, but we do not want to draw
% all of them. Let define how many rays we are needed (lower this number for
% speed up, increase for more smooth and precise picture)
N = 10;

sim = SequentialOpticalModel;
sim.setBorders([0, 1.2 * L, -H/2, H/2]);
N_initial = max(51, fix(20 * (D + a) / (2 * a)));
sim.createRaysFromTemplate('flat', N_initial);

% It is very consuming to draw all rays. Let draw only 50 of them.
sim.n_visible_rays = 50;
% Let place rays next to slits. And first 50 to be more-less representable
sim.rays(1, 1 : 50) = linspace(-D/2 - a/2, D/2 + a/2, 50);
sim.rays(1, 51 : N_initial) = linspace(-D/2 - a/2, D/2 + a/2, N_initial - 50);

% Also we want to switch off autoscale, because we are interested only in
%crossection
sim.autoscale_enabled = 0;
sim.start;

% Let collimated beam propagate a little bitand
sim.freeSpace_new(0.2 * L);

% Now we want to create double slit
sim.obstacle(-H, -D/2 - a/2);
sim.obstacle(-D/2 + a/2, D/2 - a/2);
sim.obstacle(D/2 + a/2, H);

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

