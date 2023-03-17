%%
% This demonstration is here to show how simple things may imply complex
% phenomena (all physics is about it). Here we are showing, that we able to
% simulate some difraction phenomenon: double-slit difraction (difractiion of two spherical sources). We can do it due
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

clear all;
clf; hold on;

% Distance between two sources. Default wavelength is 650nm, we like handy
% numbers
D = 0.0065;
% Distance between sources and image. Using formula for two-slit difraction
% (d = L * lambda / D) we know what are we looking for: 1mm between maximas of
% intensity
L = 10;
% We want to have enough ray to approximate integral, bu we do not want to draw
% all of them. Let define how many rays we are needed (lower this number for
% speed up, increase for more smooth and precise picture)
N = 10^4;

sim = SequentialOpticalModel;
sim.setBorders([0, L, -D, D]);
sim.createRaysFromTemplate('combined', 2, 15, 1.5*D/L);
% We want 2 point source 15 radial rays each. They radiate into angle, ajusted
% for simulation area. But we need some extra tweaking for y coordinate
sim.rays(1, 1:15) = -D/2;
sim.rays(1, 16:30) = D/2;
% Also we want to switch off autoscale, because we are interested only in
%crossection
sim.autoscale_enabled = 0;

sim.start;

% Now we want to add some rays, that would not be showed, but will be
%calculated as normal rays. We will use all these rays to compute intensity
sim.createHiddenRays(N);

% This commented out line is lens. Lens with high positive f effectively moves
% away the sources (increases L_eff)
%
% sim.freeSpace_new(L / 2);
% sim.thinLens_new(L * 2 / 3);

sim.freeSpace_new;
% Let calculate intensity and plot it! First argument is resolution. We want
% higher than default value, because defaults are to be workaround at any case,
% while we have good usage scenario
sim.calcIntensity(120);

