%%
% This simple demonstration is created to show how to specify
% wavelength-dependent refractive index

clear all;
clf; hold on;

sim = SequentialOpticalModel;
sim.setBorders([0, 6, -0.15, 0.15]);
% It is convinient here to use 'natural' or 'wavelength-dependent' options
% because they were created to represent visible wavelength
sim.ray_color_policy = 'natural';
sim.createRaysFromTemplate('rainbow', 50);
sim.start;

sim.freeSpace_new(2);
% Here we specify chromatic abberation. We use Couchy formula, where refractive
% index n = A + B/lambda^2 + C/lambda^4 + ... Second argument is zero, which
% means automatically ajusted height of lens. Third argument is array of
% coefficients for Couchy formula. Coefficients below are for fused silica glass
% (from Wikipedia)
sim.thinLens_new(4, 0, [1.458, 3.54 * 10^(-15)]);

% Now we should see, that rays aren't focused at the same point. But it hard to
% notice it. Let uncomment "reverse" operation with diverging lens with the same material.
%sim.freeSpace_new(2);
%sim.thinLens_new(-2, 0, [1.458, 3.54 * 10^(-15)]);
% Finally, let see what is happening at the long distance.
%sim.setBorders([0, 60, -0.15, 0.15]);

sim.freeSpace_new;

