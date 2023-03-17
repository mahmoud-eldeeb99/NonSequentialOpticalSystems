%%
% This simple demonstration is created to show how easy our library can be used

clear all;
clf; hold on;

sim = SequentialOpticalModel;
% Setup simulation area
sim.setBorders([0, 6, -0.15, 0.15]);
% Colorish rays! Use 'random' for rainbow. Possible templates are
% 'for-each-ray', 'natural', 'wavelength-dependent' or just specify single color
% name.
sim.ray_color_policy = 'random';
% Using template for creating rays. There are several basic templates: 'single',
% 'flat', 'spherical', 'combined', 'white'. Please, check library for more
% arguments for tweaking.
sim.createRaysFromTemplate('combined');
% This is essential to claim all previous settings were initializating process
sim.start;

% Let rays propagate 2 meters right
sim.freeSpace_new(2);
% Then go through thin lens with focal distance 2 meters
sim.thinLens_new(2);
sim.freeSpace_new(2);
% here we specify lens half height
sim.thinLens_new(2, 0.02);
% We can omit the end point, that means we will go till the right border
sim.freeSpace_new;

