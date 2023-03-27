


clear all;


a= SequentialOpticalModel; %% creating object of the class model

L1 = 2; % Distance from object to Lens 1
f1 = 2; % Focal length of Lens 1
L2 = 3; % Distance between Lens 1 and Lens 2
f2 = 2; % Focal length of Lens 2
L3 = 2; % Distance after Lens 2
figure(1);clf;hold on;
[rays0,rayColors] = a.createRays(4);
xlabel('Optic Axis')

ylabel('Image Plane')

ylabel('Object Plane')
%% Ray Tracing
% Propogate over free space distance L1
rays1 = a.freeSpace(rays0,L1);
a.drawRays(rays0,rays1,0,L1,rayColors);
% Transport through Lens 1
rays2 = a.thinLens(rays1,f1);
% Draw and label Lens 1
line([L1 L1],[1.1*min(rays1(1,:)) 1.1*max(rays1(1,:))],'LineStyle',':','Color',[0 0 0],'LineWidth',4) 
text(1.05*L1,1.05*min(rays1(1,:)),'f_1')
% Propogate over free space distance L2 (to location L1+L2)
rays3 = a.freeSpace(rays2,L2);
a.drawRays(rays2,rays3,L1,L1+L2,rayColors);
% % Transport through Lens 2
% rays4 = a.thinLens(rays3,f2);
% % Draw and label Lens2
% line([L1+L2 L1+L2],[1.1*min(rays3(1,:)) 1.1*max(rays3(1,:))],'LineStyle',':','Color',[0 0 0],'LineWidth',4) 
% text(1.025*(L1+L2),1.05*min(rays3(1,:)),'f_2')
% % Propogate over free space distance L3 (to location L1+L2+L3)
% rays5 = a.freeSpace(rays4,L3);
% a.drawRays(rays4,rays5,L1+L2,L1+L2+L3,rayColors);






