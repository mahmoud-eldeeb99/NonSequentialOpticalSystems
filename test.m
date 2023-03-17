a= SequentialOpticalModel;





figure(1);clf;hold on;
[rays0,rayColors] = a.createRays(2);
xlabel('Optic Axis')
ylabel('Image Plane')
ylabel('Object Plane')
%% Ray Tracing
LensNum=2;
d=[.5 .3];
L=[2 4 6];
Lm=cumsum(L);
n=[2 3];
p1=[5 5];
p2=[6 6];
semiDia=[1 3];

rays1 = a.freeSpace(rays0,L(1));
a.drawRays(rays0,rays1,0,L(1),rayColors);

lens(rays1,p1,p2,d,n,L,Lm,rayColors,LensNum,semiDia)


function out= lens(rays1,p1,p2,d,n,L,Lm,rayColors,LensNum,semiDia)

for i= 1:LensNum
a= SequentialOpticalModel;
rays2 = a.thickLens(rays1,p1(i),p2(i),d(i),n(i));
% Draw and label Lens 1
line([Lm(i) Lm(i)],[semiDia(i)*min(rays1(1,:)) semiDia(i)*max(rays1(1,:))],'LineStyle','-','Color',[0 0 0],'LineWidth',40*d(i)) 
% text(1.05*L1(i),1.05*min(rays1(1,:)),'L_1')
% Propogate over free space distance L2 (to location L1+L2)
rays3 = a.freeSpace(rays2,L(i+1));
a.drawRays(rays2,rays3,Lm(i),Lm(i)+L(i+1),rayColors);
rays1=rays3;

end

    
    



end


