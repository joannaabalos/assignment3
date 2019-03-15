%% 3 
rng(2)
clear
global C

C.q_0 = 1.60217653e-19;             % electron charge
C.m_0 = 9.10938215e-31;             % electron mass
C.mn = 0.26*C.m_0;                  % Effective Electron Mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.T = 300;                          % Kelvin

%New Assignment 3 constants
Vapp = 2; %Applied voltage on semiconductor
eConc = 10e15; %Concentration

vth = sqrt(2*C.kb*C.T/C.mn); %Thermal velocity
MTBC1 = 0.2e-12; %Mean time between colissions (s)

numPart = 30000; %Number of particles
xlim = 200e-9;
ylim = 100e-9;
nx = 200;
ny = 100;
xyarea = xlim*ylim;
bin = numPart/10;
dt = ylim/vth/100; %Scale time

%Random starting positions
x=rand(1,numPart)*xlim;
y=rand(1,numPart)*ylim;

%Defining box boundaries for indexing
xBoxBound = x>0.79e-7 & x<1.21e-7;
topBox = xBoxBound & y>0.59e-7;
bottomBox = xBoxBound & y<0.41e-7;
inBoxes = topBox | bottomBox; %Inside both box regions

%Draw boxes
figure(7)
boxX = [0.8e-7 0.8e-7 1.2e-7 1.2e-7];
boxBottom = [0 0.4e-7 0.4e-7 0];
plot(boxX,boxBottom,'color','k');
axis ([0 xlim 0 ylim])
hold on
boxTop = [1e-7 0.6e-7 0.6e-7 1e-7];
plot(boxX,boxTop,'color','k');

while(sum(inBoxes)>0)
    
    %Assign new random location
    x(inBoxes) = rand(1,length(x(inBoxes)))*xlim;
    y(inBoxes) = rand(1,length(y(inBoxes)))*ylim;
    
    %Re-checking box boundaries
    xBoxBound = x>0.8e-7 & x<1.2e-7;
    topBox = xBoxBound & y>0.6e-7;
    bottomBox = xBoxBound & y<0.4e-7;

    inBoxes = topBox | bottomBox; %Inside box box regions
    
end

%Random angle
hAngle = 360; %highest angle
lAngle = 0; %lowest angle
angle = (hAngle-lAngle).*rand(1,numPart) + lAngle; %Random angle within range

%Random MB velocity
MBfact = vth;
MBvx = randn(1,numPart)*MBfact;
MBvy = randn(1,numPart)*MBfact;
MBv = sqrt(MBvx.^2+MBvy.^2);

%Random MB velocities travelling at random angle
vx = MBvx.*cos(angle);
vy = MBvy.*sin(angle);

%Scatter probability
Pscat = 1-exp(-dt/MTBC1);

%Initialize Mean Free Path vector
MFPs = zeros(1,numPart);

%Subset of particles to be plotted
numPartPlot = 7;
subset = randi([1,numPart],numPartPlot,1).';

%New Assignment 3 calculations
Etot = Vapp/xlim; %E = V/d
Fapp = Etot*C.q_0; %F = Eq
a = Fapp/C.m_0; %F = ma 

%Increase applied voltage to see acceleration
Vapp = 1.5;
Etot = Vapp/xlim; %E = V/d
Fapp = Etot*C.q_0; %F = Eq
a = Fapp/C.m_0; %F = ma 

maxTime = 1000;
for time=1:maxTime
    %Scattering
    scatter = Pscat > rand(1,numPart); %Particles that will scatter
    vx(scatter) = 0;
    vy(scatter) = 0;    
    angle = (hAngle-lAngle).*rand(1,numPart) + lAngle; %Random angle within range
    
    %Random MB velocity
    MBvx = randn(1,numPart)*MBfact; %New x component velocity
    MBvy = randn(1,numPart)*MBfact; %New y component velocity
    MBv = sqrt(MBvx.^2+MBvy.^2); %New actual particle velcoity
   
    %Random MB velocities travelling at random angle
    vx(scatter) = MBvx(scatter).*cos(angle(scatter));
    vy(scatter) = MBvy(scatter).*sin(angle(scatter));

    %Box boundaries
    xBoxBound = x>=0.79e-7 & x<=1.21e-7;
    topBox = xBoxBound & y>=0.59e-7;
    bottomBox = xBoxBound & y<=0.41e-7;
    
    %Bouncing off of horizontal part of box
    boxHorTop = (vy>=0 & topBox); %Particles bouncing off horizontal part of top box
    boxHorBottom = (vy<=0 & bottomBox); %Particles bouncing off horizontal part of bottom box
    y(boxHorTop & y>=0.4e-7 & y<=0.6e-7) = 0.6e-7;
    y(boxHorBottom & y>=0.4e-7 & y<=0.6e-7) = 0.4e-7;  
    boxHorTopBottom = boxHorTop | boxHorBottom;
    vy(boxHorTopBottom) = -1.*vy(boxHorTopBottom); %Reverse velocity       
%     vy(boxHorBottom) = -1.*vy(boxHorBottom); %Reverse velocity    
    
    %y boundaries
    yBoundTop = y >= ylim;
    y(yBoundTop) = ylim;
    yBoundBottom = y<=0;
    y(yBoundBottom) = 0;
    yBound = yBoundTop | yBoundBottom;
    vy(yBound) = -1.*vy(yBound); %Reverse velocity 
    
    %Updating y position
    yPrev = y;
    y = y + vy*dt; 
    
    %x boundaries
    rightBound = (x>=xlim & vx>=0); %Positive xvelocities reaching right boundary
    x(rightBound) = 0; %Relocate particle to left side
    leftBound = (x<=0 & vx<=0); %Negative xvelocities reaching left boundary
    x(leftBound) = xlim; %Relocate particle to right side
    
    %Bouncing off of vertical part of boxes
    xBoxLeftBound = x>=0.79e-7;
    xBoxRightBound = x<=1.21e-7;
%     xBoxBound = xBoxLeftBound & xBoxRightBound;
%     topBox = xBoxBound & y>=0.59e-7;
%     bottomBox = xBoxBound & y<=0.41e-7;
%     bounceLeft = xBoxLeftBound & topBox &vx >=0;
%     bounceRight = xBoxLeftBound & topBox &vx <=0;
%     x(bounceLeft) = 0.8e-7;
%     x(bounceRight) = 1.2e-7;  
%     vx(bounceLeft) = -1.*vx(bounceLeft); %Reverse velocity
%     vx(bounceRight) = -1.*vx(bounceRight); %Reverse velocity
    
    boxVert = (topBox | bottomBox); %Particles bouncing off vertical part of either box
    leftVert = boxVert & vx>=0;
    rightVert = boxVert & vx<=0;
    x(leftVert ) = 0.79e-7;
    x(rightVert) = 1.21e-7;  
    vx(leftVert) = -1.*vx(leftVert); %Reverse velocity
    vx(rightVert) = -1.*vx(rightVert); %Reverse velocity
  
    %Updating x position
    xPrev = x;
    vx = vx + a*dt; %adding acceleration to velocity
    x = x + vx*dt + .5*a*dt^2;
    
    %Plotting subset
    for i=1:numPartPlot
        plot([xPrev(subset(i)) x(subset(i))],[yPrev(subset(i)) y(subset(i))])
    end
    title('Particles Trajectories with Rectangle Boundaries')
    axis ([0 xlim 0 ylim])
    drawnow
    hold on
    
    %Semiconductor temperature
    v = sqrt(vx.^2+vy.^2);
    overallTemp = C.mn*sum(v.^2)/(2*C.kb);
    avgTemp(time) = overallTemp/numPart;  
    
end

% The following figure displays the electron denisty map.

%Electron Density Map
figure(8)
EDM = hist3([y',x'],[20,20]); 
pcolor(EDM)
title('Electron Density Map (particles/cm^2)')

%Electron Temperature Map
%Define mesh and organize particles into regions
xTempLim = linspace(0,xlim,40);
yTempLim = linspace(0,ylim,40);
xTempReg = discretize(x,xTempLim);
yTempReg = discretize(y,yTempLim);

for xx=1:1:39
    for yy=1:1:39            
        %Total velocities in region
        inTempReg = (xTempReg == xx) & (yTempReg == yy);
        totInTempReg = sum(inTempReg);
        vxTot=sum(vx(inTempReg));
        vyTot=sum(vy(inTempReg));
        vTot = sqrt((vxTot)^2+(vyTot)^2);
        
        %Calculate Temperature
        tempMap(yy,xx) = C.mn*vTot^2/(2*C.kb);
    end
end

%Temperature Surface Plot
figure(9)
surf(tempMap);
title('Temperature Map')
xlabel('x Region')
ylabel('y Region')
zlabel('Temperature (K)')
colorbar;
view(0,90)

%Temperature Surface Plot
figure(10)
surface(tempMap);
title('Temperature Map (degrees K)')
colorbar;
view(45,45)

%Electron Density Map
figure(11)
EDM = hist3([y',x'],[50,50]); 
surface(EDM)
title('Electron Density Map')
xlabel('x Region')
ylabel('y Region')
zlabel('Density (particles/cm^2)')


%Electron Density Map
figure(12)
EDM = hist3([y',x'],[50,50]); 
surf(EDM)
view(45,45)
title('Electron Density Map')
xlabel('x Region')
ylabel('y Region')
zlabel('Density (particles/cm^2)')