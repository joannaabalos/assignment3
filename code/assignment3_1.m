%% Assignment 3 Monte-Carlo/Finite Difference Method

clear

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.m_0 = 9.10938215e-31;             % electron mass
C.mn = 0.26*C.m_0;                  % Effective Electron Mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.T = 300;                          % Kelvin

%New Assignment 3 constants
Vapp = .1; %Applied voltage on semiconductor
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

%Displaying final calcs
fprintf('Part 1:');
fprintf('\nThe electric field is %d V/m.',Etot);
fprintf('\nThe force on each electron is %d N.\n',Fapp);

%Increase applied voltage to see acceleration
Vapp = 1.5;
Etot = Vapp/xlim; %E = V/d
Fapp = Etot*C.q_0; %F = Eq
a = Fapp/C.m_0; %F = ma .

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
    
    %Updating x position
    xPrev = x;
    vx = vx + a*dt; %adding acceleration to velocity
    x = x + vx*dt + .5*a*dt^2;
    
    %Plotting
    figure(1)
    title('Particles Trajectories Using Maxwell-Boltzmann Distribution')
    for i=1:numPartPlot
        plot([xPrev(subset(i)) x(subset(i))],[yPrev(subset(i)) y(subset(i))])
    end
    axis ([0 xlim 0 ylim])
    drawnow
    hold on
    
    %Semiconductor temperature
    v = sqrt(vx.^2+vy.^2);
    overallTemp = C.mn*sum(v.^2)/(2*C.kb);
    avgTemp(time) = overallTemp/numPart;  
    
    %Mean Free Path
    MFPs(scatter) = 0;
    notScatter = ismissing(scatter,0); %Prep index for particles that have not scattered
    MFPs(notScatter) =  MFPs(notScatter) + v(notScatter)*dt; %Add distance travelled by particles if not scattered
    MFP = sum(MFPs)/numPart;
    
    %Mean Time Between Collisions
    MTBC = MFP*numPart/sum(v);
    
    %Calculating drift current
    mu = sum(v)/numPart/Etot;
    driftI(time) =C.q_0*eConc*mu*Etot/xyarea;
end

%Plotting Current
xtime = linspace(0,maxTime,maxTime);
figure(2)
plot(xtime,driftI)
title('Current')
xlabel('Time (s)')
ylabel('Current (A)')

%Plotting Temperature over time
figure(3)
plot(xtime,avgTemp)
title('Average Semiconductor Temperature Over Time')
xlabel('Time (s)')
ylabel('Temperature (K)')

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
figure(4)
surf(tempMap);
title('Temperature Map')
xlabel('x Region')
ylabel('y Region')
zlabel('Temperature (K)')
colorbar;

%Electron Density Map
figure(5)
EDM = hist3([y',x'],[50,50]); 
surf(EDM)
title('Electron Density Map')
xlabel('x Region')
ylabel('y Region')
zlabel('Density (particles/cm^2)')

