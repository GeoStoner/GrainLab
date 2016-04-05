% Simple model of ripple creation in sand bend using constant flight
% distance with calculated trajectory. 
% Ryan Stoner. March 30, 2016 for modeling in the Earth Sciences

clear
clc
%% Initialize

xmin = 0;
xmax = 500;
dx = 1;

% Setting up initial distance and height, both unitless
xarr = xmin:dx:xmax;                       
yarr = ones(1,length(xarr))+(rand(1,length(xarr))-0.9); 

% Setting up impact trajectories as straight lines with slope m. We also
% need to calculate the bmax, the maximum y intercept

m = - 0.2;
bmax = yarr(length(yarr))-m*xarr(length(yarr));

% Maximum (farthest impact) and minumum impact. For bossible plotting
% purposes only
ymax = m*xarr + bmax;
ymin = m*xarr + yarr(1);           

%
figure(1)
plot(xarr,yarr,'k-')
ylim([0 5])
hold on
plot(xarr,ymax,'b-')
ylim([0 10])
figure(1)
plot(xarr,ymin,'g-')
ylim([0 5])
hold off

%% Looping
% n - number of grains
n = 10000;

% Outside of loop define width of box and average number of grains moved to
% obtain height. Also average distance of transport.
avegrloss = 8;
grwspl = 10;
aved = 20;

% Add nplot - Number of plots, tplot - grains between plot
nplot = 1000;
tplot = n/nplot;
nframe = 0;

for i= 1:n
    % First calculate bhit - random intercept to determine line of random
    % impact. yhit - equation of that line
    bhit = (bmax-yarr(1))*rand(1) + yarr(1);
    yhit = m*xarr + bhit;
    
    % Find location of impact
    indbangrange = find(yhit-yarr<0);
    indbang = min(indbangrange);
    % Find number of grains lost, randomized, rounded to make whole no.
    % Divide by width of grains in box/width of splash for a loss of a
    % certain height.
    loss = abs(avegrloss*rand(1));
    lossh = loss/grwspl; 
    
    % Loss of material in impact
    yarr(indbang) = yarr(indbang) - lossh;
    
    
    % Add grains splashed at a certain distance, wrap around boundary
    if indbang< length(xarr)-aved
    yarr(indbang + aved)= yarr(indbang + aved)+lossh;
    
    else
    yarr(indbang-length(xarr)+aved+1) = yarr(indbang-length(xarr)+aved+1)...
        + lossh;
    
    
    end
    if(mod(n,tplot)==0)
nframe = nframe+1;

    figure(2)
plot(xarr,yarr,'k-')
ylim([0 10])

pause(0.01)
    end
% Converting time to a string to print string out
% trounded = round(n(i),0);
% strtime = num2str(trounded);

   
    
end
figure(2)
plot(xarr,yarr,'k-')
ylim([0 10])
