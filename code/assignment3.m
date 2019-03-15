%% ELEC 4700 Assignment-3 Monte-Carlo/Finite Difference Method
% Assignment 3 - Joanna Abalos 100962263

close all
clear
clc

% In this assignment, 30 000 particles are modelled to calculate
% temperatures, make models and observations using a combination of Monte-
% Carlo modeling and Finite Difference Method. 7 particles are plotted to 
% observe their trajectories.

assignment3_1

%Increasing the applied voltage to 1.5V shows the particle's curved
%trajectories as a result of the static electric field. 

%The relationship between electron drift current density and average
%carrier velocity is the following: 
%driftcurrent =q*eConc*mu*E/area where mu = velocity/#particles/E
%Over time, the current over the entire semiconductor becomes/approaches a 
%constant value.

assignment3_2
assignment3_3

%The particles are most dense at one side of the bottle-neck opening. This
%is due to the applied voltage and resulting electric field forcing the
%electrons right-ward, and bouncing off the boxes. The next step to make
%the simulation more accurate is to make the electrons react to each other.
%This would include bouncing off and repelling each other.

