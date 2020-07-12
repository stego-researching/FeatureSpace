% This compile.m is used under Windows
clc;clear;
fprintf('compiling genDistMapC ... ');

%% 64 bit Windows and Matlab
mex -output gen_distmap_single gen_distmap.cpp

fprintf('done\n');
% copyfile('gen_distortion_mapC.mexw64','../');
% fprintf('Compiled BIM MEX file was copied to ../../../ folder.\n');
