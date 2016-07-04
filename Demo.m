% demo code of "Direct 3D Pose Estimation of a Planar Target"
%
% Usage:
%   Run demo.m
%
% Disclaimer:
%   It is provided for educational/researrch purpose only.
%   Please cite the paper if you find the code useful.
%
%   Direct 3D Pose Estimation of a Planar Target
%   Hung-Yu Tseng, Po-Chen Wu, Ming-Hsuan Yang and Shao-Yi Chien
%   IEEE Winter Conference on Applications of Computer Vision, WACV 2016
%
% Contact:
%   Hung-Yu Tseng
%   hytseng@media.ee.ntu.edu.tw
clc; clear all; close all;
Marker = im2double(imread('imgs/timage.png')); % target image
tDim = 0.12; % length of the shorter size of the target. Here is 12 cm.
I = im2double(imread('imgs/cimage.jpg')); % camera image
[height, width, ~] = size(I);

f = [500.858378, 501.2506075]; % camera focal length
p = [320.645466, 179.1686375]; % camera principle point
in_mat = [f(1),0,p(1),0;0, f(2),p(2),0;0,0,1,0;0,0,0,1]; % camera intrinsic matrix
exmat = Test_DPE(Marker, I, in_mat, tDim/2, 0.2, 0.7, 0.25, 1, 1, 1); % main program
%exmat = Test_APE(Marker, I, in_mat, tDim/2, 0.2, 0.7, 0.25, 1, 1, 1); % main program

% render x, y, z axis on the camera image
f = figure('Position', [150 150 width height]);
[corner_x, corner_y] = draw_coordinate(exmat, in_mat, tDim/2);
imagesc(I); 
truesize;
axis off; hold on;
plot([corner_x(1);corner_x(2)], [corner_y(1);corner_y(2)], 'r', 'LineWidth', 3);
plot([corner_x(1);corner_x(3)], [corner_y(1);corner_y(3)], 'g', 'LineWidth', 3);
plot([corner_x(1);corner_x(4)], [corner_y(1);corner_y(4)], 'b', 'LineWidth', 3);