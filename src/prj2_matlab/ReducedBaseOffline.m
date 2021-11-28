using DrWatson
@quickactivate "MIT6.339"
using DelimitedFiles
using MAT
using JLD2

% -------------------------------------------------------------------------
%
%   ReducedBasisOffline.m
%
% -------------------------------------------------------------------------
%   Reduced basis off-line: computes the reduced basis matrices and source 
%   vector from a sample.
% -------------------------------------------------------------------------


% Load grid triangulations
matdata = matread(srcdir("prj2_matlab","grids.mat"))

% Mesh level to be used
mesh = matdata["coarse"];

% Load sample matrix and get the number of samples
load sn.dat;
... 

% Solve the sample cases and construct Z
Z = zeros(mesh["nodes"],N);
....

% Initialization of reduced-base matrix and vector
ANq = cell(6,1);
FN = zeros(N);

% Calculate Anq and FN
....

% Save the reduced matrix and vector
@save "ReducedBasis.jld2" ANq FN