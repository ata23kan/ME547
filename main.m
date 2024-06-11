%% ME 547 Term Project
% Authors: Serkan Can, Onur Ata, Atakan Aygün
clear all;
clc;
close all;

% Optimizing Diani Rey Model equation for uniaxial tension 

% data_fit = 'uniaxial_fit';
data_fit = 'equibiaxial_fit';
% data_fit = 'biaxial_fit';

switch data_fit
  case 'uniaxial_fit'
    ut_fit
  case 'equibiaxial_fit'
    et_fit
  case 'biaxial_fit'
    bt_fit
end

