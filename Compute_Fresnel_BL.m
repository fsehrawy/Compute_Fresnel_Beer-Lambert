% This script implements Fresnel equations for direct incidence (angle = 0)
% Written on June 8th, 2020
% Reference calculator used for comparison:
% https://edavies.me.uk/2008/10/fresnel/

clear
close all
clc

lambda=[600]; % Wavelength [nm]
layer_names=["ARC" "ARC" "Textile" "Silicon"]; % Layers in a studied structure starting with the incidence layer
layer_th=[100, 100, 100, 500]; % Layer thickness [nm]
RI=[1.5, 2, 2.5+0.01*1i, 3+0.25*1i]; %]; % Complex refractive index

n_in=1; % Input medium (air)
n_out=1; % Output medium (air)

RI=[n_in,RI,n_out]; % Define index of input and output mediums
layer_th=[0,layer_th,0]; % Define thickness of input and output mediums

IF=struct([]); % Define interfaces data structure
layer=struct([]); % Define layer data structure

% Define coefficients of the first interface
IF(1).r=((n_in-real(RI(2)))^2)/((n_in+real(RI(2)))^2); % Reflection coeff.
IF(1).t=(4*n_in*real(RI(2)))/((n_in+real(RI(2)))^2); % Transmission coeff.
% Note: reflection and transmission coefficents are in lowercase r & t
% while the reflected and transmitted power are in uppercase R & T

% Calculate the coefficients for the rest of the interfaces between layers
for i=2:length(layer_names)+1
    IF(i).r=((real(RI(i))-real(RI(i+1)))^2)/((real(RI(i))+real(RI(i+1)))^2);
    IF(i).t=(4*real(RI(i))*real(RI(i+1)))/((real(RI(i))+real(RI(i+1)))^2);
end

% Initialize the reflected, absorbed, and transmitted power of each layer
for i=1:length(layer_names)+2
    layer(i).R=0;layer(i).A=0;layer(i).T=0;
end

% Initialize the reflected, absorbed, and transmitted power at input (from air)
layer(1).R=IF(1).r; layer(1).A=0; layer(1).T=IF(1).t;

% Iteratively computes the R/A/T power of each layer to simulate the propagation of light
for j=1:length(layer_names)+1

    % Loop through the different layers
    for i=2:length(layer_names)+1
        % Implement the Beer-Lambert model to calculate the single pass absorbance in a bulk material
        % The material is characterized by an absorption coefficient, alpha [cm^-1]        
        Abs_SP=1-exp(-4*pi*imag(RI(i))/(lambda*1e-7)*(layer_th(i)*1e-7));

        % Calculate the absorbed power in a layer which contains three components:
        % 1. Single pass of light through a material: "Abs_SP*layer(i-1).T"
        % 2. Double pass of the remaining light reflected at the interface with the next layer, IF(i+1): "(1-Abs_SP)*layer(i-1)*Abs_SP.T*IF(i).r"
        % 3. The reflected power from the next layer, R(i+1): "Abs_SP*layer(i+1).R"
        layer(i).A = layer(i).A + Abs_SP*layer(i-1).T + (1-Abs_SP)*layer(i-1).T*Abs_SP*IF(i).r + Abs_SP*layer(i+1).R;

        % Calculate the power reflected from a layer, which is the light that
        % goes through a layer twice without being absorbed
        layer(i).R=layer(i).R + (1-Abs_SP)^2*layer(i-1).T*IF(i).r; 

        % Calculate the power transmitted to the next layer
        layer(i).T= layer(i).T + (1-Abs_SP)*layer(i-1).T*IF(i).t;
        layer(i-1).T=0;
    end
end

% Ensure that the total R + A + T power = 1
disp(struct2table(layer))
disp(['R + A + T = ' num2str(sum(sum(table2array(struct2table(layer)))))]);
