clear all; close all; clc;

%% Colors and parameters
blue = [114 147 203]./255;
red = [211 94 96]./255;
black = [128 133 133]./255;
green = [132 186 91]./255;
brown = [171 104 87]./255;
d = 3e-3; % Thickness of 3 mm
h = 6.626e-34; % Planck's constant [J*s].
k_B = 1.381e-23; % Boltzmann constant [J/K].
c = 3e8;
T = 300;

%% Get filename
warning('off','all');
[filename,folder] = uigetfile('*.csv','MultiSelect','on');
if iscell(filename)==0
    if filename==0
        disp('No .csv file selected');
        return;
    end
    num = 1;
    filename = {filename};
else
    num = numel(filename);
end
% Preallocate the cells
wavelength_cell_refl = cell(1,num);
wavelength_cell_trans = cell(1,num);
R_cell = cell(1,num);
T_cell = cell(1,num);
Rf_cell = cell(1,num);
k_cell = cell(1,num);
n_cell = cell(1,num);
sample_name = cell(1,num);
num_set = 0;
num_single = 0;
x_by_l = linspace(0,1,1000);
J_by_J_0_cell = cell(1,num);
B_cell = cell(1,num);
A_cell = cell(1,num);
B_SI_cell = cell(1,num);
A_SI_cell = cell(1,num);
emissivity_cell = cell(1,num);
emissivity_SI_cell = cell(1,num);

%% Plotting reflectance
for i = 1:num
    name_cell = strsplit(filename{i},'.csv');
    sample_name{i} = name_cell{1};
    name_split = strsplit(name_cell{1},'_');
    if strcmp(name_split{end},'refl')
        if any(strcmp(filename,strcat(erase(sample_name{i},'_refl'),'_trans.csv')))
            num_set = num_set+1;
            T_refl = readtable([folder,filename{i}]);
            T_trans = readtable([folder,strcat(erase(sample_name{i},'_refl'),'_trans.csv')]);
            wavelength_cell_refl{i} = T_refl{:,1}./1e9;
            wavelength_cell_trans{i} = T_trans{:,1}./1e9;
            R_cell{i} = T_refl{:,2}./1e2;
            T_cell{i} = T_trans{:,2}./1e2;
            if ~isempty(wavelength_cell_refl{i}) && ~isempty(wavelength_cell_trans{i})
                if wavelength_cell_refl{i}==wavelength_cell_trans{i}
                    % Plot the reflectance and transmittance
                    A_cell{i} = 1-R_cell{i}-T_cell{i};
                    fig1 = figure(3*num_set-2+num_single);
                    subplot(1,2,1);
                    plot(wavelength_cell_refl{i}./1e-9,R_cell{i}./1e-2,'LineWidth',2,'Color',blue);
                    xlim('tight');
                    xlabel('Wavelength [nm]');
                    ylabel('Transmittance [%]');
                    grid on;
                    title('Transmittance vs Wavenumber')
                    subplot(1,2,2);
                    plot(wavelength_cell_trans{i}./1e-9,T_cell{i}./1e-2,'LineWidth',2,'Color',red);
                    xlim('tight');
                    xlabel('Wavelength [nm]');
                    ylabel('Transmittance [%]');
                    grid on;
                    title('Transmittance vs Wavelength');
                    sgtitle({'Reflectance and Transmittance Measured Using Integrating Sphere',erase(sample_name{i},'_refl')},'Interpreter','None');
                    pos = get(fig1,'Position');
                    set(fig1,'Position',pos.*[1/2,1,2,1]);
                    % Plot the complex RI and semi-infinite reflectance
                    Rf_cell{i} = (2+T_cell{i}.^2-(1-R_cell{i}).^2-((2+T_cell{i}.^2-(1-R_cell{i}).^2).^2-4.*R_cell{i}.*(2-R_cell{i})).^(1/2)).*(2.*(2-R_cell{i})).^(-1);
                    k_cell{i} = wavelength_cell_refl{i}./(4.*pi.*d).*log(Rf_cell{i}.*T_cell{i}./(R_cell{i}-Rf_cell{i}));
                    n_cell{i} = (1+Rf_cell{i})./(1-Rf_cell{i})+(4.*Rf_cell{i}./(1-Rf_cell{i}).^2-k_cell{i}.^2).^(1/2);
                    A_SI_cell{i} = 1-Rf_cell{i};
                    fig2 = figure(3*num_set-1+num_single);
                    subplot(1,2,1);
                    hold on;
                    plot(wavelength_cell_refl{i}./1e-9,n_cell{i},'LineWidth',2,'Color',red);
                    plot(wavelength_cell_refl{i}./1e-9,k_cell{i},'LineWidth',2,'Color',blue);
                    xlim('tight');
                    xlabel('Wavelength [nm]');
                    ylabel('n and k');
                    legend('Refractive index n','Extinction coefficient k','Location','Best');
                    grid on;
                    title('Refractive Index and Extinction Coefficient');
                    hold off;
                    subplot(1,2,2);
                    plot(wavelength_cell_refl{i}./1e-9,Rf_cell{i}./1e-2,'LineWidth',2,'Color',black);
                    xlim('tight');
                    xlabel('Wavelength [nm]');
                    ylabel('Reflectance [%]');
                    grid on;
                    title('Semi-Infinite Reflectance vs Wavelength');
                    sgtitle({'Complex Refractive Index and Semi-Infinite Reflectance Extracted',erase(sample_name{i},'_refl')},'Interpreter','None');
                    pos = get(fig2,'Position');
                    set(fig2,'Position',pos.*[1/2,1,2,1]);
                    % Plot the visulization of irradiance penetration
                    J_by_J_0_cell{i} = zeros(numel(wavelength_cell_refl{i}),numel(x_by_l));
                    for j = 1:numel(x_by_l)
                        J_by_J_0_cell{i}(:,j) = exp(log(T_cell{i}./(1-R_cell{i})).*x_by_l(j));
                    end
                    fig3 = figure(3*num_set+num_single);
                    imagesc(wavelength_cell_refl{i}./1e-9,x_by_l,real(transpose(J_by_J_0_cell{i})),'CDataMapping','scaled','Interpolation', 'bilinear');
                    %set(gca,'YDir','normal');
                    colormap('gray');
                    c1 = colorbar;
                    title({'Visualization of the Penetration of Irradiation',erase(sample_name{i},'_refl')},'Interpreter','None');
                    xlabel('Wavelength [nm]');
                    ylabel('Nomalized depth x/d');
                    ylabel(c1,'Normalized Intensity J/J_0');
                else
                    disp('Incompatible wavelength range');
                end
            end
        else
            num_single = num_single+1;
            T_refl = readtable([folder,filename{i}]);
            wavelength_cell_refl{i} = T_refl{:,1}./1e9;
            R_cell{i} = T_refl{:,2}./1e2;
            if ~isempty(wavelength_cell_refl{i})
                fig1 = figure(3*num_set+num_single);
                plot(wavelength_cell_refl{i}./1e-9,R_cell{i}./1e-2,'LineWidth',2,'Color',red);
                xlim('tight');
                xlabel('Wavelength [nm]');
                ylabel('Reflectance [%]');
                grid on;
                title({'Reflectance and Transmittance Measured Using Integrating Sphere',erase(sample_name{i},'_refl')},'Interpreter','None');
                pos = get(fig1,'Position');
                set(fig1,'Position',pos.*[1,1,1,1]);
            end
        end
    end
end

% for i = 1:num
%     name_cell = strsplit(filename{i},'.csv');
%     sample_name{i} = name_cell{1};
%     T = readtable([folder,filename{i}]);
%     wavenumber_cell{i} = T{:,1};
%     wavelength_cell{i} = 0.01./wavenumber_cell{i};
%     R_cell{i} = T{:,2};
%     if ~isempty(wavenumber_cell{i})
%         figure;
%         subplot(2,1,1);
%         plot(wavenumber_cell{i},R_cell{i},'LineWidth',2,'Color',blue);
%         xlim('tight');
%         xlabel('Wavenumber [cm^-^1]');
%         ylabel('Reflectance [%]');
%         grid on;
%         title('Reflectance vs Wavenumber')
%         subplot(2,1,2);
%         plot(wavelength_cell{i}./1e-6,R_cell{i},'LineWidth',2,'Color',red);
%         xlim('tight');
%         xlabel('Wavelength [μm]');
%         ylabel('Reflectance [%]');
%         grid on;
%         title('Reflectance vs Wavelength');
%         sgtitle({'IR Reflectance Measured Using Integrating Sphere',sample_name{i}},'Interpreter','None');
%     end
% end