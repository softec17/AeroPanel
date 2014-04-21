% Runs coeffgen to generate hypersonic aerodynamic database for a given
% geometry
clear all; close all; clc;

%% User input for geometry
shape_choice = 'Generate';          % 'Load' or 'Generate'
shape_def = 'Sphere';             % Optional; only need if 'Generate'
                                    % 'Sphere-Cone', 'Biconic', 'Torus',
                                    % 'Hyperboloid', 'Ogive', 'Sharp-Cone',
                                    % 'Sphere', 'Cylinder', 'Disk'
shape_geom_file = 'shape_file_generic.txt';

% Aerodynamics
aerodynamic_flag = 'Free-Molecular';
options.gamma = 1.4;
options.sigmaN = 1;     % Diffuse reflection (interaction of air w/ engineering surfaces)
options.sigmaT = 1;     % Diffuse reflection (interaction of air w/ engineering surfaces)
options.TwTinf = 1;     % Cold surface
cg_offset = zeros(3,1);

%% Define rotation matrix about x axis
R_1 = @(phi) [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

%% Load geometry
switch shape_choice
    case 'Load'
        
        % Loads axisymmetric 2D shape
        geom_data = importdata(shape_geom_file);
        Xlist = geom_data.data(:,1);
        Rlist = geom_data.data(:,2);
        
        % Copies of the profile pts (x,r,z) coordinates
        rp = Rlist;
        xp = Xlist;
        zp = zeros(size(xp));
        
        % Generate axisymetric 3D shape by rotating 2D shape
        num_meridians = 360;
        phi = linspace(0,2*pi,num_meridians);
        % Preallocate space
        r = zeros(length(rp),length(phi));
        x = r; z = r;
        % Initialize
        r(:,1) = rp; x(:,1) = xp; z(:,1) = zp;
        % Revolves profile 360 degrees
        for jj = 1:1:length(phi)
            mat = R_1(phi(jj));
            for ii = 1:1:length(rp)
                holder = mat*[xp(ii);rp(ii);zp(ii)];
                x(ii,jj) = holder(1); r(ii,jj) = holder(2); z(ii,jj) = holder(3);
            end
        end
        
    case 'Generate'
        % Loads axisymmetric 2D shape
        geom_data = importdata(shape_geom_file);
        
        theta = deg2rad(geom_data(1));
        rn = geom_data(3);
        rd = geom_data(2);
        num_r = geom_data(4);
        num_phi = geom_data(5);
        flare_angle = deg2rad(geom_data(6));
        flare_D = 2*geom_data(7);
        
        % Predefine all flags to zero
        flag = 0; hyp_flag = 0; torus_flag = 0; ogive_flag = 0; sharp_cone_flag = 0; sphere_flag = 0; cylinder_flag = 0; disk_flag = 0;
        switch shape_def
            case 'Sphere-Cone'
                flag = 0;
            case 'Biconic'
                flag = 1;
            case 'Hyperboloid'
                hyp_flag = 1;
            case 'Torus'
                torus_flag = 1;
            case 'Ogive'
                ogive_flag = 1;
            case 'Sharp-Cone'
                sharp_cone_flag = 1;
            case 'Sphere'
                sphere_flag = 1;
            case 'Cylinder'
                cylinder_flag = 1;
            case 'Disk'
                disk_flag = 1;
        end
        [x,r,z]=creategeom(theta,rn,rd,num_r,num_phi,flag,flare_angle,flare_D,hyp_flag,torus_flag,ogive_flag,sharp_cone_flag,sphere_flag,cylinder_flag,disk_flag);
end

%% Set aerodynamic parameters and options
MList = [50,40,30];
MachColor = {'b','k','r'};

CLList = zeros(length(MList),1);
CDList = CLList;
CAList = CLList;
CNList = CLList;
CZList = CLList;
CMList = CLList;
for ii = 1:length(MList)
    options.Mach = MList(ii);
    %% Generate aerodynamics using coeffgen.m
    AOAList = 0; % Angle of Attack
    AOAListrad = deg2rad(AOAList);
    for jj = 1:length(AOAList)
        AOA = AOAListrad(jj);
        [CL,CD,CA,CN,CZ,CM]=coeffgen(x,r,z,aerodynamic_flag,options,cg_offset,AOA);
    end
    CLList(ii) = CL;
    CDList(ii) = CD;
    CAList(ii) = CA;
    CNList(ii) = CN;
    CZList(ii) = CZ;
    CMList(ii) = CM;
end

% Display CA's
disp('Mach List')
disp(MList')
disp('CA List')
disp(CAList)

% Display CD's (Will be different for parallel flow on a disk). -Kyle EDIT
disp('CD List (Parallel Disk)')
disp(CDList)

% Plots mesh
% figure(7)
% col = zeros(size(x));
% h = mesh(x,r,z,col);
% axis equal