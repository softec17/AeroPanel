%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample user input deck for AeroPanel code
%
% Necessary inputs:
% choice            : 'Load' or 'Generate'
% aerodynamic_flag  : 'Free-Molecular' 'Newtonian' 'Modified-Newtonian'
% cg_offset         : 3 x 1 array with offset from the 0,0,0 pt. of the grid. Expects quantities in non-dimensionalized form with max(x), max(y) and max(r) for normalization 
% MachList          : L x 1 -- L size list of Mach number
% AOAList           : N x 1 -- N size list of total angle of attack (deg)
% shape_def         : (Optional -- needed for 'Generate')
%                       'Sphere-Cone', 'Biconic', 'Torus',
%                       'Hyperboloid', 'Ogive', 'Sharp-Cone',
%                       'Sphere', 'Cylinder', 'Disk'
% options.gamma     : (Optional -- need for 'Free-Molecular' and 'Modified-Newtonian) specific heat ratio (air = 1.4)
% options.sigmaN    : (Optional -- need for 'Free-Molecular') diffuse reflection (interaction of air w/ engineering surfaces)
% options.sigmaT    : (Optional -- need for 'Free-Molecular') diffuse reflection (interaction of air w/ engineering surfaces)
% options.TwTinf    : (Optional -- need for 'Free-Molecular') ratio of wall to freestream temperature (1 = cold surface)
%
% Outputs:
% CLList, CDList, CAList, CNList, CZList, CMList = L x N list of data
% CL = lift force coefficient, CD = drag force coefficient, CN = normal force coefficient
% CZ = side force coefficient, CA = axial force coefficient, CM = pitching moment coefficient
%
% Assumptions:
% Let r = points in the radial direction
% refarea = pi*max(max(r))^2;         % Reference area for all aero coeffs
% lref = 2*max(max(r));               % Reference length for CM calculations
% So for a sphere-cone, refarea = max base area and lref = max diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice = 'Generate';         
shape_def = 'Sphere';                                         
shape_geom_file = 'genericshapefile.txt';
aerodynamic_flag = 'Free-Molecular';
options.gamma = 1.4;
options.sigmaN = 1;     
options.sigmaT = 1;     
options.TwTinf = 1;    
cg_offset = zeros(3,1);
AOAList = [0:1:10]';
MachList = [40,50]';

% Run AeroPanel -- do not change
[CLList, CDList, CAList, CNList, CZList, CMList] = ...
    aeropanel(aerodynamic_flag,choice,shape_def,shape_geom_file,options,cg_offset,AOAList,MachList);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the shape_geom_file format
% --------------------------------------------
% Updated: April 20, 2014 (Soumyo Dutta)
%
% Two options:
%
% 1. Input geometry (use with 'Load' option) and then evaluate. 
%    Assumes pts. of the 1/4 profile.  Col. 1 = x pts. Col. 2 = r pts.
%
% 2. Generate geometry and then evaluate shape
%    Assumes entry in this order
%    Row 1: Half-cone angle, deg
%    Row 2: Radius, m
%    Row 3: Nose radius, m
%    Row 4: Num_r -- number of points in mesh in radial direction (grid description)
%    Row 5: Num_phi -- number of meridians generated (angular direction spacing) (grid description)
%    Row 6: Flare angle, deg., Optional; needed for biconic to define the second cone half angle (ignored for all other shapes)
%    Row 7: Flare radius, m, Optional; needed for biconic to define the second cone's max. radius (ignored for all other shapes)
% See below for more discussion on the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Half-cone angle, deg
% -----------------------
% Sphere-Cone: Cone half angle
% Biconic's first half cone angle
% Sharp-Cone: Cone half angle
% Ignored by all other shapes
%
%
% Radius, m
% -----------------------
% Sphere-Cone: Max. Radius
% Biconic: First max. radius
% Sharp-Cone: Max. Radius
% Ogive: r = D/2 for ogive profile
% Sphere: r = D/2 --- radius of the sphere
% Circular Cylinder: length = 2*r --- length of the cylinder
% Torus: large radius
% Hyperboloid: Length = r --- length of the hyperbola
% Disk: r = radius of the circular section of disk
%
%
% Nose Radius, m
% -----------------------
% Sphere-Cone: Nose Radius
% Biconic: Nose radius
% Sharp-Cone: Ignored
% Ogive: fineness ratio
% Sphere: Ignored
% Circular Cylinder: radius of the spherical caps
% Torus: little radius
% Hyperboloid: radius/nose radius used to generate profile
% Disk: thickness of the disk