function [CL,CD,CA,CN,CZ,CM]=coeffgen(x,r,z,aerodynamic_flag,options,cg_offset,AOA)
% For a given geometry, mesh quality, cg offset and angle of attack (alpha)
% this function gives the lift, drag and pitching moment coefficient
%
% Inputs:
% x,r,z = body geometry --- x and r are the outermold line and z is the
% triad direction
% aerodynamic_flag = 'Newtonian', 'Modified-Newtonian', or 'Free-Molecular'
% options = variable structure with the data for different aerodynamic_flag
%           None for 'Newtonian'
%           options.gamma and options.Mach for 'Modified-Newtonian'
%           options.sigmaN, options.sigmaT, and options.TwTinf also for
%           'Free-Molecular'
% cg_offset = 3 x 1 vector showing offset from nose in x, r, and z
% AOA = angle of attack in radians
%
% Outputs:
% CL,CD,CA,CN,CZ,CM = lift, drag, axial, normal, and side force and
% pitching moment coefficients
%
% By: Soumyo Dutta and Milad Mahzari (originally March 2011)
% Updated: Novemeber 18, 2012 (changed by Soumyo Dutta)

% Determine what type of aerodynamics is being run
switch aerodynamic_flag
    case 'Newtonian'
        mod_newton_flag = 0;
        free_mol_flag = 0;
        gamma = NaN;    % Ignored, set arbitarily here
        Mach = NaN;     % Ignored, set arbitarily here 
        sigmaN = NaN;   % Ignored, set arbitarily here
        sigmaT = NaN;   % Ignored, set arbitarily here
        TwTinf = NaN;   % Ignored, set arbitarily here
    case 'Modified-Newtonian'
        mod_newton_flag = 1;
        free_mol_flag = 0;
        gamma = options.gamma;
        Mach = options.Mach;
        sigmaN = NaN;   % Ignored, set arbitarily here
        sigmaT = NaN;   % Ignored, set arbitarily here
        TwTinf = NaN;   % Ignored, set arbitarily here
    case 'Free-Molecular'
        mod_newton_flag = 0;
        free_mol_flag = 1;
        gamma = options.gamma;
        Mach = options.Mach;
        sigmaN = options.sigmaN;
        sigmaT = options.sigmaT;
        TwTinf = options.TwTinf;
    otherwise
        error('Incorrect aerodynamic flag chosen')
end

% Define velocity vector with angle of attack as the angular displacement
unitvel = -[cos(AOA);sin(AOA);0];

% Total area
refarea = pi*max(max(r))^2;         % Reference area for all aero coeffs
lref = 2*max(max(r));               % Reference length for CM calculations

% CG coordinates
cg_loc = cg_offset.*[max(max(x));max(max(r));max(max(z))];

% Preallocate for center of grid arrays
xc = zeros(size(x));
rc = xc;
zc = xc;
% Intialize force and moment coefficents
CA = 0; CN = 0; CZ = 0; CM = 0;
% Loops through all the pts. of the mesh
[rlength, philength]=size(x);
for jj = 1:1:philength-1
    for ii = 1:1:rlength-1
        % Define four (or three) corners of the mesh
        A = [x(ii,jj);r(ii,jj);z(ii,jj)];
        B = [x(ii+1,jj);r(ii+1,jj);z(ii+1,jj)];
        C = [x(ii+1,jj+1);r(ii+1,jj+1);z(ii+1,jj+1)];
        D = [x(ii,jj+1);r(ii,jj+1);z(ii,jj+1)];
        % Find unit normal, area and mesh center
        [unitnormal,area,center]=normal_func(A,B,C,D);
        xc(ii,jj) = center(1); rc(ii,jj) = center(2); zc(ii,jj)= center(3);
        % Local Incidence angle
        theta = acos(sum(unitvel.*unitnormal));
        % Calculate pressure and shear coefficients
        [Cpi,Cti] = Cp_calc(theta,mod_newton_flag,gamma,Mach,free_mol_flag,sigmaN,sigmaT,TwTinf);
        % Force coefficient vector
        shearvec = (cross_analytic(unitnormal,cross_analytic(unitvel,unitnormal)));
        CFvec = Cpi.*(-unitnormal) + Cti*shearvec./norm(shearvec);
        % Force and moment coefficients (component)
        CAi = CFvec(1); CNi = CFvec(2); CZi = CFvec(3);
        deltar = cg_loc - center; % Moment arm
        CMoment = cross_analytic(deltar,CFvec);
        CMi = CMoment(3);       % Note that pitching direction is z axis
        % Add differential force and moment coefficient to total
        % coefficient
        CA = CA + CAi*(area/refarea);
        CN = CN + CNi*(area/refarea);
        CZ = CZ + CZi*(area/refarea);
        CM = CM + CMi*(area)/(refarea*lref);
    end
end

% Calculate drag and lift coefficient from axial and normal coefficients
rotatedForces = body2stability(AOA)*[CN;CA;CZ];
CL = rotatedForces(1);
CD = rotatedForces(2);

return

function [Cpi,Cti] = Cp_calc(theta,mod_newton_flag,gamma,M,free_mol_flag,sigmaN,sigmaT,TwTinf)
% Uses Newtonian aerodynamics with shadowing (but not body shadowing)
% Also has flag free molecular
% Note theta is the angle between the unit normal and the velocity vector

if free_mol_flag == 0        % Continuum flow
    if mod_newton_flag == 0  % Classical Newtonian
        if abs(theta)<=pi/2  % Not in shadow
            Cpi = 2*cos(theta)^2;
        else                 % In shadow
            Cpi = 0;
        end
    else                % This is modified Newtonian
        Cpmax = (2/(gamma*M^2))*(((gamma+1)^2*M^2/(4*gamma*M^2-2*...
            (gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M^2)/(gamma+1))-1);
        if abs(theta)<=pi/2  % Not in shadow
            Cpi = Cpmax*cos(theta)^2;
        else                 % In shadow
            Cpi = 0;
        end
    end
    Cti = 0;
else                    % This is free molecular
    if abs(theta)<pi/2  % Not in shadow
        psi = pi/2 - theta;     % This is incidence angle
        Sinf= M*sqrt(gamma/2);
        Cpi=(1/Sinf^2)*(((2-sigmaN)*Sinf*sin(psi)/sqrt(pi)+(sigmaN/2)*...
            sqrt(TwTinf))*exp(-Sinf^2*sin(psi)^2)+((2-sigmaN)*...
            (0.5+Sinf^2*sin(psi)^2)+(sigmaN/2)*sqrt(TwTinf)*...
            sqrt(pi)*Sinf*sin(psi))*(1+erf(Sinf*sin(psi))));
        Cti=-(sigmaT*cos(psi)/(Sinf*(sqrt(pi))))*(exp(-Sinf^2*sin(psi)^2)+...
            sqrt(pi)*Sinf*sin(psi)*(1+erf(Sinf*sin(psi))));
    else
        Cpi = 0;
        Cti = 0;
    end
end

return

function vec = cross_analytic(vec1,vec2)

skewmat = [0 -vec1(3) vec1(2);vec1(3) 0 -vec1(1);-vec1(2) vec1(1) 0];
vec = skewmat*vec2;

return

function [normal_line,area,center]=normal_func(A,B,C,D)
% Using the boundary pts. of a grid, this finds the grid normal, surface
% area of the grid and the geometric grid center
% Automatically deals with situations where the grid is triangular as then
% pt. A and D will be the same

% Vectors between pts
vec1 = B-A;
vec2 = C-B;
vec3 = D-A;
vec4 = D-C;
% Vector lengths
vec4size = norm(vec4);
vec3size = norm(vec3);
vec2size = norm(vec2);
vec1size = norm(vec1);
% Different calculations if the area is triangular instead of a
% quadrilateral (note only vec2 or vec3 can be 0 length i.e. vertex of a
% triangle; vec1 and vec4 will always be non-zero length--assumed)
if vec2size == 0
    % Normal to the plane (unit vector)
    normal_line = cross_analytic(vec3,vec1);
    normal_line = normal_line./norm(normal_line);
    % Angle between the vectors that form a corner of the quadrilateral
    angle = acos(sum(vec3.*vec1)./(vec3size*vec1size));
    % Height of the quadrilateral
    h = vec1size*sin(angle);
    % Area of teh quadrilateral (and the triangle when A = D)
    area = 0.5*(vec2size+vec3size)*h;
elseif vec3size == 0
    % Normal to the plane (unit vector)
    normal_line = cross_analytic(vec2,vec1);
    normal_line = normal_line./norm(normal_line);
    % Angle between the vectors that form a corner of the quadrilateral
    angle = acos(sum(vec2.*vec1)./(vec2size*vec1size));
    % Height of the quadrilateral
    h = vec1size*sin(angle);
    % Area of teh quadrilateral (and the triangle when A = D)
    area = 0.5*(vec2size+vec3size)*h;
else
    % Normal to the plane (unit vector)
    normal_line = cross_analytic(vec2,vec1);
    normal_line = normal_line./norm(normal_line);
    % Angle between the vectors that form a corner of the quadrilateral
    angle = acos(sum(vec1.*vec2)./(vec2size*vec1size));
    % Height of the quadrilateral
    h = vec1size*sin(angle);
    % Area of teh quadrilateral (and the triangle when A = D)
    area = 0.5*(vec2size+vec3size)*h;
end

% Geometric center of the quadrilateral
center = 0.25.*(A+D+B+C);

return

function mat = body2stability(phi)
% Angle rotation in the 3-direction (z-direction)

mat = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];

return

function fun = centerfun(t,A,B,C,D)
% Defines the root-finding problem to find the center of a quadrilateral or
% triangular shape

x = t(1);
y = t(2);
z = t(3);

xm1 = (A(1)+B(1))/2;
xm2 = (C(1)+B(1))/2;
xm3 = (A(1)+D(1))/2;
xm4 = (C(1)+D(1))/2;
ym1 = (A(2)+B(2))/2;
ym2 = (C(2)+B(2))/2;
ym3 = (A(2)+D(2))/2;
ym4 = (C(2)+D(2))/2;
zm1 = (A(3)+B(3))/2;
zm2 = (C(3)+B(3))/2;
zm3 = (A(3)+D(3))/2;
zm4 = (C(3)+D(3))/2;
keyboard
fun = zeros(3,1);
if abs(xm3)<1e-9
    fun(1) = ym1 + ((ym4 - ym1)/(xm4 - xm1))*(x - xm1) + ((ym4 - ym1)/(zm4 - zm1))*(z - zm1) - ym3 - ((ym2 - ym3)/(xm2 - xm3))*(x - xm3) - ((ym2 - ym3)/(zm2 - zm3))*(z - zm3);
    fun(2) = xm1 + ((xm4 - xm1)/(ym4 - ym1))*(y - ym1) + ((xm4 - xm1)/(zm4 - zm1))*(z - zm1) - xm3 - ((xm2 - xm3)/(ym2 - ym3))*(y - ym3) - ((xm2 - xm3)/(zm2 - zm3))*(z - zm3);
    fun(3) = zm1 + ((zm4 - zm1)/(ym4 - ym1))*(y - ym1) + ((zm4 - zm1)/(xm4 - xm1))*(x - xm1) - zm3 - ((zm2 - zm3)/(ym2 - ym3))*(y - ym3) - ((zm2 - zm3)/(xm2 - xm3))*(x - xm3);
elseif abs(xm2)<1e-9
    
else
    fun(1) = ym1 + ((ym4 - ym1)/(xm4 - xm1))*(x - xm1) + ((ym4 - ym1)/(zm4 - zm1))*(z - zm1) - ym3 - ((ym2 - ym3)/(xm2 - xm3))*(x - xm3) - ((ym2 - ym3)/(zm2 - zm3))*(z - zm3);
    fun(2) = xm1 + ((xm4 - xm1)/(ym4 - ym1))*(y - ym1) + ((xm4 - xm1)/(zm4 - zm1))*(z - zm1) - xm3 - ((xm2 - xm3)/(ym2 - ym3))*(y - ym3) - ((xm2 - xm3)/(zm2 - zm3))*(z - zm3);
    fun(3) = zm1 + ((zm4 - zm1)/(ym4 - ym1))*(y - ym1) + ((zm4 - zm1)/(xm4 - xm1))*(x - xm1) - zm3 - ((zm2 - zm3)/(ym2 - ym3))*(y - ym3) - ((zm2 - zm3)/(xm2 - xm3))*(x - xm3);
end

return