function [x,r,z]=creategeom(theta,rn,rd,num_r,num_phi,flag,flare_angle,flare_D,hyp_flag,torus_flag,ogive_flag,sharp_cone_flag,sphere_flag,cylinder_flag,disk_flag)
% Generates mesh for a given geometry
% Theta will be in radians

% Torus
if (torus_flag == 1)
    %         [x,r,z]=
    %         creategeom(theta,rn,rmax,num_r,num_phi,flag,flare_angle,flare_D,hyp_flag,torus_flag);
    r = rn;
    R = rd;
    Xlist = [linspace(-r,r,num_r/2),linspace(-r,r,num_r/2)];
    Rlist = [sqrt(r^2-Xlist(1:length(Xlist)/2).^2)+R,-sqrt(r^2-Xlist(length(Xlist)/2+1:end).^2)+R];

% Ogive
elseif (ogive_flag == 1)
    D = rd;
    fineness = rn;      % Fineness is l/d
    xmax = fineness*D;
    Xlist = linspace(0,xmax,num_r);
    Rlist = (sqrt((fineness.^2+0.25).^2-(Xlist./D-fineness).^2)-(fineness.^2-0.25)).*D;
    
% Sharp cone
elseif (sharp_cone_flag==1)
    % Prellocating space and initializing radius direction
    Rlist = linspace(0,rd,num_r);
    TANT = tan(theta);
    % Calculates x-location of pts. along the contour of the axisymmetric body
    Xlist = Rlist./TANT;
    
% Sphere
elseif (sphere_flag==1)
    % Prellocating space and initializing x direction
    D = 2*rd;
    Xlist = linspace(0,D,num_r);
    Rlist = sqrt(rd.^2 - (Xlist - rd).^2);
    
% Circular Cylinder
elseif (cylinder_flag==1)
    % Prellocating space and initializing x direction
    D = 2*rd;
    xmax = 2*rn + D;
    Xlist = linspace(0,xmax,num_r); Xlist = Xlist(:);
    num_x = length(Xlist);
    Rlist = zeros(num_x,1);
    for ii = 1:1:num_x
        if Xlist(ii) <=rn
            Rlist(ii) = sqrt(rn.^2 - (rn - Xlist(ii)).^2);
        elseif Xlist(ii)<=(rn+D)
            Rlist(ii) = Rlist(ii-1);
        else
            Rlist(ii) = sqrt(rn.^2 - (Xlist(ii) - (rn+D)).^2);
        end
    end
    
% Disk
elseif (disk_flag==1)
    % Prellocating space and initializing x direction
    xmax = rn;
    Xlist = linspace(0,xmax,num_r); Xlist = Xlist(:);
    Rlist = zeros(size(Xlist));
    Rlist(1) = 0; Rlist(end) = 0;
    Rlist(2:end-1) = rd;
         
else
    
    % Sphere-cone
    if (hyp_flag == 0)
        % Prellocating space and initializing radius direction
        Rlist = linspace(0,rd,num_r);
        Xlist = zeros(size(Rlist));
        SINT = sin(theta);
        COST = cos(theta);
        TANT = tan(theta);
        len_r = length(Rlist);  % Number of pts. in the r direction
        % Limits in the r and x direction for the sphere
        rlim = rn*COST;
        xlim = rn*(1-SINT);
        % Calculates x-location of pts. along the contour of the axisymmetric body
        for ii = 2:len_r
            if Rlist(ii)>rlim       % If the pt. is on the cone
                Xlist(ii) = xlim + (Rlist(ii)-rlim)/TANT;
            else                    % If the pt. is on the sphere
                Xlist(ii) = rn - sqrt(rn.^2-Rlist(ii).^2);
            end
        end
        
        rmax_new = flare_D/2;
        % Bi-conic (flared)
        if (flag==1)
            Rlist_new=zeros(1,10);
            Xlist_new=zeros(1,10);
            delta_r=(rmax_new-Rlist(end))/10;
            for i=1:10
                Rlist_new(i)=delta_r*i;
                Xlist_new(i)=Rlist_new(i)/tan(flare_angle);
            end
            
            Rlist_new=Rlist_new+Rlist(end);
            Rlist=[Rlist Rlist_new];
            Xlist_new=Xlist_new+Xlist(end);
            Xlist=[Xlist Xlist_new];
        end
        
    % Hyperboloid 
    else
        Xlist=linspace(0,rd,num_r);
        Rlist = sqrt(2.*rn.*Xlist + Xlist.^2.*tan(theta)^2);
    end
end

% Copies of the profile pts (x,r,z) coordinates
rp = Rlist;
xp = Xlist;
zp = zeros(size(xp));

% Phi = clock angle (about x axis) for the axisymmetric body (from 0 to
% 2pi)
phi = linspace(0,2*pi,num_phi);
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

return

function mat = R_1(phi)
% Angle rotation in the 1-direction (x-direction)

mat = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

return