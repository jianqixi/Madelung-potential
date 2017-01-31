%% To calculate madelung constant from ewald summation 
%Based on Comput. Mater. Sci. 123, 131-138 (2016)
%If use this code, please cite the above Reference. Thanks very much!

n =20;
kn = 20;
sigma = 0.8; % parameter for gaussian

V = 0;           %total potential
V_sr = 0;        %short-range potential
V_lr = 0;        %long-range potential
V_self = 0;      %self-interaction potential
V_neutral = 0;   %neutral background potential
S = 0;           %structure factor

r_tmp = [];
k_tmp = [];
r_list = []; % direct lattice
k_list = []; % reciprocal lattice 

epsilon = 1*[1 0 0;0 1 0;0 0 1]; 
% dielectric constant from experimental result, 
% if dielectric constant is isotorpic epsilon=1, otherwise, epsilon=[]

a=15.48;       %the length of unit cell in hexagonal, the unit is Angstrom
c=20.224;      %the length of unit cell in hexagonal, the unit is Angstrom
Vc=sqrt(3)/2*a*a*c;
%Vc=a.^3;
%fprintf('the volume of supercell is %f\n',Vc);
li=[1 0 0];
lj=[0 1 0];
lk=[0 0 1];
%%%%%%%%%%%%%%%%%%%for the hexagonal cell%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%the direct vectors
ai=a*li;
aj=-sin(pi/6)*a*li+cos(pi/6)*a*lj;
ak=c*lk;
%%%%the reciprocal vectors
ki=2*pi/(cos(pi/6)*a)*(cos(pi/6)*li+sin(pi/6)*lj);
kj=2*pi/(cos(pi/6)*a)*lj;
kk=2*pi/c*lk;
%%%%%%%%%%%%%%for the simple cubic cell%%%%%%%%%%%%%%
%ai=a*li;
%aj=a*lj;
%ak=a*lk;

%ki=2*pi/a*li;
%kj=2*pi/a*lj;
%kk=2*pi/a*lk;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This part is the description of geometry
% direct lattice
    for ni = -n:n
        for nj = -n:n
            for nk = -n:n               
                r_list = ni*ai+nj*aj+nk*ak;
                r = r_list;  % input all the element of the 
                             % ir row in matrix r_list into vector r sequently 
                %fprintf('the value of r %f %f %f\n',r);
                rii= sqrt(r*inv(epsilon)*r');    
                %fprintf('the value of rii %f\n',rii);
                if(rii ~= 0)
                    V_sr = V_sr + 1/sqrt(det(epsilon)) *1/rii * erfc(rii * sigma);   
                end 
                %fprintf('the value of short-range part:%f \n',V_sr);
            end
        end
    end
    fprintf('the value of short-range part:%f \n',V_sr);

% reciprocal lattice
    for bi = -kn:kn
        for bj = -kn:kn
            for bk = -kn:kn
                k_list = bi*ki+bj*kj+bk*kk; 
                k=k_list;
                %fprintf('the value of k %f %f %f\n',k);
                kii=k*epsilon*k';
                %fprintf('the value of kii %f\n',kii);
                if(kii~=0)
                    V_lr = V_lr + real(4 * pi / Vc  * exp(-kii/(4*sigma*sigma))/kii);
                end
               %fprintf('the value of long-range part:%f \n',V_lr);
            end
        end
    end
fprintf('the value of long-range part:%f \n',V_lr);

% self interaction potential
V_self = -2 * sigma/ (sqrt(pi*det(epsilon)));
fprintf('the value of self part:%f \n',V_self);

% neutral background potential
V_neutral = -pi/(Vc*sigma*sigma);
fprintf('the value of neutral part:%f \n',V_neutral);

% total electrostatic potential
V = V_sr + V_lr + V_self+V_neutral;
fprintf('the value of total part:%f \n',V);


