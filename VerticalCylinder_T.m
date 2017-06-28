function [Txx,Tyy,Txy,Txz,Tyz,Tzz,Tzz1,Tzz2] = VerticalCylinder_T(obsxyz,ModelPara)

% Calculate the gravity and gravity gradient anomalies of a vertical cylinder
%
% In calculation, the parameters needed are the radius of the cylinder 'a',
% thickness 'd', residual density 'rho', and the depth to the top 'h', 
% the lateral coordinate of the center (x0,y0), coordinate of the observation
% point is (x,y,z)
%
% Input---
%         obsxyz: coordinates of observation points
%           A cell array, obsxyz{1} and obsxyz{2} are the gridding points, 
%           a matrix, while obsxyz{3} is the calculation plane, usually 
%           equals 0, the default value is used if [] is inputted.

%         ModelPara: parameters of the cylinder disk, including
%           a,d,h,rho,x0, y0. It is an array contains four elements, each 
%           element must be assigned, i.e. the length must be 0 or 6
% Output---
%         Gravity gradient components of the cylinder and Tzz of the two
%         prisms, i.e. Tzz1, Tzz2
%
% Example---
%         x=0:99;y=0:99;
%         [obsxyz{1},obsxyz{2}] = meshgrid(x,y); obsxyz{3} = 0;
%         ModelPara = [10,30,20,1e3,50,50];    % Units are meter and kg/m^3
%         [Txz,Tyz,Tzz] = VerticalCylinder_T(obsxyz,ModelPara);
%      or [Txz,Tyz,Tzz] = VerticalCylinder_T([],[]);  
%
%
% Author: Dailei Zhang
% Date  : 2016/12/23
% Email : zhangdailei@hotmail.com
%-------------------------------------------------------------------------%


%---------------------Check the inputs----------------------------%

numObs = length(obsxyz);
if numObs == 0
    [x,y] = meshgrid(0:99,0:99);
    z = 0;
elseif numObs == 2
    x = obsxyz{1}; y = obsxyz{2};
    z = 0;
elseif numObs == 3
    x = obsxyz{1}; y = obsxyz{2}; z = obsxyz{3};
end

numPara = length(ModelPara);
if numPara == 0
     ModelPara = [10,30,20,1e3,50,50];         % Default value
end

a = ModelPara(1); d = ModelPara(2);  
h = ModelPara(3); rho = ModelPara(4);          % Residual density
x0 = ModelPara(5); y0 = ModelPara(6);

%-------------------------------Calculation----------------------%
G = 6.67e-11;               % m^3/(kg*s^2)

r = sqrt((x-x0).^2+(y-y0).^2);
A = sqrt(4*a*r./((h+z)^2+(a+r).^2));
B = sqrt(4*a*r./((h+d+z)^2+(a+r).^2));

% Complete elliptic integrals of the first and the second kind
[KA,EA] = ellipke(A);
[KB,EB] = ellipke(B);


%temp = 1./r.^2;
%[m,n] = find(temp==Inf);
%temp(m,n) = (temp(m-1,n)+temp(m+1,n)+temp(m,n-1)+temp(m,n+1))/4;

%Txz = 1e9*G*rho*(x-x0).*(sqrt((h+z)^2+(a+r).^2).*((2-A).*EA-2*KA)-...
%    sqrt((h+d+z)^2+(a+r).^2).*((2-B).*EB-2*KB));
%Txz = Txz.*temp;

%Tyz = 1e9*G*rho*(y-y0).*(sqrt((h+z)^2+(a+r).^2).*((2-A).*EA-2*KA)-...
%    sqrt((h+d+z)^2+(a+r).^2).*((2-B).*EB-2*KB));
%Tyz = Tyz.*temp;

Tzz = 2*1e9*G*rho*(h+z)*a*(-((r.^2+(h+z)^2-a^2).*EA-((r-a).^2+(h+z)^2).*KA)./...
    (((r-a).^2+(h+z)^2).*((r+a).^2+(h+z)^2))+((r.^2+(h+d+z)^2-a^2).*EB-...
    ((r-a).^2+(h+d+z)^2).*KB)./(((r-a).^2+(h+d+z)^2).*((r+a).^2+(h+d+z)^2)));

% find the zero value point in r and the corresponding position in Tzz,
% replace it with nan, then fill nan with function fillnan

[pr,pc] = find(r==0);
Tzz(pr,pc) = nan;
Tzz = fillnan(Tzz);

%Calculate Tzzz using equation (6) in the above paper
%for i=length(x)-1
%    for j=length(y)-1
%        temp(i,j) = Txz(i,j)-Txz(i+1,j)+Tyz(i,j)-Tyz(i,j+1);
%    end
%end

%Tzzz = temp;

%-------------------Components transforming-----------------------%
%

[nr,nc] = size(Tzz);         
mn1 = 2^(nextpow2(nr)+1);    
mn2 = 2^(nextpow2(nc)+1);   
rdiff = floor((mn1-nr)/2);   
cdiff = floor((mn2-nc)/2);   

data = extend(Tzz,[mn1 mn2]);      % edge extension
[row, col] = size(data);

fz = fft2(data);                   % 2D-FFT
fz = fftshift(fz);           

dx = x(1,2)-x(1,1); dy = y(2,1)-y(1,1);

wnx = 2*pi/(dx*(row-1));           % wavenumber
wny = 2*pi/(dy*(col-1));
cx = row/2+1;cy = col/2+1;

Txx = zeros(row,col); Tyy = zeros(row,col);
Txy = zeros(row,col); Txz = zeros(row,col);
Tyz = zeros(row,col);

for I=1:row
    freqx = (I-cx)*wnx;            % angle frequency, freqxu
    for J=1:col
        freqy = (J-cy)*wny;
        freq = sqrt(freqx^2+freqy^2);  
        Txx(I,J) = -freqx^2/freq^2*fz(I,J);  
        Tyy(I,J) = -freqy^2/freq^2*fz(I,J);
        Txy(I,J) = -freqx*freqy/freq^2*fz(I,J);
        Txz(I,J) = 1i*freqx/abs(freq)*fz(I,J);
        Tyz(I,J) = 1i*freqy/abs(freq)*fz(I,J);
    end             
end

Txx = fftshift(Txx); Tyy = fftshift(Tyy);
Txy = fftshift(Txy); Txz = fftshift(Txz); 
Tyz = fftshift(Tyz); 

Txx = ifft2(fillnan(Txx)); Tyy = ifft2(fillnan(Tyy));
Txy = ifft2(fillnan(Txy)); Txz = ifft2(fillnan(Txz));
Tyz = ifft2(fillnan(Tyz));

Txx = Txx(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff);   % delete the other parts
Tyy = Tyy(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff);
Txy = Txy(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff);
Txz = Txz(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff);
Tyz = Tyz(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff);

Txx = real(Txx);                                % keep real part
Tyy = real(Tyy); Txy = real(Txy);
Txz = real(Txz); Tyz = real(Tyz);

%---------------Model parameters of the two prisms---------------%

modxyz1 = [x0-a/sqrt(2) x0+a/sqrt(2);y0-a/sqrt(2) y0+a/sqrt(2);h h+d];
modxyz2 = [x0-a x0+a;y0-a y0+a;h h+d];

obsxyz = [{x},{y}];
Tzz1 = prsm1(obsxyz,modxyz1,rho);
Tzz2 = prsm2(obsxyz,modxyz2,rho);

%-------------------------------Plot------------------------------%

xx = 0:50:50*99;                     % Used for plotting
yy = 0:50:50*99;

figure(1),

subplot(331), imagesc(xx,yy',flipud(rot90(Txx))), axis xy,axis equal,axis tight,title('Txx'),colorbar
subplot(332), imagesc(xx,yy',flipud(rot90(Txy))), axis xy,axis equal,axis tight,title('Txy'),colorbar
subplot(333), imagesc(xx,yy',flipud(rot90(Txz))), axis xy,axis equal,axis tight,title('Txz'),colorbar
subplot(335), imagesc(xx,yy',flipud(rot90(Tyy))), axis xy,axis equal,axis tight,title('Tyy'),colorbar
subplot(336), imagesc(xx,yy',flipud(rot90(Tyz))), axis xy,axis equal,axis tight,title('Tyz'),colorbar
subplot(339), imagesc(xx,yy',flipud(rot90(Tzz))), axis xy,axis equal,axis tight,title('Tzz'),colorbar

colormap(geocolor)

figure(2),

subplot(131), imagesc(xx,yy',flipud(rot90(Tzz1))), axis xy,axis equal,axis tight,title('Small Prism'),colorbar
subplot(132), imagesc(xx,yy',flipud(rot90(Tzz2))), axis xy,axis equal,axis tight,title('Large Prism'),colorbar
subplot(133), imagesc(xx,yy',flipud(rot90(Tzz))), axis xy,axis equal,axis tight,title('Cylinder'),colorbar

colormap(geocolor) 

end

%------------------------Results verification---------------------%
% gravity and gravity gradient anomalies of two prisms are calculated here
% to verify the results of this vertical cylinder.
%
%  Author: Dailei Zhang & Danian Huang
%  Date  : 2016/12/20
%  Email : zhangdailei@hotmail.com 

function Y = prsm1(obsxyz,modxyz1,rho)

x = obsxyz{1}; [m,n] = size(x); x = x(1,:);
y = obsxyz{2}; y = y(:,1);

Tzz = zeros(m,n);

for i=1:n
    for j=1:n
        obs = [x(i) y(j) 0];       % Observation is on the ground, z=0
        Tzz(i,j) = Prsm2Tij('Tzz',obs,modxyz1,rho);
    end
end

Y = Tzz;

end

function Y = prsm2(obsxyz,modxyz2,rho)

x = obsxyz{1}; [m,n] = size(x); x = x(1,:);
y = obsxyz{2}; y = y(:,1);

Tzz = zeros(m,n);

for i=1:n
    for j=1:n
        obs = [x(i) y(j) 0];       % Observation is on the ground, z=0
        Tzz(i,j) = Prsm2Tij('Tzz',obs,modxyz2,rho);
    end
end

Y = Tzz;

end

%------------------------Gravity and gradient calculation------------------
function Tij=Prsm2Tij(action,obsxyz,modxyz,rho)

% function called to calculate Tzz of a prism
%
% Input---
%         rho: residual density, with unit of kg/m^3
%         action: type of calculation of gravity gradient Tij                                         
%         obsxyz: ordinates of a survey point, an array has three elements
%         modxyz: ordinates of two points of the prism, a matrix of dimension of 3x2
%
% Output---
%         Tij: gravity gradient component
%
% Example---
%         Tij = Prsm2Tij('Tzz',[3,4,0],[1,2;1,2;1,2],250)
%  

G=6.67*10^(-11); 
Tij=0;

for i=1:2
   x = obsxyz(1)-modxyz(1,i);
   x2 = x*x;
   
   for j=1:2
      y = obsxyz(2)-modxyz(2,j);
      y2 = y*y;
      
      for k=1:2
         z = -(obsxyz(3)-modxyz(3,k));      % z-point down-ward
         z2 = z*z;
         r = sqrt(x2+y2+z2);
         uijk = (-1)^(i+j+k);
         
         switch action
          
            case 'Txz'
               Ttemp = -log(y+r);
               
            case 'Tyz'
               Ttemp = -log(x+r);
               
            case 'Tzz'
               if(z==0)
                  Ttemp = -pi/2;
               else   
                  Ttemp = -atan((x*y)/(z*r));
               end   
               
            case 'Tyy'
               if (y==0)
                  Ttemp = -pi/2;
               else   
                  Ttemp = -atan((x*z)/(y*r));
               end   
            case 'Txx'
                
               if (x==0)
                  Ttemp = -pi/2;
               else   
                  Ttemp = -atan((y*z)/(x*r));
               end   
               
            case 'Txy'
               Ttemp = log(z+r);
         end                            % end switch action
         Tij = Tij+uijk*Ttemp;
      end
   end
end

Tij = 1e9*G*rho*Tij;                    % 1e9 is the unit conversion factor, from s^-2 to E

end                    

%-------------------------------Edge extension-----------------------------
function y = extend(x,siz)

% function called to do edge extension of data x using cosine method
%
% Input---
%   x   : a time series or a matrix
%   siz : the size of y, i.e., the result of extension;
%         siz should be a scalar if x is a sequence and should be an array
%         with two elements if x is a matrix, i.e., siz = [row col];

% Output---
%   y   : result of extension, whose size is determined by Input siz
%
%  Author: Dailei Zhang
%  Date  : 2015/12/15
%  Email : zhangdailei@hotmail.com 


num = length(size(x));
if num>2
    error('x must not be over two dimension')
end

minSiz = min(size(x));  % check dimension of x, minSiz>1 2D; minSiz=1, 1D

if minSiz==1
    if rem(log2(siz),1)==0  % check if siz is power of 2 by an integer
        % doing 1D edge extension
        xlen = length(x);
        Lx1 = floor((siz-xlen)/2);
        Lx2 = siz-xlen-Lx1;  % extension length on left and right direction
        ext = zeros(siz,1);
        for i=1:Lx1          % left extension
            ext(i) = x(1)*(1+cos(pi*(Lx1+1-i)/Lx1))/2; 
        end
        for j=1:Lx2          % right extension
            ext(xlen+Lx1+j) = x(xlen)*(1+cos((pi*j)/(Lx2)))/2;
        end
        ext(Lx1+1:Lx1+xlen) = x;
    else
        error('siz should be power of 2 by an integer')
    end
    y = ext;
end


if minSiz>1
    if ~(rem(log2(siz(1)),1))&&~(rem(log2(siz(2)),1))==1
       [m,n] = size(x);
       M = siz(1);N=siz(2);
       ext = zeros(M,N);
       Ly1 = floor((M-m)/2);
       Ly2 = M-m-Ly1;
       Lx1 = floor((N-n)/2);
       Lx2 = N-n-Lx1;
       % extension along rows 
       ext_row = zeros(m,N);
       for i=1:m
           for k=1:Lx1            % left extension
               ext_row(i,k) = x(i,1)*(1+cos(pi*(Lx1+1-k)/Lx1))/2;   
           end
           for k=1:Lx2            % right extension
               ext_row(i,n+Lx1+k) = x(i,n)*(1+cos((pi*k)/(Lx2)))/2;
           end
       end
       ext_row(1:m,Lx1+1:Lx1+n) = x;
       
       % extension along columns
       for i=1:N
           for k=1:Ly1            % left extension
               ext(k,i) = ext_row(1,i)*(1+cos(pi*(Ly1+1-k)/Ly1))/2;   
           end
           for k=1:Ly2            % right extension
               ext(m+Lx1+k,i) = ext_row(m,i)*(1+cos((pi*k)/(Ly2)))/2;
           end
       end
       ext(Ly1+1:Ly1+m,:) = ext_row;
    else
        error('siz should both be power of 2 by an integer')
    end
    y = ext;
end

end

%--------------------------------fillnan---------------------------%
function [newim, mask] = fillnan(im)
    
    % Generate distance transform from non NaN regions of the matrix. 
    % L will contain indices of closest non NaN points in the image
    
    mask = ~isnan(im);   
    
    if all(isnan(im(:)))
        newim = im;
        warning('All elements are NaN, no filling possible\n');
        return
    end
    
    [~,L] = bwdist(mask);   
    
    ind = find(isnan(im));  % Indices of points that are NaN
    
    % Fill NaN locations with value of closest non NaN pixel
    newim = im;
    newim(ind) = im(L(ind));
end

%
function [scale]=geocolor()

% Colormap conversion to the matlab format
% The color template used is read from Geosoft
% The read values are 39, the remaining are obtained using interp1 function
%
%  Author: Dailei Zhang
%  Date  : 2016/11/28
%  Email : zhangdailei@hotmail.com 

scale = ...
   [0,0,0.992156862745098;
    0,0.166666666666667,0.996078431372549;
    0,0.333333333333333,1;
    0,0.415686274509804,1;
    0,0.498039215686275,1;
    0,0.582352941176471,1;
    0,0.666666666666667,1;
    0,0.745098039215686,0.996078431372549;
    0,0.823529411764706,0.992156862745098;
    0,0.864705882352941,0.992156862745098;
    0,0.905882352941177,0.992156862745098;
    0,0.949019607843137,0.992156862745098;
    0,0.992156862745098,0.992156862745098;
    0,0.996078431372549,0.888235294117647;
    0,1,0.784313725490196;
    0,0.996078431372549,0.666666666666667;
    0,0.992156862745098,0.549019607843137;
    0,0.992156862745098,0.398039215686275;
    0,0.992156862745098,0.247058823529412;
    0,0.992156862745098,0.215686274509804;
    0,0.992156862745098,0.184313725490196;
    0,0.992156862745098,0.158823529411765;
    0,0.992156862745098,0.133333333333333;
    0,0.992156862745098,0.0666666666666667;
    0,0.992156862745098,0;
    0.137254901960784,0.992156862745098,0;
    0.274509803921569,0.992156862745098,0;
    0.331372549019608,0.992156862745098,0;
    0.388235294117647,0.992156862745098,0;
    0.403921568627451,0.992156862745098,0;
    0.419607843137255,0.992156862745098,0;
    0.490196078431373,0.992156862745098,0;
    0.560784313725490,0.992156862745098,0;
    0.633333333333333,0.992156862745098,0;
    0.705882352941177,0.992156862745098,0;
    0.776470588235294,0.992156862745098,0;
    0.847058823529412,0.992156862745098,0;
    0.919607843137255,0.992156862745098,0;
    0.992156862745098,0.992156862745098,0;
    0.992156862745098,0.949019607843137,0;
    0.992156862745098,0.905882352941177,0;
    0.992156862745098,0.864705882352941,0;
    0.992156862745098,0.823529411764706,0;
    0.996078431372549,0.786274509803922,0;
    1,0.749019607843137,0;
    0.996078431372549,0.723529411764706,0;
    0.992156862745098,0.698039215686275,0;
    0.996078431372549,0.682352941176471,0;
    1,0.666666666666667,0;
    1,0.619607843137255,0;
    1,0.572549019607843,0;
    1,0.554901960784314,0;
    1,0.537254901960784,0;
    0.996078431372549,0.513725490196078,0;
    0.992156862745098,0.490196078431373,0;
    0.992156862745098,0.449019607843137,0;
    0.992156862745098,0.407843137254902,0;
    0.992156862745098,0.366666666666667,0;
    0.992156862745098,0.325490196078431,0;
    0.992156862745098,0.262745098039216,0;
    0.992156862745098,0.200000000000000,0;
    0.992156862745098,0.137254901960784,0;
    0.992156862745098,0.0745098039215686,0;
    0.992156862745098,0.0372549019607843,0;
    0.992156862745098,0,0;
    0.996078431372549,0,0.107843137254902;
    1,0,0.215686274509804;
    1,0,0.321568627450980;
    1,0,0.427450980392157;
    1,0,0.570588235294118;
    1,0,0.713725490196078;
    1,0.0215686274509804,0.784313725490196;
    1,0.0431372549019608,0.854901960784314;
    1,0.258823529411765,0.927450980392157;
    1,0.474509803921569,1;
    0.996078431372549,0.545098039215686,0.996078431372549;
    0.992156862745098,0.615686274509804,0.992156862745098];

end