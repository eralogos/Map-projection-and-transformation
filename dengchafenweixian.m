function varargout = dengchafenweixian(varargin)
%VGRINT1  Van Der Grinten I Polyconic Projection
%
%  In this projection, the world is enclosed in a circle.  Scale is true
%  along the Equator and increases rapidly away from the Equator. Area
%  distortion is extreme near the poles.  This projection is neither
%  conformal nor equal area.
%
%  This projection was presented by Alphons J. van der Grinten 1898. He
%  obtained a U.S. patent for it in 1904.  It is also known simply as the
%  Van der Grinten projection (without a "I").

% Copyright 1996-2011 The MathWorks, Inc.

mproj.default = @dengchafenweixianDefault;
mproj.forward = @dengchafenweixianFwd;
mproj.inverse = @dengchafenweixianInv;
mproj.auxiliaryLatitudeType = 'geodetic';
mproj.classCode = 'Poly';

varargout = applyProjection(mproj, varargin{:});

%--------------------------------------------------------------------------

function mstruct = dengchafenweixianDefault(mstruct)

% The trimlon values below are pulled in by eps(180) to keep
% diff(trimlon) < 360, which ensures that subfunction adjustFrameLimits
% in private/resetmstruct.m will clamp the frame longitude limit to
% trimlon.  That's appropriate for this special projection that is
% intended to display the entire earth.  If the FLonLimit interval is
% not forced to be a subset of [-180 180], the projection may not be
% one-to-one.  Near one limit or the other, two different points
% (in the geographic system) could project to the same point in the map
% plane, resulting in a map that appears to fold back onto itself.

mstruct.mapparallels = [];
mstruct.nparallels   = 0;
mstruct.fixedorient  = [];
[mstruct.trimlat, mstruct.trimlon] = fromDegrees( ...
    mstruct.angleunits, [-90 90], [-180, 180]);

%--------------------------------------------------------------------------

function [x, y] = dengchafenweixianFwd(mstruct, lat, lon)

%epsilon = 10*epsm('radians');
a = ellipsoidprops(mstruct);

%定义基本方程和常数 
xn=100.8721.*lat-18.3126.*lat.^2;
yn=165-9.4729.*lat-28.8595.*lat.^2;
lonn=pi;
b=1.1;
c=0.0005050505;

%中央经线x0函数式
x0=63.41514.*lat+0.9404646.*lat.^3;
%计算各纬线的动半径p 中央经线,右边经线与各纬度线交点的动径角an 投影之后的动径角ai
p=(yn.^2+(xn-x0).^2)./2./(xn-x0);
an=asin(yn./abs(p));
ai=an./lonn.*b.*(1-c.*lon).*lon;
%计算赤道的特殊情况
yl=yn./lonn.*b.*(1-c.*lon).*lon;
%求x和y注意调转xy
q=p+x0;

if lat==0
    x=yl;
    y=0;
else
    y=q-p.*cos(ai);
    x=abs(p).*sin(ai);
end



%--------------------------------------------------------------------------

function [lat, lon] = dengchafenweixianInv(mstruct, x, y)

a = ellipsoidprops(mstruct);

% Normalize by the radius

x = x / (pi*a);
y = y / (pi*a);

% Pick up NaN place holders and points at (0,0)

lon = x;
lat = y;

% Process points not at (0,0)

indx1 = find(x ~= 0 | y ~= 0);
indx2 = find(x ~= 0);

% Inverse transformation

if ~isempty(indx1)

    fact1 = x(indx1).^2 + y(indx1).^2;
    c1 = -abs(y(indx1)) .* (1 + fact1);
    c2 = c1 - 2* y(indx1).^2 + x(indx1).^2;
    c3 = -2*c1 + 1 + 2* y(indx1).^2 + fact1.^2;

    d = y(indx1).^2 ./ c3 + (2*c2.^3./c3.^3 - 9*c1.*c2./c3.^2)/27;
    a1 = (c1 - c2.^2./(3*c3))./c3;
    m1 = 2*sqrt(-a1/3);
    cos_theta1_times3 = 3*d ./ (a1.*m1);
    % Correct for possible round off/noise
    cos_theta1_times3(cos_theta1_times3 < -1) = -1;
    cos_theta1_times3(cos_theta1_times3 >  1) =  1;
    theta1 = acos(cos_theta1_times3)/3;
    lat(indx1) = pi * sign(y(indx1)) .* ...
        (-m1.*cos(theta1+pi/3)-c2./(3*c3));

    % Points at non-zero longitude

    if ~isempty(indx2)
        c1 = x(indx2).^2 + y(indx2).^2;
        c2 = x(indx2).^2 - y(indx2).^2;
        c3 = (c1 - 1 + sqrt(1 + 2*c2 + c1.^2));
        lon(indx2) = pi * c3 ./ (2*x(indx2));
    end
end
