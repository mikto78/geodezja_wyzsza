clear;
% Queenstown ---> Nowa Zelandia
%phi = -45.0302300;
%lambda = 168.6627100;

% Libreville ---> Gabon
%phi = 0.3924100
%lambda = 9.4535600;

% Rovaniemi ---> Finlandia
phi = 66.5000000;
lambda = 25.7166700;

% rektascensja i deklinacja Alfa Arietis z gwiazdozbioru Barana
r = 2 + 7/60 +10.406/3600
d = 23 +27/60 + 44.7/60

%kąt w kolejnych godzinach
h = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
kat_h = katgodz(2001, 4, 7, h, lambda, r)

for i = 1:size(kat_h, 1)
    if kat_h(i) > 360
        kat_h(i) = kat_h(i) - 360;
    end
end

zenit = odl_zenitalna(phi, d, kat_h)
Az = azymut(phi, d, kat_h)

%wsp xyz
x = sind(zenit).*cosd(Az);
y = sind(zenit).*sind(Az);
z = cosd(zenit);

%wizualizacja położenia gwiazdy na kuli
[X,Y,Z] = sphere(16);
X = X(1:end,:);
Y = Y(1:end,:);
Z = Z(1:end,:);
surf(X,Y,Z,'FaceColor','blue','FaceAlpha',0.4)
axis equal, hold on;
scatter3(x,y,z, 150, 'yellow', '*')

function g = GMST(JD)
T = (JD - 2451545)/36525;
g = 280.46061837 + 360.98564736629 * (JD - 2451545) + 0.000387933*T.^2 - T.^3/38710000;
g = mod(g,360);
end

function t = katgodz(y,m,d,h,lambda,alfa) 
time = datetime(y,m,d);
jd = juliandate(time); 
g = GMST(jd); 
UT1 = h * 1.002737909350795; 
S = UT1 * 15 + lambda + g; 
t = S - alfa * 15; 
end  

%odległość zenitalna
function z = odl_zenitalna(fi, dek, t)
cosz = sind(fi) .* sind(dek) + cosd(fi) .* cosd(dek) .* cosd(t);
z = acosd(cosz)
end

%azymuty gwiazd
function Az = azymut(fi, dek, t)
x = -cosd(dek) .* sind(t)
y = cosd(fi) .* sind(dek) - sind(fi) .* cosd(dek) .* cosd(t)
Az = atan2d(x, y) 
if Az > 360
    Az = Az - 360
elseif Az < 0
    Az = Az + 360
end
end


