beginfig(-1)

lw := 0.5 pt;
ahw := 2 * ahlength * sind(ahangle / 2);
ahl := ahlength * cosd(ahangle / 2);

red := (0.7, 0, 0);
green := (0, 0.5, 0);
blue := (0, 0, 0.7);
color magenta, cyan, orange, violet;
magenta := (0.7, 0, 0.7);
cyan := (0, 0.5, 0.7);
orange := (0.91, 0.494, 0);
violet := (0.349, 0, 0.702);
color lred, lgreen, lblue, lmagenta, lcyan, lorange, lviolet;
lred        := 0.7 white + 0.3 red;
lgreen      := 0.7 white + 0.3 green;
lblue       := 0.7 white + 0.3 blue;
lmagenta    := 0.7 white + 0.3 magenta;
lcyan       := 0.7 white + 0.3 cyan;
lorange     := 0.7 white + 0.3 orange;
lviolet     := 0.7 white + 0.3 violet;

color blenderorange, lblenderorange;
blenderorange := (1.000, 0.522, 0);
lblenderorange := 0.5 white + 0.5 blenderorange;

evenscale := 0.5;
dotscale := 0.3;

vardef arcsind primary x = angle((1+-+x,x)) enddef;
vardef arccosd primary x = angle((x,1+-+x)) enddef;

% CONFIG VARS
u := 3 cm;
l := 0.3 u;
w := 0.05 u;

% DRAW STATEMENTS

pair pa, pb, pc, pd, pe;
pa := (-0.339 u, -1.397 u);
pb := (+1.466 u, -1.052 u);
pc := (+1.466 u, +1.052 u);
pd := (-0.339 u, +1.397 u);
pe := (-1.327 u, +0.954 u);
path poly;
poly := pa--pb--pc--pd--cycle;
fill poly withcolor lblenderorange;
pickup pencircle scaled 2 pt;
draw poly;
draw pd--pe;

draw (((-l / 2, 0)--(+l / 2, 0)) shifted (0, +w / 2)) rotated 45 shifted 0.5[pa, pd];
draw (((-l / 2, 0)--(+l / 2, 0)) shifted (0, -w / 2)) rotated 45 shifted 0.5[pa, pd];
draw (((-l / 2, 0)--(+l / 2, 0)) shifted (0, +w / 2)) rotated 45 shifted 0.5[pb, pc];
draw (((-l / 2, 0)--(+l / 2, 0)) shifted (0, -w / 2)) rotated 45 shifted 0.5[pb, pc];

fill fullcircle scaled (8 pt) shifted pa;
fill fullcircle scaled (8 pt) shifted pb;
fill fullcircle scaled (8 pt) shifted pc;
fill fullcircle scaled (8 pt) shifted pd;
fill fullcircle scaled (8 pt) shifted pe;

% END

endfig;
end
