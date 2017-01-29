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

% DRAW STATEMENTS

pair pa, pb, pc, pd, pe, pf;
pa := (-0.366 u, -0.576 u);
pb := (+1.589 u, -0.944 u);
pc := (+1.346 u, +0.981 u);
pd := (-0.295 u, +1.788 u);
pe := (-1.238 u, +0.735 u);
pf := (-1.441 u, -1.050 u);
path poly;
poly := pa--pb--pc--pd--cycle;
fill poly withcolor lblenderorange;
pickup pencircle scaled 2 pt;
draw poly;
draw pa--pf;
draw pd--pe;

fill fullcircle scaled (8 pt) shifted pa;
fill fullcircle scaled (8 pt) shifted pb;
fill fullcircle scaled (8 pt) shifted pc;
fill fullcircle scaled (8 pt) shifted pd;
fill fullcircle scaled (8 pt) shifted pe;
fill fullcircle scaled (8 pt) shifted pf;

% END

endfig;
end
