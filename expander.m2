p = 14
R = ZZ[x,y] / (ideal(x^a, y^b))

expandboxpos = (a,b,d) -> (
f = 1;
for i from 1 to a do f = f * (1 - x^i);

for i from 1 to b do
  for j from ceiling(i * sqrt(d)) to a do
    f = f * (1 - x^j * y^i);
f
)

expansion = expand_diag_allpos(a,b,2)
listForm(expansion)

-- for a box, use these --
R = ZZ[x,y] / (ideal(x, y))^p
expand_diag_allpos = (k,d) -> expand_box_allpos(k,k,d)

-- For Ok+ only --

expandbox = (a,b,d) -> (
R := QQ[x,y, Inverses=>true,MonomialOrder=>Lex];
use R;
f := 1;

for i from 1 to a do (
  f = f * (1 - x^i);
  f = part(0,a-1,{1,0},f)
);

for i from 1 to b do (
  for j from ceiling(i * sqrt(d)) to a do (
    f = f * (1 - (x^j)*(y^i))*(1 - (x^j)*(y^-i));
    f = part(0,a-1,{1,0},f);
    bnd = b + ceiling(a * sqrt(d));
    f = part(-(bnd-1),bnd-1,{0,1},f); 
));

f = part(-(b-1),b-1,{0,1},f);

f
)
