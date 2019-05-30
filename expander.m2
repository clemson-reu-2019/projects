
-- This expands the coefficients for Ok++ --

p = 14 -- or put whatever you want here
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

-- for a diag, use these --
R = ZZ[x,y] / (ideal(x, y))^p
expand_diag_allpos = (k,d) -> expand_box_allpos(k,k,d)

-- The following expands the coefficients for Ok+ --

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

-- This calculates the number of partitions with less than a certain number of terms --

m = 30 -- or put whatever you want
n = 20
R = ZZ[x,y] / (ideal(x^m, y^n))



-- expandlessthan:
-- calculate the expansion for p_r(n) (where r = a + b*sqrt(d))
--   that is, the number of partitions with parts
--   less than r.
--
-- Parameters:
-- a - the integer part of r
-- b - the quadratic part of r
-- maxA - how many integer-part terms to compute
-- maxB - how many quadratic-part terms to compute
-- d - the square root of the quadratic field
--
-- Example: 
-- expandlessthan(3,1,30,20,2)
-- would caluclate the expansion for p_r(n)
--   with terms up to 30 + 20*sqrt(2), 
--   with each part being less than 3 + 1*sqrt(2)

--- START COPYING HERE

expandgeomseries = (k,x,N) -> (
  g = 1;
  for i from 1 to N do (
    g = g + x^i;
  );
  k*g
)

expandlessthan = (a,b,maxA,maxB,d) -> (
  R := ZZ[x,y] / (ideal(x^(maxA+1),y^(maxB+1)));
  f = 1;
  for i from 1 to a do (
    f = f * expandgeomseries(1,x^i,max(maxA,maxB));
  );

  for i from 1 to b do (
    for j from ceiling(i * sqrt(d)) to a do (
      f = f * expandgeomseries(1,(x^j)*(y^i),max(maxA,maxB));
    );
  );
  f
)

--- STOP COPYING HERE
