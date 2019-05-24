p = 14
R = ZZ[x,y] / ideal(x^p, y^p)

expand_part = (k,d) -> (
f = 1;
for i from 1 to k do f = f * (1 - x^i);

for i from 1 to k do
  for j from ceiling(i * sqrt(d)) to k do
    f = f * (1 - x^j * y^i);
f
)