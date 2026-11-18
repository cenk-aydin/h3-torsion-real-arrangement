-- M2 --script is_decomposable.m2 N Q   (Q a prime with Q = 1 mod N)
needsPackage "HyperplaneArrangements"

n = value scriptCommandLine#1;
q = value scriptCommandLine#2;
R = ZZ/q[x, y, z];
g = 2;
while #(set apply(n, k -> (g^((q - 1) // n * k)) % q)) < n do g = g + 1;
w = (g^((q - 1) // n))_R;

-- the points (u, u^2, 1 + u^3) of x^3 + y^3 = xyz: three are collinear iff u1 u2 u3 = -1
A = arrangement(apply(n, k -> (u := -w^k; u * x + u^2 * y + (1 + u^3) * z)), R);

triples = set apply(select(flats(2, A), F -> #toList F == 3), F -> sort toList F);
zeroSum = set select(subsets(n, 3), s -> sum s % n == 0);
print("n = " | toString n | " over ZZ/" | toString q
    | ": flat sizes " | toString tally apply(flats(2, A), F -> #toList F)
    | ", triple points = zero-sum triples: " | toString(triples === zeroSum));

I = orlikSolomon(A, QQ);
b = betti res(coker vars((ring I) / I), LengthLimit => 3);
phi3 = (3 * b_(3, {3}, 3) - 3 * b_(1, {1}, 1) * b_(2, {2}, 2) + b_(1, {1}, 1)^3 - b_(1, {1}, 1)) // 3;
print("n = " | toString n | ": phi_3 = " | toString phi3 | ", local sum 2t = " | toString(2 * #zeroSum)
    | ", isDecomposable: " | toString isDecomposable A);
