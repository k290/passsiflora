\\ valid test polynomials like genPoly() of test_pdfilter_linux.cpp: f = G(x^2), G = x^5 + c4 x^4 + ... + c0, irreducible,
\\ disc(f) = (2^a 5^b) * square; x^8 coefficient zeroed in a quarter, x -> x+t (|t| <= 30) in three quarters.
\\ usage: gp -q genpolys.gp  (N and the seed are set below)   -> writes polys.txt, one polynomial per line, leading coefficient first
default(parisize, 200000000);
setrand(20261025);
N = 20000;

valid(f) = {
  my(D = abs(poldisc(f)));
  if (D == 0, return(0));
  while (D % 2 == 0, D \= 2);
  while (D % 5 == 0, D \= 5);
  issquare(D)
}

gen() = {
  my(f, c, c0, t);
  until (polisirreducible(f),
    c = vector(5, i, random(100));            \\ c[1] = c0, ..., c[5] = c4
    if (random(4) == 0, c[5] = 0);
    c0 = c[1];
    if (c0 >= 55 && c0 < 70, c0 = -1, if (c0 >= 70 && c0 < 85, c0 = 2, if (c0 >= 85, c0 = 5, c0 = c0^2)));
    f = x^10 + c[5]*x^8 + c[4]*x^6 + c[3]*x^4 + c[2]*x^2 + c0);
  if (random(4) != 0, t = 1 + random(30); if (random(2), t = -t); f = subst(f, x, x + t));
  f
}

main() = {
  my(f, bad = 0, v);
  for (i = 1, N,
    f = gen();
    if (!valid(f), bad++);
    v = Vec(f);
    write("polys.txt", concat(vector(11, k, Str(v[k], " ")))));
  print(N, " polynomials written to polys.txt, ", bad, " not valid");
}
main();
