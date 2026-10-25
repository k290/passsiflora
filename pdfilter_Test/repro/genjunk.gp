\\ junk: random monic degree-10 polynomials that do NOT pass (the filter should reject nearly all of them)
default(parisize, 100000000);
setrand(77);
N = 5000;
passes(f) = {
  my(D = abs(poldisc(f)));
  if (D == 0, return(1));
  while (D % 2 == 0, D \= 2);
  while (D % 5 == 0, D \= 5);
  issquare(D)
}
main() = {
  my(f, v, k = 0);
  while (k < N,
    f = x^10 + sum(i = 0, 9, (random(201) - 100) * x^i);
    if (passes(f), next);
    k++;
    v = Vec(f);
    write("junk.txt", concat(vector(11, j, Str(v[j], " ")))));
  print(N, " junk polynomials written to junk.txt");
}
main();
