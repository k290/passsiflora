\\ The real results of the 31 work units as coefficient vectors, and their shifted copies x -> x+t, |t| <= 6.
\\ Reads hits.txt (one polynomial per line, as in the truth files), writes hits_vec.txt and hits_shift_vec.txt
\\ (11 integers per line, leading coefficient first).      gp -q hits.gp < /dev/null
default(parisize, 100000000);
{
my(P = readstr("hits.txt"), f, v);
for (i = 1, #P,
  f = eval(P[i]);
  v = Vec(f);
  write("hits_vec.txt", strjoin(apply(c -> Str(c), v), " "));
  for (t = -6, 6,
    if (t == 0, next);
    v = Vec(subst(f, x, x + t));
    write("hits_shift_vec.txt", strjoin(apply(c -> Str(c), v), " "))));
print(#P, " results, ", 12 * #P, " shifted copies");
}
