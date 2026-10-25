# classify.awk -- what the x87 filter computes for lane 1, counted from the standard error of the instrumented feeder
# (repro.sh, part inst) run with FEED_ALONE=1, so that lane 1 is the polynomial at every stage:
#   LDS/LDD line   the ten loaded coefficients: "exact" if all ten are integers, "noisy" otherwise
#   EU line        the Euler loop: the first squaring at which base is an integer (-1 never, 0 the final constant itself),
#                  acc at that moment, the final acc
#   CBS/CBD/CG     the chain's last line; a chain_block run with ok = F is redone by chain_general and is not counted
# Columns: chain runs; final acc = -1; final acc = +1; of the runs ending at -1 or +1: final constant an integer, base an
# integer only from squaring 1 on with acc still 1 at that moment, neither; runs that reject the polynomial.
#     awk -f classify.awk trace.txt
/^LD[SD]/ { prec = (substr($0, 3, 1) == "S") ? "single" : "double"; load = "exact"
            for (i = 3; i <= NF; i++) if ($i + 0 != int($i + 0)) load = "noisy"
            next }
/^EU/     { first = $4 + 0; accthen = $5 + 0; fin = $6 + 0; next }
/^CB[SD]/ { if (substr($0, 41, 1) != "T") next; rej = (fin == -1); count(); next }
/^CG  /   { rej = ($6 + 0 > 0 && $7 + 0 == 0 && fin == -1); count(); next }
function count(  k) {
  k = prec " " load; seen[k] = 1; runs[k]++
  if (fin == -1) m1[k]++
  if (fin == 1) p1[k]++
  if (fin == 1 || fin == -1) {
    if (first == 0) cint[k]++
    else if (first >= 1 && accthen == 1) late[k]++
    else other[k]++
  }
  if (rej) rj[k]++
}
END {
  printf "%-6s %-5s %6s %6s %6s %10s %10s %6s %7s\n", "prec", "load", "runs", "acc=-1", "acc=+1", "const.int", "base.later", "other", "rejects"
  split("single exact,single noisy,double exact,double noisy", order, ",")
  for (i = 1; i <= 4; i++) { k = order[i]; if (!(k in seen)) continue; split(k, w, " ")
    printf "%-6s %-5s %6d %6d %6d %10d %10d %6d %7d\n", w[1], w[2], runs[k], m1[k], p1[k], cint[k], late[k], other[k], rj[k] }
}
