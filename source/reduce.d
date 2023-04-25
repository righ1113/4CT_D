module reduce;

import std.file   : readText;
import std.json   : JSONValue, parseJSON;
import std.stdio  : writef, writefln;
import std.string : format;

const int    VERTS   = 27;	/* max number of vertices in a free completion + 1 */
const int    DEG     = 13;	/* max degree of a vertex in a free completion + 1, must be at least 13 because of row 0 */
const int    EDGES   = 62;	/* max number of edges in a free completion + 1 */
const int    MAXRING = 14;	/* max ring-size */
const long[] SIMATCHNUMBER
  = [0L, 0L, 1L, 3L, 10L, 30L, 95L, 301L, 980L, 3_228L, 10_797L, 36_487L, 124_542L, 428_506L, 1_485_003L];

//alias tp_gconf  = long[DEG][VERTS];
alias tp_angle  = long[5][EDGES];
alias tp_edgeno = long[EDGES][EDGES];


// --------------------------------------------------------------------------------------------------------------------
void reduce()
{
  long verts, ring, nlive, ncodes, i, nbyte, count, edges, nReal;
  long[17] power;
  long[EDGES + 1] contract;
  tp_angle angle, diffangle, sameangle;
  byte[] live, real2;
  tp_edgeno edgeno;
  long[string] gConfMajor;

  power[1] = 1;
  for (i = 2; i < 17; i++) power[i] = 3 * power[i - 1];	/* power[i] = 3^(i-1) for i>0 */
  ncodes = (power[MAXRING] + 1) / 2;       live  = new byte[ncodes];
  nbyte  = SIMATCHNUMBER[MAXRING] / 8 + 2; real2 = new byte[nbyte];

	string fileName = "data/r_good_confs.json";
  auto jv2 = readText(fileName).parseJSON();
  // writeln(jv2[1][0][1]); // => 122
  JSONValue[] jarr = jv2.array();

  for (count = 0; count < jarr.length; count++) {
    writefln("%d", count);
    if (count == 5) break;

    // 変数のリセット
    contract[]  = 0;
    angle[]     = [0, 0, 0, 0 ,0];
    diffangle[] = [0, 0, 0, 0 ,0];
    sameangle[] = [0, 0, 0, 0 ,0];
    live[]      = 0;
    real2[]     = 0;
    edgeno[] = [0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L,
      0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L ,0L, 0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L, 0L ,0L, 0L ,0L];

    // gConf でよく使うものをまとめた
    auto gConf = jarr[count];
    verts = gConf[0+1][0].integer;
    ring  = gConf[0+1][1].integer;
    assert(ring <= MAXRING, "E43 Ring-size bigger than");
    edges = 3 * verts - 3 - ring;
    gConfMajor["verts"]  = verts;
    gConfMajor["ring"]   = ring;
    gConfMajor["term"]   = 3 * (verts - 1) - ring;
    gConfMajor["edges"]  = edges;
    gConfMajor["claim"]  = gConf[0+1][2].integer;
    gConfMajor["cont0"]  = gConf[0+2][0].integer;
    gConfMajor["contE"]  = gConf[0+1][3].integer;
    gConfMajor["bigno"]  = (power[ring + 1] - 1) / 2;
    gConfMajor["ncodes"] = (power[ring] + 1) / 2;

    // 1.
    getEdgeNo(gConf, edgeno, gConfMajor);

    // 2.
    findangles(gConf, angle, diffangle, sameangle, contract, edgeno, gConfMajor);

    // 3.
    ncodes = gConfMajor["ncodes"];
    for (i = 0; i < ncodes; i++) live[i] = 1;
    nlive = findlive(live, angle, power, gConfMajor);

    // 4.
    nbyte = SIMATCHNUMBER[ring] / 8 + 1;
    for (i = 0; i <= nbyte; i++) real2[i] = -1;
    do
      nReal = dReduceTestMatch(ring, real2, power, live, nbyte);
    while (isUpdateLive(live, ncodes, nlive, nReal));

    // 5.
    checkCReduce(live, nlive, diffangle, sameangle, contract, power, gConfMajor);
  }
  writefln("Reducibility of %d configurations verified", count);

}
private: // ここ以下すべて private


// --------------------------------------------------------------------------------------------------------------------
/* Numbers edges from 1 up, so that each edge has as many later edges in
 * triangles as possible; the ring edges are first.  edgeno[u][v] will be the
 * number of the edge with ends u,v if there is such an edge and 0 otherwise. */
void getEdgeNo(JSONValue gConf, ref tp_edgeno edgeno, long[string] gConfMajor) pure {
  long d, h, u, v, w, x, verts = gConfMajor["verts"], ring = gConfMajor["ring"], term = gConfMajor["term"];
  long maxint, maxes;
  long[VERTS] max, done;
  long inter, maxdeg, best, first, previous;
  auto grav = gConf[best];

  for (v = 1; v <= ring;         v++) { u = (v > 1) ? v - 1 : ring; edgeno[u][v] = edgeno[v][u] = v; }
  for (x = ring + 1; x <= verts; x++) {
    // First we find all vertices from the interior that meet the "done" vertices in an interval,
    // and write them in max[1] .. max[maxes]
    for (maxint = maxes = 0, v = ring + 1; v <= verts; v++) {
      if      (done[v])         { continue; } inter = inInterval(gConf[v+2].array, done);
      if      (inter >  maxint) { maxint = inter; maxes = 1; max[1]       = v; }
      else if (inter == maxint) {                            max[++maxes] = v; }
    }
    // From the terms in max we choose the one of maximum degree
    for (maxdeg = 0, h = 1; h <= maxes; h++) {
      d = gConf[max[h]+2][0].integer; if (d > maxdeg) { maxdeg = d; best = max[h]; }
    }
    grav = gConf[best+2];
    d = grav[0+1].integer; first = 1; previous = done[grav[d+1].integer];
    while ((previous) || (!done[grav[first+1].integer])) {
      previous = done[grav[1+first++].integer]; if (first > d) { first = 1; break; }
    }
    for (h = first; done[grav[h+1].integer]; h++) {
      edgeno[best][grav[h+1].integer] = term;
      edgeno[grav[h+1].integer][best] = term--; if (h == d)    { if (first == 1) break; h = 0; }
    }
    done[best] = 1;
  } // This eventually lists all the internal edges of the configuration

  // Now we must list the edges between the interior and the ring
  for (x = 1; x <= ring; x++) {
    maxint = 0;
    for (v = 1; v <= ring; v++) {
      if (done[v])        { continue; }
      u = (v > 1) ? v - 1 : ring; w = (v < ring) ? v + 1 : 1;
      inter = 3 * gConf[v+2][0+1].integer + 4 * (done[u] + done[w]);
      if (inter > maxint) { maxint = inter; best = v; }
    }
    grav = gConf[best+2]; u = (best > 1) ? best - 1 : ring;
    if (done[u]) { // C では、 int型 の 0 のみが偽となり、それ以外が全て真として扱われる。
      for (h = grav[0+1].integer - 1; h >= 2; h--) {
        edgeno[best][grav[h+1].integer] = term;
        edgeno[grav[h+1].integer][best] = term--;
      }
    } else {
      for (h = 2; h < grav[0+1].integer; h++) {
        edgeno[best][grav[h+1].integer] = term;
        edgeno[grav[h+1].integer][best] = term--;
      }
    }
    done[best] = 1;
  }
}

// if grav meets the done vertices in an interval of length >=1, it returns the length of the interval,
// and otherwise returns 0
long inInterval(JSONValue[] grav, long[] done) pure {
  long d = grav[0+1].integer, j, first, last, worried = 0, length;

  for (first = 1; (first < d) && (!done[grav[first+1].integer]); first++){}
  if  (first == d) { return (done[grav[d+1].integer]); }
  for (last = first; (last < d) && (done[grav[last + 1+1].integer]); last++){} length = last - first + 1;
  if  (last  == d) { return (length); }
  if  (first > 1)  { for (j = last + 2; j <= d; j++) { if (done[grav[j+1].integer]) return (0L); } return (length); }
  for (j = last + 2; j <= d; j++) {
    if (done[grav[j+1].integer]) { length++; worried = 1; } else if (worried) { return (0L); }
  }
  return (length);
}


// --------------------------------------------------------------------------------------------------------------------
/* writes into angle[i] all edges with number >i on a common triangle T say
 * with edge i; and if there is a contract X given, and i is not in X, writes
 * into diffangle[i] all such edges such that no edge of T is in X, and
 * writes into sameangle[i] all such edges not in X so that the third edge of
 * T is in X. Sets contract[i] to 1 if edge number i is in X and to zero
 * otherwise, checks that X is sparse, and if |X|=4 checks that X has a triad */
void findangles(JSONValue gConf, ref tp_angle angle, ref tp_angle diffangle, ref tp_angle sameangle,
  ref long[EDGES+1] contract, tp_edgeno edgeno, long[string] gConfMajor) pure {

  long a, b, c, h, i, j, u, v, w, verts = gConfMajor["verts"], ring = gConfMajor["ring"], edges = gConfMajor["edges"];
  long[VERTS] neighbour;

  assert(edges < EDGES, format("Configuration has more than %d edges", EDGES - 1));

  contract[0] = gConfMajor["cont0"];	// 4 -> 0, number of edges in contract
  assert((contract[0] >= 0 && contract[0] <= 4), "        ***  ERROR27: INVALID CONTRACT  ***");
  // for (i = 5; i <= 2 * contract[0] + 4; i++)
  //   if (gConf[0+1][i].integer < 1 || gConf[0+1][i].integer > gConf[0+1][0].integer)
  //     assert(0, "  ERROR29: ILLEGAL CONTRACT  ***");
  contract[EDGES] = gConfMajor["contE"];

  for (i = 1; i <= contract[0]; i++) {
    u = gConf[0+2][2 * i - 1].integer;
    v = gConf[0+2][2 * i    ].integer;
    if (edgeno[u][v] < 1) { assert(0, "         ***  ERROR29: CONTRACT CONTAINS NON-EDGE  ***"); }
    contract[edgeno[u][v]] = 1;
  }
  for (i = 1; i <= ring; i++)
    if (contract[i]) { assert(0, "         ***  ERROR21: CONTRACT IS NOT SPARSE  ***"); }

  diffangle[0][0] = angle[0][0] = verts;
  diffangle[0][1] = angle[0][1] = ring;
  diffangle[0][2] = angle[0][2] = edges;
  for (v = 1; v <= verts; v++) {
    for (h = 1; h <= gConf[v+2][0+1].integer; h++) {
      if ((v <= ring) && (h == gConf[v+2][0+1].integer)) continue;
      i = (h < gConf[v+2][0+1].integer) ? h + 1 : 1; u = gConf[v+2][h+1].integer; w = gConf[v+2][i+1].integer;
      a = edgeno[v][w]; b = edgeno[u][w]; c = edgeno[u][v];
      if (contract[a] && contract[b]) { assert(0, "         ***  ERROR22: CONTRACT IS NOT SPARSE  ***"); }
      findAnglesSub(a, b, c, angle, diffangle, sameangle, contract);
      findAnglesSub(b, a, c, angle, diffangle, sameangle, contract);
    }
  }

  /* checking that there is a triad */
  if  (contract[0] < 4) return;
  for (v = ring + 1; v <= verts; v++) { /* v is a candidate triad */
    for (a = 0, i = 1; i <= gConf[v+2][0+1].integer; i++) { u = gConf[v+2][i+1].integer;
      for (j = 1; j <= 8; j++) { if (u == gConf[0+2][j].integer) { a++; break; } }
    }
    if  (a < 3)                                    { continue; }
    if  (gConf[v+2][0].integer >= 6)               { return; }
    for (u = 1; u <= verts; u++)                   { neighbour[u] = 0; }
    for (i = 1; i <= gConf[v+2][0+1].integer; i++) { neighbour[gConf[v+2][i].integer] = 1; }
    for (j = 1; j <= 8;          j++)              { if (!neighbour[gConf[0+2][j].integer]) return; }
  }
  assert(0, "         ***  ERROR28: CONTRACT HAS NO TRIAD  ***");
}

void findAnglesSub(long x, long y, long c, ref tp_angle angle, ref tp_angle diffangle, ref tp_angle sameangle,
  long[EDGES+1] contract) pure {
    long d, e;

  if (x > c) {
    d = angle[c][0] >= 4 ? 4 : ++angle[c][0];
    angle[c][d] = x;
    if ((!contract[x]) && (!contract[y]) && (!contract[c])) {
      e = diffangle[c][0] >= 4 ? 4 : ++diffangle[c][0];
      diffangle[c][e] = x;
    }
    if (contract[y]) {
      e = sameangle[c][0] >= 4 ? 4 : ++sameangle[c][0];
      sameangle[c][e] = x;
    }
  }
}


// --------------------------------------------------------------------------------------------------------------------
/* computes {\cal C}_0 and stores it in live. That is, computes codes of
 * colorings of the ring that are not restrictions of tri-colorings of the
 * free extension. Returns the number of such codes */
long findlive(byte[] live, tp_angle angle, long[] power, long[string] gConfMajor) {
  long j, i, u, ring = gConfMajor["ring"], edges = gConfMajor["edges"], extentclaim = gConfMajor["claim"];
  long[] am;
  long extent, bigno = gConfMajor["bigno"], ncodes = gConfMajor["ncodes"], ret;	/* needed in "chgLive" */
  long[EDGES] c, forbidden;	/* called F in the notes */

  c[edges] = 1; j = edges - 1; c[j] = 2; forbidden[j] = 5;
  for (extent = 0;;) {
    while (0 != (forbidden[j] & c[j])) {
      ret = isEndFL(j, c, edges, ring, ncodes, extent, extentclaim); if (ret) return ret;
    }
    if (j == ring + 1) {
      chgLive(c, power, ring, angle, live, extent, bigno);
      ret = isEndFL(j, c, edges, ring, ncodes, extent, extentclaim); if (ret) return ret;
    } else {
      am = angle[--j]; // 前置きの場合は、先にデクリメントをおこなう
      c[j] = 1;
      for (u = 0, i = 1; i <= am[0]; i++) {
        u |= c[am[i]];
      }
      forbidden[j] = u;
    }
  }
}

long isEndFL(ref long j, ref long[EDGES] c, long edges, long ring, long ncodes, long extent, long extentclaim) {
  c[j] <<= 1;
  while (c[j] & 8) {
    if (j >= edges - 1) {
      printStatus(ring, ncodes, extent, extentclaim);
      return (ncodes - extent); // 0 にはならないはず
    }
    c[++j] <<= 1;
  }
  return (0L);
}

void printStatus(long ring, long totalcols, long extent, long extentclaim) {
  writef("\n\n   This has ring-size %d, so there are %d colourings total,\n", ring, totalcols);
  writef("   and %d balanced signed matchings.\n", SIMATCHNUMBER[ring]);
  writef("\n   There are %d colourings that extend to the configuration.", extent);
  assert(extent == extentclaim, "   *** ERROR31: DISCREPANCY IN NUMBER OF EXTENDING COLOURINGS ***");
  writef("\n\n            remaining               remaining balanced\n");
  writef("           colourings               signed matchings\n");
  writef("\n              %7d", totalcols - extent);
}

/* Given a colouring specified by a 1,2,4-valued function "col", it computes
 * the corresponding number, checks if it is in live, and if so removes it. */
void chgLive(long[] col, long[] power, long ring, tp_angle angle, byte[] live, ref long p, long bigno) pure {
  long colno, sum, i, min, max, w;
  long[5] weight;

  for (i = 1; i <= ring; i++) {
    sum = 7 - col[angle[i][1]] - col[angle[i][2]];
    //if (sum < 0) sum = 0;
    //if (sum > 4) sum = 4;
    weight[sum] += power[i];
  }
  min = max = weight[4];
  for (i = 1; i <= 2;    i++) {
    w = weight[i];
    if      (w < min) min = w;
    else if (w > max) max = w;
  }
  colno = bigno - 2 * min - max;
  if (live[colno]) { p++; live[colno] = 0; }
}


// --------------------------------------------------------------------------------------------------------------------
/* runs through "live" to see which colourings still have `real' signed
 * matchings sitting on all three pairs of colour classes, and updates "live"
 * accordingly; returns 1 if nlive got smaller and stayed >0, and 0 otherwise */
bool isUpdateLive(ref byte[] live, long ncols, ref long p, long nReal) {
  long i, nlive = p, newnlive = 0;

  if  (live[0] > 1)           { live[0] = 15; }
  for (i = 0; i < ncols; i++) { if (live[i] != 15) { live[i] = 0; } else { newnlive++; live[i] = 1; } } p = newnlive;
  writef("               %d\n", nReal); // right
  writef("            %9d", newnlive);  // left
  if ((newnlive < nlive) && (newnlive > 0)) return (true);
  if (!newnlive)
    writef("\n\n\n                  ***  D-reducible  ***\n\n");
  else
    writef("\n\n\n                ***  Not D-reducible  ***\n");
  return (false);
}

/* This generates all balanced signed matchings, and for each one, tests
 * whether all associated colourings belong to "live". It writes the answers
 * in the bits of the byteacters of "real". */
long dReduceTestMatch(long ring, ref byte[] real2, long[] power, ref byte[] live, long nbyte) pure {
  long a, b, n, nReal, realterm;
  long[10] interval;
  long[4][10] weight;
  long[4][16][16] matchweight;
  byte bit;

  /* "nReal" will be the number of balanced signed matchings such that all associated colourings belong to "live",
   * ie the total number of nonzero bits in the entries of "real" */
  nReal = 0; bit = 1; realterm = 0;
  /* First, it generates the matchings not incident with the last ring edge */

  for (a = 2; a <= ring; a++) {
    for (b = 1; b < a; b++) {
      matchweight[a][b][0] = 2 * (power[a] + power[b]);
      matchweight[a][b][1] = 2 * (power[a] - power[b]);
      matchweight[a][b][2] =      power[a] + power[b];
      matchweight[a][b][3] =      power[a] - power[b];
    }
  }
  for (a = 2; a < ring; a++) {
    for (b = 1; b < a; b++) { n = 0; weight[1] = matchweight[a][b];
      if (b >= 3)     { n = 1; interval[1] = 1;             interval[2]     = b - 1; }
      if (a >= b + 3) { n++;   interval[2 * n - 1] = b + 1; interval[2 * n] = a - 1; }
      augment(n, interval, 1L, weight, matchweight, live, real2, nReal, ring, 0L, 0L, bit, realterm, nbyte);
    }
  }

  /* now, the matchings using an edge incident with "ring" */
  for (a = 2; a <= ring; a++) {
    for (b = 1; b < a; b++) {
      matchweight[a][b][0] =  power[a] +     power[b];
      matchweight[a][b][1] =  power[a] -     power[b];
      matchweight[a][b][2] = -power[a] -     power[b];
      matchweight[a][b][3] = -power[a] - 2 * power[b];
    }
  }
  for (a = 2; a <= 2; a++) {
    for (b = 1; b < ring; b++) { n = 0; weight[1] = matchweight[ring][b];
      if (b >= 3)        { n = 1; interval[1] = 1;             interval[2]     = b - 1;    }
      if (ring >= b + 3) { n++;   interval[2 * n - 1] = b + 1; interval[2 * n] = ring - 1; }
      augment(n, interval, 1L, weight, matchweight, live, real2, nReal, ring,
        (power[ring + 1] - 1) / 2, 1L, bit, realterm, nbyte);
    }
  }

  return (nReal);
}

/* Finds all matchings such that every match is from one of the given
 * intervals. (The intervals should be disjoint, and ordered with smallest
 * first, and lower end given first.) For each such matching it examines all
 * signings of it, and adjusts the corresponding entries in "real" and
 * "live". */
void augment(long n, long[] interval, long depth, ref long[4][10] weight, long[4][16][16] matchweight, ref byte[] live,
  ref byte[] real2, ref long pnReal, long ring, long basecol, long on,
  ref byte pbit, ref long prealterm, long nbyte) pure {
    long h, i, j, r, newn, lower, upper;
    long[10] newinterval;

  checkreality(depth, weight, live, real2, pnReal, ring, basecol, on, pbit, prealterm, nbyte);
  depth++;
  for (r = 1; r <= n; r++) { lower = interval[2 * r - 1]; upper = interval[2 * r];
    for (i = lower + 1; i <= upper; i++) {
      for (j = lower; j < i; j++) { weight[depth] = matchweight[i][j];
        for (h = 1; h < 2 * r - 1; h++) { newinterval[h] = interval[h];                               } newn = r - 1;
        if  (j > lower + 1)             { newn++; newinterval[h++] = lower; newinterval[h++] = j - 1; }
        if  (i > j + 1)                 { newn++; newinterval[h++] = j + 1; newinterval[h++] = i - 1; }
        augment(newn, newinterval, depth, weight, matchweight, live, real2, pnReal, ring,
          basecol, on, pbit, prealterm, nbyte);
      }
    }
  }
}

/* For a given matching M, it runs through all signings, and checks which of
 * them have the property that all associated colourings belong to "live". It
 * writes the answers into bits of "real", starting at the point specified by
 * "bit" and "realterm". "basecol" is for convenience in computing the
 * associated colourings; it is zero for matchings not incident with "ring".
 * "on" is nonzero iff the matching is incident with "ring". */
void checkreality(long depth, long[4][] weight, ref byte[] live, ref byte[] real2, ref long pnReal, long ring,
  long basecol, long on, ref byte pbit, ref long prealterm, long nbyte) pure {
    long i, k, nbits = 1 << (depth - 1), col, parity;
    long[8] choice;
    ulong left;

  /* k will run through all subsets of M minus the first match */
  for (k = 0; k < nbits; k++, pbit = cast(byte)(cast(int)pbit << 1)) {
    if (!pbit) {
      pbit = 1; ++prealterm;
      assert(prealterm <= nbyte, "More than entries in real are needed");
    }
    if  (0 == (pbit & real2[prealterm])) { continue; } col = basecol; parity = ring & 1;
    for (i = 1, left = k; i < depth; i++, left = left >> 1) {
      /* i.e. if a_i=1, where k=a_1+2a_2+4a_3+... */
      if (0 != (left & 1)) {	parity = parity ^ 1; choice[i] = weight[i][1]; col += weight[i][3]; }
      else                 {                       choice[i] = weight[i][0]; col += weight[i][2]; }
    }
    if (parity) { choice[depth] = weight[depth][1]; col += weight[depth][3]; }
    else        { choice[depth] = weight[depth][0]; col += weight[depth][2]; }
    if (!isStillReal(col, choice, depth, live, on)) { real2[prealterm] = real2[prealterm] ^ pbit; }
    else                                            { pnReal++; }
  }
}

/* Given a signed matching, this checks if all associated colourings are in
 * "live", and, if so, records that fact on the bits of the corresponding
 * entries of "live". */
bool isStillReal(long col, long[] choice, long depth, ref byte[] live, long on) pure {
  long mark, i, j, twopower, b, c, ntwisted, nuntwisted;
  long[64] sum, twisted, untwisted;

  if (col < 0) { if (0 == live[-col]) return (false); twisted[ntwisted++] = -col; sum[0] = col; }
  else         { if (0 == live[col])  return (false); untwisted[nuntwisted++] = sum[0] = col; }
  for (i = 2, twopower = 1, mark = 1; i <= depth; i++, twopower = twopower << 1) { c = choice[i];
    for (j = 0; j < twopower; j++, mark++) { b = sum[j] - c;
      if (b < 0) { if (0 == live[-b]) return (false); twisted[ntwisted++] = -b; sum[mark] = b; }
      else       { if (0 == live[b])  return (false); untwisted[nuntwisted++] = sum[mark] = b; }
    }
  }
  /* Now we know that every coloring that theta-fits M has its code in
   * "live". We mark the corresponding entry of "live" by theta, that is,
   * set its second, third or fourth bit to 1 */
  if (on) {
    for (i = 0; i < ntwisted;   i++) live[twisted[i]]   |= 8;
    for (i = 0; i < nuntwisted; i++) live[untwisted[i]] |= 4;
  } else {
    for (i = 0; i < ntwisted;   i++) live[twisted[i]]   |= 2;
    for (i = 0; i < nuntwisted; i++) live[untwisted[i]] |= 2;
  }

  return (true);
}


// --------------------------------------------------------------------------------------------------------------------
/* checks that no colouring in live is the restriction to E(R) of a
 * tri-coloring of the free extension modulo the specified contract */
void checkCReduce(byte[] live, long nlive, tp_angle diffangle, tp_angle sameangle, long[] contract, long[] power,
  long[string] gConfMajor) {
    long j, i, u;
    long[] dm, sm;
    long[EDGES] c, forbidden;	/* called F in the notes */
    long ring = gConfMajor["ring"], bigno = gConfMajor["bigno"];
    long start = gConfMajor["edges"];	/* called s in the notes */

  if (0 == nlive) {
    if (0 == contract[0]) { writef("\n"); return; }
    else                  { assert(0, "         ***  ERROR23: CONTRACT PROPOSED  ***"); }
  }
  assert((0 != contract[0]),       "       ***  ERROR24: NO CONTRACT PROPOSED  ***");
  assert(nlive == contract[EDGES], "       ***  ERROR25: DISCREPANCY IN EXTERIOR SIZE  ***");
  while (contract[start]) { start--; } c[start] = 1; j = start;
  while (contract[--j])   {          } dm = diffangle[j]; sm = sameangle[j]; c[j] = 1;
  for (u = 4, i = 1; i <= dm[0]; i++) { u |=  c[dm[i]]; }
  for (i = 1; i <= sm[0]; i++)        { u |= ~c[sm[i]]; } forbidden[j] = u;

  for (;;) {
    while (forbidden[j] & c[j]) {
      if (isEndCCR(j, c, contract, start)) return;
    }
    if (j == 1) {
      assert(!isInLive(c, power, ring, live, bigno), "       ***  ERROR26: INPUT CONTRACT IS INCORRECT  ***");
      if (isEndCCR(j, c, contract, start)) return;
      continue;
    }
    while (contract[--j])                 {} dm = diffangle[j]; sm = sameangle[j]; c[j] = 1;
    for   (u = 0, i = 1; i <= dm[0]; i++) { u |=  c[dm[i]]; }
    for   (i = 1;        i <= sm[0]; i++) { u |= ~c[sm[i]]; } forbidden[j] = u;
  }
}

bool isEndCCR(ref long j, ref long[EDGES] c, long[] contract, long start) {
  c[j] <<= 1;
  while (c[j] & 8) {
    while (contract[++j]) {}
    if (j >= start) {
      writef("               ***  Contract confirmed  ***\n\n");
      return (true);
    }
    c[j] <<= 1;
  }
  return (false);
}

// Same as "chgLive" above, except now it returns whether the colouring is in live, and does not change live.
bool isInLive(long[] col, long[] power, long ring, byte[] live, long bigno) pure {
  long colno, i, min, max, w;
  long[5] weight;

  for (i = 1; i <= ring; i++) { weight[col[i]] += power[i]; } min = max = weight[4];
  for (i = 1; i <= 2;    i++) {
    w = weight[i];
    if      (w < min) min = w;
    else if (w > max) max = w;
  }
  colno = bigno - 2 * min - max;
  return (0 != live[colno]);
}


// --------------------------------------------------------------------------------------------------------------------
// これは四色定理の可約性を調べるプログラムです
// 0
//    This has ring-size 6, so there are 122 colourings total,
//    and 95 balanced signed matchings.
//    There are 16 colourings that extend to the configuration.
//             remaining               remaining balanced
//            colourings               signed matchings
//               106                       25
//               10                       15
//               8                       11
//               4                       3
//               1                       0
//               0
//                   ***  D-reducible  ***
// プログラムは正常終了しました

// 0x00 … 0 … 0
// 0x01 … +1 … 1
// ：
// 0x7e … +126 … 126
// 0x7f … +127 … 127
// 0x80 … -128 … 128
// 0x81 … -127 … 129
// 0x82 … -126 … 130
// ：
// 0xfe … -2 … 254
// 0xff … -1 … 255



