module discharge;

import std.file   : readText;
import std.json   : JSONValue, parseJSON;
import std.stdio  : writeln, writef, writefln, File;
import std.string : format, chomp, split;
import std.conv   : to;

const int VERTS      = 27;           /* max number of vertices in a free completion + 1 */
const int DEG        = 13;           /* max degree of a vertex in a free completion + 1, must be at least 13 because of row 0 */
const int CONFS      = 640;          /* max number of configurations		*/
const int MAXVAL     = 12;
const int CARTVERT   = (5*MAXVAL+2); /* domain of l_A, u_A, where A is an axle */
const int INFTY      = 12;           /* the "12" in the definition of limited part  */
const int MAXOUTLETS = 110;          /* max number of outlets */
const int MAXSTR     = 256;          /* max length of an input string */
const int MAXSYM     = 50;           /* max number of symmetries */
const int MAXELIST   = 134;          /* length of edgelist[a][b] */
const int MAXASTACK  = 5;            /* max height of Astack (see "reduce") */
const int MAXLEV     = 12;           /* max level of an input line + 1 */

alias tp_vertices = int[CARTVERT];
alias tp_adjmat   = int[CARTVERT][CARTVERT];
alias tp_confmat  = long[DEG][VERTS];
struct Tp_axle {
  tp_vertices low;
  tp_vertices upp;
}
struct Tp_outlet {
  int number;  /* +/-n for outlet corresponding to rule n, always !=0 */
  int nolines; /* |M(T)| */
  int value;
  int[17] pos;
  int[17] low;
  int[17] upp;
}
struct Tp_cond {
  int n;
  int m;
}
struct Tp_posout {
  Tp_outlet T;
  int x;
}
struct Tp_query {
  int u;
  int v;
  int z;
  int xi;
}
alias tp_question = Tp_query[VERTS];
alias tp_edgelist = int[MAXELIST][9][12];


// --------------------------------------------------------------------------------------------------------------------
void discharge(int deg) {
  writefln("hoge. %d", deg);

  string fileName1 = "data/d_good_confs.json";
  auto jv1 = readText(fileName1).parseJSON();
  JSONValue[] jarr1 = jv1.array();
  writeln(jarr1[0]["a"]);

  string fileName2 = "data/d_rules.json";
  auto jv2 = readText(fileName2).parseJSON();
  JSONValue[] jarr2 = jv2.array();
  writeln(jarr2[10]["z"]);

  string filename3 = "data/d_tactics07.txt";
	auto fin = File(filename3,"r");
  int cnt = 0;
	foreach (line; fin.byLine) {
    if (cnt == 10) break;
		writeln(line.chomp.split);
    ++cnt;
  }
  writeln("33".to!(int) + 1);

}
private: // ここ以下すべて private


// --------------------------------------------------------------------------------------------------------------------
/*************************************************************************
      apply   Verifies symmetry line as described in [D]
*************************************************************************/
void apply(int k, int epsilon, int level, int line, Tp_axle A, Tp_outlet[] sym, int nosym, int lineno) pure {
  int i;

  // if  (sscanf(S, "%*s%d%d%d%d", &k, &epsilon, &level, &line) != 4) error("Syntax error", lineno);
  assert((k >= 0 && k <= A.low[0] && epsilon >= 0 && epsilon <= 1), "Illegal symmetry");
  for (i = 0; i < nosym; i++) { if (sym[i].number == line) break; }
  assert((i < nosym),                                               "No symmetry as requested");
  assert((sym[i].nolines == level + 1),                             "Level mismatch");
  if  (epsilon == 0) {
    assert((0 != outletForced(A, sym[i], k + 1)),                   "Invalid symmetry");
  } else {
    assert((0 != reflForced(A, sym[i], k + 1)),                     "Invalid reflected symmetry");
  }
}

/*********************************************************************
        outletForced
If (T,x) is enforced by A, then returns the value of T, otherwise 0
*********************************************************************/
int outletForced(Tp_axle A, Tp_outlet T, int x) pure {
  int i, p, deg = A.low[0], xx = x - 1;

  for (i = 0; i < T.nolines; ++i) {
    p = T.pos[i];
    p = xx + (p - 1) % deg < deg ? p + xx : p + xx - deg;
    if (T.low[i] > A.low[p] || T.upp[i] < A.upp[p]) return (0);
  }
  return (T.value);
}

/*********************************************************************
        outletPermitted
If (T,x) is permitted by A, then returns the value of T, otherwise 0
*********************************************************************/
int outletPermitted(Tp_axle A, Tp_outlet T, int x) pure {
  int deg = A.low[0], i, p, xx = x - 1;

  for (i = 0; i < T.nolines; ++i) {
    p = T.pos[i]; p = xx + (p - 1) % deg < deg ? p + xx : p + xx - deg;
    if (T.low[i] > A.upp[p] || T.upp[i] < A.low[p]) return (0);
  }
  return (T.value);
}

/************************************************************************
        reflForced
Returns the value of T if M is fan-free and every cartwheel compatible
with A is compatible with tau^(x-1)sigma M, where M is the axle
corresponding to T
************************************************************************/
int reflForced(Tp_axle A, Tp_outlet T, int x) pure {
  int deg = A.low[0], i, p, q, xx = x - 1;

  for (i = 0; i < T.nolines; ++i) {
    p = T.pos[i];
    p = xx + (p - 1) % deg < deg ? p + xx : p + xx - deg;
    if      (p < 1 || p > 2 * deg) return (0);
    else if (p <= deg)             q = deg - p + 1;
    else if (p < 2 * deg)          q = 3 * deg - p;
    else                           q = 2 * deg;
    if (T.low[i] > A.low[q] || T.upp[i] < A.upp[q]) return (0);
  }
  return (T.value);
}


// --------------------------------------------------------------------------------------------------------------------
void caseSplit() {
}


// --------------------------------------------------------------------------------------------------------------------
void libDischarge() {
}


// --------------------------------------------------------------------------------------------------------------------
void reduce() {
}



