module discharge;

import std.file   : readText;
import std.json   : JSONValue, parseJSON;
import std.stdio  : writeln, writef, writefln, File;
import std.string : format, chomp, split;
import std.conv   : to;
import std.range  : iota;

const int VERTS      = 27;           /* max number of vertices in a free completion + 1 */
const int DEG        = 13;           /* max degree of a vertex in a free completion + 1,
                                      * must be at least 13 because of row 0 */
const int CONFS      = 640;          /* max number of configurations */
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
struct Tp_xyv {
  int x;
  int y;
  int v;
}
alias tp_question = Tp_query[VERTS];
alias tp_edgelist = int[MAXELIST][9][12];


// --------------------------------------------------------------------------------------------------------------------
void discharge(int deg) {
  //string[50] ch;
  int lev;	/* level of line being processed */
  int prtline;	/* print details about line number "prtline" */
  int nosym;	/* number of symmetries, see "sym" below */
  //char str[MAXSTR];	/* holds input line */
  //char fname[MAXSTR];	/* name of file to be tested */
  Tp_axle[MAXLEV + 1] axles;	/* axles[l] is A_l of [D] */
  //Tp_outlet[2 * MAXOUTLETS] sym;	/* sym[i] (i=0,..,nosym) are T_i (i=0,..,t-1) of [D] */
  int a, printmode, print = 4;
  //char *ch;
  Tp_xyv[10] xyv;

  writefln("deg = %d", deg);

  string fileName1 = "data/d_good_confs.json";
  auto jv1 = readText(fileName1).parseJSON();
  JSONValue[] jarr1 = jv1.array();
  writeln(jarr1[0]["a"]);

  string fileName2 = "data/d_rules.json";
  auto jv2 = readText(fileName2).parseJSON();
  JSONValue[] jarr2 = jv2.array();
  writeln(jarr2[10]["z"]);

  Tp_outlet[206] sym;

  writef("\n\n\n");

  // axles の初期化
  axles[0].low[0] = deg;
  foreach (n; iota(5 * deg)) { axles[0].low[n + 1] = 5; }
  axles[0].upp[0] = deg;
  foreach (n; iota(5 * deg)) { axles[0].upp[n + 1] = INFTY; }

  string filename3 = "data/d_tactics07.txt";
  auto fin = File(filename3,"r");
  int cnt = 0;
  foreach (line; fin.byLine) {
    auto ch = line.chomp.split;
    writef("%d, nosym=%d  ", cnt, nosym); writeln(ch);
    if (ch[0] == "Degree") { ++cnt; continue; }
    if (ch[0] == "Q.E.D.") { break; }
    //if (cnt == 20) break;

    switch (ch[1]) {
    case "S":
      apply(ch[2].to!(int), ch[3].to!(int), ch[4].to!(int), ch[5].to!(int), axles[lev], sym, nosym, cnt + 1);
      nosym = delSym(nosym, sym, lev);
      lev--;
      break;
    case "R":
      assert(reduce(axles[lev], cnt + 1, 1), "Reducibility failed");
      nosym = delSym(nosym, sym, lev);
      lev--;
      break;
    case "H":
      foreach (i, temp; ch[2 .. $]) {
        auto temp2 = split(temp[1 .. $ - 1], ",");
        xyv[i].x = temp2[0].to!(int);
        xyv[i].y = temp2[1].to!(int);
        xyv[i].v = temp2[2].to!(int);
      }
      libDischarge(xyv, axles[lev], cnt + 1, print, sym);
      nosym = delSym(nosym, sym, lev);
      lev--;
      break;
    case "C":
      caseSplit(ch[2].to!(int), ch[3].to!(int), axles[lev], axles[lev + 1], sym, nosym, lev, cnt + 1, print, deg);
      lev++;
      break;
    default:
      assert(0, "Invalid instruction");
    }

    cnt++;
  }
  writefln("end.");

}
private: // ここ以下すべて private


// --------------------------------------------------------------------------------------------------------------------
/*************************************************************************
      apply   Verifies symmetry line as described in [D]
*************************************************************************/
void apply(int k, int epsilon, int level, int line, Tp_axle A, Tp_outlet[] sym, int nosym, int lineno) {
  int i;

  // if  (sscanf(S, "%*s%d%d%d%d", &k, &epsilon, &level, &line) != 4) error("Syntax error", lineno);
  assert((k >= 0 && k <= A.low[0] && epsilon >= 0 && epsilon <= 1), "Illegal symmetry");
  writefln("    %d %d %d %d", sym[0].number, sym[1].number, sym[2].number, sym[3].number);
  for (i = 0; i < nosym; i++) {
    if (sym[i].number == line) {
      break;
    }
  }
  assert((i < nosym),                                               "No symmetry as requested");
  assert((sym[i].nolines == level + 1),                             "Level mismatch");
  if  (epsilon == 0) {
    assert((0 != outletForced(A, sym[i], k + 1)),                   "Invalid symmetry");
  } else {
    assert((0 != reflForced(  A, sym[i], k + 1)),                   "Invalid reflected symmetry");
  }
}

int delSym(int nosym, Tp_outlet[] sym, int lev) {
  if (nosym < 1 || sym[nosym - 1].nolines - 1 < lev) return nosym;
  else                                               return delSym((nosym - 1), sym, lev);
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
/*************************************************************************
      caseSplit  Verifies condition line as described in [D]
*************************************************************************/
void caseSplit(int n, int m, ref Tp_axle A, ref Tp_axle A2, ref Tp_outlet[206] sym,
  ref int pnosym, int lev, int lineno, int print, int deg) {
    int i, j, good;
    static Tp_cond[MAXLEV] cond;

  //if (sscanf(S, "%*s%d%d", &n, &m) != 2)    error("Syntax error", lineno);
  /* check condition and compatibility with A */
  assert((n >= 1 && n <= 5 * deg),        "Invalid vertex in condition");
  assert((m >= -8 && m <= 9 && (m <= -5 || m >= 6)), "Invalid condition");
  j = (n - 1) / deg;
  i = (n - 1) % deg + 1;
  assert(((n <= 2 * deg) || (A.low[i] == A.upp[i] && A.low[i] >= j + 4)),
                                            "Condition not compatible with A");
  A2 = A; //copyAxle(A + 1, A);
  if (m > 0) { /* new lower bound */
    assert((A.low[n] < m && m <= A.upp[n]),  "Invalid lower bound in condition");
    A      .upp[n] = m - 1;
    A2.low[n] = m;
  } else {   /* new upper bound */
    assert((A.low[n] <= -m && -m < A.upp[n]), "Invalid upper bound in condition");
    A      .low[n] = 1 - m;
    A2.upp[n] = -m;
  }

  /* remember symmetry unless contains a fan vertex */
  for (i = 0, good = 1; i <= lev; i++)
    if (cond[i].n > 2 * deg || cond[i].n < 1) good = 0;
  if (good) { /* remember symmetry */
    assert((pnosym < MAXSYM), "Too many symmetries");
    if (print >= 0)   writef("Adding symmetry:");
    sym[pnosym].number = lineno;
    sym[pnosym].value = 1;
    sym[pnosym].nolines = lev + 1;
    for (i = 0; i <= lev; ++i) {
      sym[pnosym].pos[i] = cond[i].n;
      if (cond[i].m > 0) {
        sym[pnosym].low[i] = cond[i].m;
        sym[pnosym].upp[i] = INFTY;
      } else {
        sym[pnosym].low[i] = 5;
        sym[pnosym].upp[i] = -cond[i].m;
      }
      if (print >= 0) writef(" (%d,%d,%d)", sym[pnosym].pos[i], sym[pnosym].low[i], sym[pnosym].upp[i]);
    }
    if (print >= 0) { writef("\n"); }
    pnosym++;
  } else if (print >= 0) {
    writef("Symmetry not added\n");
  }
  cond[lev].n = n; cond[lev].m = m; cond[lev + 1].n = 0; cond[lev + 1].m = 0;
}


// --------------------------------------------------------------------------------------------------------------------
/*************************************************************************
    libDischarge
If str==NULL it assumes that A is the trivial axle of degree deg, where
deg=A.low[0]. It reads rules and computes the corresponding outlets. It
writes the outlets into the file specified by OUTLETFILE so they can be
verified for accuracy.
If str!=NULL it verifies hubcap line as described in [D]
**************************************************************************/
void libDischarge(Tp_xyv[] xyv, Tp_axle A, int lineno, int print, Tp_outlet[] sym) {
  int i, j, a, total = 0, deg = A.low[0];
  int[2 * MAXOUTLETS + 1] s;
  static int nouts = -7;
  static Tp_posout[2 * MAXOUTLETS] posout;

  if (nouts == -7) {
    for (i = 0; i < 206; i++) posout[i].T = sym[i];
    nouts = 0;
  }

  for (i = 0; i < xyv.length; i++) {
    if  (xyv[i].x == 0) break;
    if  (print >= 3)       writef("\n-.Checking hubcap member (%d,%d,%d)\n", xyv[i].x, xyv[i].y, xyv[i].v);
    for (j = 0; j < nouts; j++) { posout[j].x = xyv[i].x; s[j] = 0; }
    if  (xyv[i].x != xyv[i].y)          { for (; j < 2 * nouts; j++) { posout[j].x = xyv[i].y; s[j] = 0; } }
    s[j] = 99; /* to indicate end of list */
    libDischargeCore(A, posout, s, xyv[i].v, 0, 0, lineno, print);
  }
  if (print >= 3) { writef("\n"); }
}

/*************************************************************************
    libDischargeCore    Verifies (H1)
*************************************************************************/
void libDischargeCore(Tp_axle A, Tp_posout[] posout, int[] s, int maxch, int pos, int depth, int lineno, int print) {
  int deg = A.low[0], i, x, forcedch, allowedch;
  bool good;
  int[] sprime;
  Tp_axle AA;

  /* 1. compute forced and permitted rules, allowedch, forcedch, update s */
  forcedch = allowedch = 0;
  for (i = 0; s[i] < 99; i++) {
    if      (s[i] > 0)                          forcedch += posout[i].T.value;
    if      (s[i])                              continue;
    if      (outletForced(A, posout[i].T, posout[i].x))     { s[i] = 1; forcedch += posout[i].T.value; }
    else if (!outletPermitted(A, posout[i].T, posout[i].x)) s[i] = -1;
    else if (posout[i].T.value > 0)                  allowedch += posout[i].T.value;
  }

  /* 2. print */
  if (print >= 3) {
    writef("POs: ");
    for (i = 0; s[i] < 99; i++) {
      if (s[i] < 0)  continue;
      if (s[i] == 0) writef("?");
      writef("%d,%d ", posout[i].T.number, posout[i].x);
    }
    writef("\n");
  }

  /* 3. check if inequality holds */
  if (forcedch + allowedch <= maxch) {
    if (print >= 3) writef("Inequality holds. Case done.\n");
    return; // true end 1
  }

  /* 4. check reducibility */
  if (forcedch > maxch) {
    assert(reduce(A, lineno, print >= 4 ? 1 : 0), "Incorrect hubcap upper bound");
    if (print >= 3 && print < 4)               writef("Reducible. Case done.\n");
    return; // true end 2
  }

  /* 5. */
  for (; s[pos] < 99; pos++, good = true) {
    if (s[pos] || posout[pos].T.value < 0) continue;

    /* accepting positioned outlet PO, computing AA */
    AA = A;
    assert321(pos, posout[pos].x, AA, posout, deg);

    /* recursion with PO forced */
    if (isGood(pos, posout[pos].x, AA, posout, s, print)) {
      sprime = s;
      sprime[pos] = 1;
      if (print >= 3) { writef("Starting recursion with "); writef(",%d forced\n", x); }
      libDischargeCore(AA, posout, sprime, maxch, pos + 1, depth + 1, lineno, print);
    }

    /* rejecting positioned outlet PO */
    if (print >= 3) { writef("Rejecting positioned outlet "); writef(",%d. ", x); }
    s[pos] = -1; allowedch -= posout[pos].T.value;
    if (allowedch + forcedch <= maxch) {
      if (print >= 3) writef("Inequality holds.\n");
      return; // true end 3
    }
    else if (print >= 3) writef("\n");
  }

  /* 6. error */
  assert(0, "Unexpected error 101");
}

void assert321(int pos, int x, Tp_axle AA, Tp_posout[] posout, int deg) {
  int p, i;
  for (i = 0; i < posout[pos].T.nolines; ++i) {
    p = posout[pos].T.pos[i];
    p = x - 1 + (p - 1) % deg < deg ? p + x - 1 : p + x - 1 - deg;
    if (posout[pos].T.low[i] > AA.low[p]) AA.low[p] = posout[pos].T.low[i];
    if (posout[pos].T.upp[i] < AA.upp[p]) AA.upp[p] = posout[pos].T.upp[i];
    assert((AA.low[p] <= AA.upp[p]), "Unexpected error 321");
  }
}

/* Check if a previously rejected positioned outlet is forced to apply */
bool isGood(int pos, int x, Tp_axle AA, Tp_posout[] posout, int[] s, int print) {
  int i;
  for (i = 0; i < pos; i++) {
    if (s[i] == -1 && outletForced(AA, posout[i].T, posout[i].x)) {
      if (print >= 3) {
        writef("Positioned outlet ");
        writef(",%d can't be forced, because it forces %d,%d\n", x, posout[i].T.number, posout[i].x);
      }
      return false;
    }
  }
  return true;
}

// --------------------------------------------------------------------------------------------------------------------
bool reduce(Tp_axle A, int lineno, int print) {
  return (true);
}



