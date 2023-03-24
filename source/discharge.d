module discharge;

import std.file   : readText;
import std.json   : JSONValue, parseJSON;
import std.stdio  : writeln, writef, writefln, File;
import std.string : format, chomp, split;
import std.conv   : to;


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
void apply() {
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



