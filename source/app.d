// $ cd 4CT_D
// reduce or discharge をコメントアウトで切り替えて下さい
// $ dub もしくはデバッグ時に dub build --arch=x86_64 --build=debug --compiler=ldc2
import reduce    : reduce;
import discharge : discharge;


// --------------------------------------------------------------------------------------------------------------------
void main()
{
  // reduce;
  discharge(7);
}



