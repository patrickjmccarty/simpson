/* Automatically generated. Do not edit ! */
#include "tclutil.h"
/* This file has been generated from the following Tcl source file(s):
 *    simpson.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char simpson_tcl[] = 
    "\n"
    "\n"
    "proc markvars {} {\n"
    "      uplevel #0 {\n"
    "        set _simpson_omitvars(_simpson_omitvars) 1\n"
    "        set _simpson_omitvars(_slave_state) 1\n"
    "        set _simpson_omitvars(ssval) 1\n"
    "        set _simpson_omitvars(spinsys) 1\n"
    "#        set _simpson_omitvars(spinsysres) 1\n"
    "        set _simpson_omitvars(_globs) 1\n"
    "        set _simpson_omitvars(__name) 1\n"
    "        foreach __name [info vars] {\n"
    "           set _simpson_omitvars($__name) 1\n"
    "        }\n"
    "      }\n"
    "} \n"
    "    \n"
    "proc savestate {} {\n"
    "      uplevel #0 {\n"
    "      if [info exists _slave_state] {unset _slave_state}\n"
    "      set _slave_state {}\n"
    "      set _globs \"global\"\n"
    "        foreach __name [info globals]  {\n"
    "          if [info exists _simpson_omitvars($__name)] continue \n"
    "          lappend _globs $__name\n"
    "          if [array exists $__name] {\n"
    "             lappend _slave_state  [list array set $__name [array get $__name]]\n"
    "          } else {\n"
    "               lappend _slave_state  [list set $__name [set $__name]]\n"
    "          }\n"
    "        }\n"
    "        return [join [list [join $_globs \" \"] [join $_slave_state \\n] ] \\n]\n"
    "      }\n"
    "}\n"
    "\n"
    "proc signalhandler {code message} {\n"
    "  global stop sigwascalled par\n"
    "\n"
    "  if [info exists sigwascalled] {\n"
    "     puts stderr \"Program is aborted.\"\n"
    "     exit\n"
    "  }\n"
    "  set stop 1\n"
    "  puts stderr \"Received '$message', setting flag to terminate.\"\n"
    "  puts stderr \"Repeat action to abort program immediately.\"\n"
    "  set sigwascalled 1\n"
    "}\n"
    "\n"
    "\n"
    "proc setvar {varname data {overwrite 0}} {\n"
    "  global $varname\n"
    "  set var $varname\n"
    "\n"
    "  set data [split $data \"\\n\"]\n"
    "\n"
    "  foreach lst $data {\n"
    "    if [string match #* [string trimleft $lst]] continue\n"
    "    set nam [lindex $lst 0]\n"
    "    if ![string length $nam] continue\n"
    "    if !$overwrite {\n"
    "      if [info exists ${var}($nam)] {\n"
    "        puts stderr \"error: '$nam' already exists in array ${var}\"\n"
    "        exit\n"
    "      }\n"
    "    }\n"
    "    set ${var}($nam) [lrange $lst 1 9999]\n"
    "  }\n"
    "}\n"
    "\n"
    "\n"
    "setvar ssnam {\n"
    "  shift N iso aniso eta alpha beta gamma\n"
    "  jcoupling N N iso aniso eta alpha beta gamma\n"
    "  quadrupole N order aniso eta alpha beta gamma\n"
    "  dipole N N aniso alpha beta gamma\n"
    "  dipole_ave N N aniso eta alpha beta gamma\n"
    "  mixing N N\n"
    "}\n"
    "\n"
    "\n"
    "proc ssSubstExpr {ex resfreq} {\n"
    "  global ssval par\n"
    "\n"
    "  if [regexp {^[0-9.Ee+-]+$} $ex] {\n"
    "    return $ex\n"
    "  }\n"
    "  set origex $ex\n"
    "  if [regsub {([0-9.Ee+-]+)p} $ex \"(\\\\1*($resfreq))\" dummy] {\n"
    "     if {$resfreq == 0} {\n"
    "        puts stderr \"error in expression '$origex'. Can only use ppm to hz convertion for chemical shift.\"\n"
    "        exit\n"
    "     }\n"
    "     if ![info exists par(proton_frequency)] {\n"
    "        puts stderr \"error in expression '$origex'. 'proton_frequency' must be set when converting ppm to Hz values.\"\n"
    "        exit\n"
    "     }\n"
    "     set freq [expr abs(($par(proton_frequency)/1.0e6)*($resfreq/1.0e6))]\n"
    "     if ![regsub {([0-9.Ee+-]+)p} $ex \"(\\\\1*($freq))\" ex] {\n"
    "        puts stderr \"error: illegal value of 'proton_frequency' = $par(proton_frequency) \"\n"
    "        exit\n"
    "     }\n"
    "  }\n"
    "  set i 0\n"
    "  while {[regsub -all {([a-z]+(_[0-9]+)+_[a-z]+)} $ex {$ssval(\\1)} ex] != 0} {\n"
    "    set ex [subst $ex]\n"
    "    if {[incr i] > 100} {\n"
    "     puts stderr \"error: substitution over 100 times af the expression '$origex'\"\n"
    "     puts stderr \"       have you made any circular references ?\"\n"
    "     exit\n"
    "    }\n"
    "  }    \n"
    "  if [catch {set ex [expr $ex]} res] {\n"
    "    puts stderr \"error: $res\"\n"
    "    exit  \n"
    "  }\n"
    "  return $ex\n"
    "}\n"
    "\n"
    "\n"
    "proc ssSetValues {} {\n"
    "  global spinsys ssval ssnam\n"
    "\n"
    "  foreach iact [array names spinsys] {\n"
    "    set val $spinsys($iact)\n"
    "    if {[string compare $iact \"nuclei\"] == 0} {\n"
    "      set ssval(nuclei) $val\n"
    "      continue\n"
    "    }\n"
    "    if {[string compare $iact \"channels\"] == 0} {\n"
    "      set ssval(channels) $val\n"
    "      continue\n"
    "    }\n"
    "    regsub {[X]+$} $iact {} iact\n"
    "\n"
    "    if ![info exists ssnam($iact)] {\n"
    "      puts stderr \"error: in spinsys section at '$iact', correct syntax is:\\n  spinsys {\"\n"
    "      foreach i [array names ssnam] {\n"
    "        puts stderr \"    $i $ssnam($i)\"\n"
    "      }\n"
    "      puts stderr \"  }\"\n"
    "      exit\n"
    "    }\n"
    "    set nam $ssnam($iact)\n"
    "\n"
    "    if {[llength $val] != [llength $nam]} {\n"
    "       puts stderr \"error: in spinsys: input list    '$iact $val'\"\n"
    "       puts stderr \"       does not match definition '$iact $nam'\"\n"
    "       exit\n"
    "    }\n"
    "    set nuc \"\"\n"
    "    set k 1\n"
    "    foreach j $nam {\n"
    "      if {[string compare $j \"N\"] == 0} {\n"
    "        set nuc \"${nuc}_[lindex $val [expr $k - 1]]\"\n"
    "      } else {\n"
    "        set ssval($iact${nuc}_$j) [lindex $val [expr $k - 1]]\n"
    "      }\n"
    "      incr k\n"
    "    }\n"
    "  }\n"
    "}\n"
    "\n"
    "\n"
    "proc ssSetSpinsys {} {\n"
    "  global spinsys spinsysres ssval ssnam\n"
    "\n"
    "  foreach iact [array names spinsys] {\n"
    "    set orignam $iact\n"
    "    set val $spinsys($orignam)\n"
    "\n"
    "    if {[string compare $iact \"nuclei\"] == 0} {\n"
    "      set spinsysres(nuclei) $ssval(nuclei)\n"
    "      continue\n"
    "    }\n"
    "    if {[string compare $iact \"channels\"] == 0} {\n"
    "      set spinsysres(channels) $ssval(channels)\n"
    "      continue\n"
    "    }\n"
    "    regsub {[X]+$} $iact {} iact\n"
    "    set nam $ssnam($iact)\n"
    "    if {[llength $val] != [llength $nam]} {\n"
    "       puts stderr \"error: input list '$val' for interaction '$iact' doesn't match definition '$nam'\"\n"
    "       exit\n"
    "    }\n"
    "    set lst {}\n"
    "    set nuc \"\"\n"
    "    set k 1\n"
    "    foreach j $nam {\n"
    "      set ival [lindex $val [expr $k - 1]]\n"
    "      if {[string compare $j \"N\"] == 0} {\n"
    "        set nuc \"${nuc}_$ival\"\n"
    "        lappend lst $ival\n"
    "      } else {\n"
    "        lappend lst $ssval($iact${nuc}_$j)\n"
    "      }\n"
    "      incr k\n"
    "    }\n"
    "    set spinsysres($orignam) $lst\n"
    "  }\n"
    "}\n"
    "\n"
    "\n"
    "\n"
    "proc spinsys_resolve { { fitval {} } } {\n"
    "  global spinsys ssval spinsysres par\n"
    "\n"
    "  ssSetValues\n"
    "  foreach i $fitval {\n"
    "    set nam [lindex $i 0]\n"
    "    set val [lindex $i 1]\n"
    "    if ![info exists ssval($nam)] {\n"
    "      puts stderr \"error: variable '$nam' does not correspond to any variable in section 'spinsys'\"\n"
    "      exit\n"
    "    }\n"
    "    set ssval($nam) $val      \n"
    "  }\n"
    "\n"
    "  ssSetSpinsys  \n"
    "  if ![info exists spinsysres(nuclei)] {\n"
    "    puts stderr \"error: spinsys has to include 'nuclei' definition\"\n"
    "    exit\n"
    "  }\n"
    "  set k 1\n"
    "  foreach j $spinsysres(nuclei) {\n"
    "    set nuc($k) $j\n"
    "    incr k\n"
    "  }\n"
    "  foreach i [array names spinsysres] {  \n"
    "    if ![string compare $i nuclei] {\n"
    "      continue\n"
    "    }\n"
    "    if ![string compare $i channels] {\n"
    "      continue\n"
    "    }\n"
    "    set lst {}\n"
    "    set resfreq 0\n"
    "    if ![string first shift $i ] {\n"
    "       set resfreq [resfreq $nuc([lindex $spinsysres($i) 0]) 1e6]\n"
    "    }\n"
    "    set k 1\n"
    "    foreach j $spinsysres($i) {\n"
    "      if {$k == 2 || $k == 3} {\n"
    "        lappend lst [ssSubstExpr $j $resfreq]\n"
    "      } else {\n"
    "        lappend lst [ssSubstExpr $j 0]\n"
    "      }\n"
    "      incr k\n"
    "    }\n"
    "    set spinsysres($i) $lst\n"
    "  }\n"
    "}\n"
    "\n"
    "proc spinsys { data } {\n"
    "  global spinsys\n"
    "\n"
    "  set okpar { channels nuclei dipole quadrupole shift jcoupling mixing dipole_ave }\n"
    "\n"
    "  set data [split $data \"\\n\"]\n"
    "  foreach lst $data {\n"
    "    if [string match #* [string trimleft $lst]] continue\n"
    "    set nam [lindex $lst 0]\n"
    "    if ![string length $nam] continue\n"
    "\n"
    "    if {[lsearch -exact $okpar $nam] == -1} {\n"
    "      puts stderr \"error: unknown name '$nam' in spinsys section, must be one\"\n"
    "      set f [join $okpar {, }]\n"
    "      puts stderr \"       of $f\"\n"
    "      exit\n"
    "    }\n"
    "    \n"
    "    while {[info exists spinsys($nam)]} {\n"
    "      set nam \"${nam}X\"\n"
    "    }\n"
    "    set spinsys($nam) [lrange $lst 1 9999]\n"
    "  }\n"
    "}\n"
    "\n"
    "\n"
    "\n"
    "proc cmp_length {a b} {\n"
    "  return [expr [string length $a] < [string length $b]]\n"
    "}\n"
    "\n"
    "proc par {data} {\n"
    "  global par\n"
    "\n"
    "  set okpar {\n"
    "    proton_frequency spin_rate sw sw1 np ni method rotor_angle\n"
    "    gamma_angles fixed_rep real_spec block_diag detect_operator\n"
    "    crystal_file start_operator name verbose various variable pulse_sequence\n"
    "    conjugate_fid dipole_check gamma_zero use_cluster new_cluster cluster_port\n"
    "    inner_rotor_angle outer_rotor_angle inner_spin_rate outer_spin_rate dor\n"
    "    string oc_tol_cg oc_tol_ls oc_mnbrak_step oc_max_iter oc_cutoff \n"
    "    oc_cutoff_iter oc_var_save_iter oc_var_save_proc oc_cg_min_step oc_max_brack_eval\n"
    "    oc_max_brent_eval oc_verbose oc_grad_level oc_method\n"
    "    rfprof_file use_3_angle_set acq_adjoint\n"
    "    zprofile zvals relax prop_method use_sparse num_cores averaging_file\n"
	"    points_per_cycle ED_symmetry\n"
	"    oc_lbfgs_eps oc_lbfgs_tol_ls oc_lbfgs_max_ls_eval oc_lbfgs_m\n"
	"    sparsity sparse_tol maxfulldim maxdimdiagonalize do_avg\n"
	"  }\n"
    "\n"
    "  set allowsubst {\n"
    "    proton_frequency spin_rate sw sw1 np ni rotor_angle\n"
    "    gamma_angles fixed_rep real_spec block_diag gamma_zero\n"
    "    inner_rotor_angle outer_rotor_angle inner_spin_rate outer_spin_rate\n"
    "    variable\n"
    "    oc_tol_cg oc_tol_ls oc_mnbrak_step oc_max_iter oc_cutoff \n"
    "    oc_cutoff_iter oc_var_save_iter oc_cg_min_step oc_max_brack_eval oc_max_brent_eval\n"
	"    points_per_cycle\n"
    "  }\n"
    "  set data [split $data \"\\n\"]\n"
    "  foreach lst $data {\n"
    "    if [string match #* [string trimleft $lst]] continue\n"
    "    set nam [lindex $lst 0]\n"
    "    set orignam $nam\n"
    "    if ![string length $nam] continue\n"
    "    if {[lsearch -exact $okpar $nam] == -1} {\n"
    "      puts stderr \"error: unknown name '$nam' in parameter section\"\n"
    "      exit\n"
    "    }\n"
    "    if ![string compare $nam variable] {\n"
    "      set nam [lindex $lst 1]\n"
    "      set lst [lrange $lst 1 9999]    \n"
    "    }\n"
    "    if ![string compare $nam string] {\n"
    "      set nam [lindex $lst 1]\n"
    "      if {![string length $nam]} {\n"
    "      	puts stderr \"error: incomplete definition of string variable\"\n"
    "	exit\n"
    "      }\n"
    "      set lst [join [lrange $lst 1 9999]]    \n"
    "    }\n"
    "    if [info exists par($nam)] {\n"
    "      puts stderr \"error: '$nam' already exists in array par\"\n"
    "      exit\n"
    "    }\n"
    "    set ex [lrange $lst 1 9999]\n"
    "    if {[lsearch -exact $allowsubst $orignam] != -1} {\n"
    "      set origex $ex\n"
    "      regsub -all {([A-Za-z_][A-Za-z0-9_]*\\(?)} $ex { \\1 } ex\n"
    "      regsub -all {([0-9\\.]+) ([Ee])} $ex {\\1\\2} ex\n"
    "      regsub -all {([0-9\\.]+[Ee]) ([+-][0-9\\.]+)} $ex {\\1\\2} ex\n"
    "      set ex \"$ex \"\n"
    "      foreach n [lsort -command cmp_length [array names par]] {\n"
    "        set i 0\n"
    "        if {[string first \" $n \" $ex] != -1} {        \n"
    "          if ![regsub -all \" \\($n\\) \" $ex { $par(\\1) } ex] {\n"
    "             puts stderr \"error: unable to match $n in $ex\"\n"
    "             exit\n"
    "          }\n"
    "          set ex [subst \"$ex\"]\n"
    "        }\n"
    "      }\n"
    "      if [catch {set ex [expr $ex]} res] {\n"
    "        puts stderr \"error: $res\"\n"
    "        exit  \n"
    "      }\n"
    "    }\n"
    "    set par($nam) $ex\n"
    "  }\n"
    "}\n"
    "\n"
    "proc fsimpson {{fitpar {}}} {\n"
    "  global par _fd\n"
    "\n"
    "  if [info exists par(verbose)] {\n"
    "    if {[string index $par(verbose) 0] == 1}  {\n"
    "      puts \"Parameters\"\n"
    "      foreach k [lsort [array names par]] {\n"
    "         puts [format \"  %-20s %s\" $k $par($k)]\n"
    "      }\n"
    "    }\n"
    "  }\n"
    "  spinsys_resolve $fitpar\n"
    "  set par(tcalc) [lindex [time {uplevel #0 {set _fd [internalsimpson]}}] 0]\n"
    "  set fd $_fd\n"
    "  unset _fd\n"
    "  return $fd\n"
    "}\n"
    "\n"
    "proc main {} {\n"
    "   puts \"There was no main section in the input file!\"\n"
    "}\n"
    "\n"
    "proc pulseq {} {\n"
    "   puts \"There was no pulseq section in the input file!\"\n"
    "}\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    "\n"
    ;
/* EOF */
/* This file has been generated from the following Tcl source file(s):
 *    ftools.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char ftools_tcl[] = 
    "\n"
    "proc fexpr {desc reexpr imexpr} {  \n"
    "   set reex [concat expr $reexpr]\n"
    "   set imex [concat expr $imexpr]\n"
    "   set np [fget $desc -np]\n"
    "   set ni [fget $desc -ni]\n"
    "   if {$ni > 1} {set np [expr $np*$ni]}\n"
    "   for {set i 1} {$i <= $np} {incr i} {\n"
    "      set c [findex $desc $i]\n"
    "      set re [lindex $c 0]\n"
    "      set im [lindex $c 1]\n"
    "      set kre [eval $reex]\n"
    "      set kim [eval $imex]\n"
    "      fsetindex $desc $i $kre $kim\n"
    "   }\n"
    "}\n"
    "\n"
    ;
/* EOF */
/* This file has been generated from the following Tcl source file(s):
 *    misc.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char misc_tcl[] = 
    "\n"
    "\n"
    "proc csapar {s1 s2 s3} {\n"
    "  set iso [expr ($s1+$s2+$s3)/3.0]\n"
    "\n"
    "  set v [list [expr abs($s1-$iso)] [expr abs($s2-$iso)] [expr abs($s3-$iso)]]\n"
    "  set a [list $s1 $s2 $s3]\n"
    "    \n"
    "  if [expr [lindex $v 0] > [lindex $v 1] ] {\n"
    "    set v [list [lindex $v 1] [lindex $v 0] [lindex $v 2]]\n"
    "    set a [list [lindex $a 1] [lindex $a 0] [lindex $a 2]]\n"
    "  }\n"
    "\n"
    "  if [expr [lindex $v 1] > [lindex $v 2] ] {\n"
    "    set v [list [lindex $v 0] [lindex $v 2] [lindex $v 1]]\n"
    "    set a [list [lindex $a 0] [lindex $a 2] [lindex $a 1]]\n"
    "  }\n"
    "  if [expr [lindex $v 0] > [lindex $v 1] ] {\n"
    "    set v [list [lindex $v 1] [lindex $v 0] [lindex $v 2]]\n"
    "    set a [list [lindex $a 1] [lindex $a 0] [lindex $a 2]]\n"
    "  }\n"
    "\n"
    "  if [expr [lindex $v 1] > [lindex $v 2] ] {\n"
    "    set v [list [lindex $v 0] [lindex $v 2] [lindex $v 1]]\n"
    "    set a [list [lindex $a 0] [lindex $a 2] [lindex $a 1]]\n"
    "  }\n"
    "  set z [lindex $a 2]\n"
    "  set y [lindex $a 0]\n"
    "  set x [lindex $a 1]\n"
    "  \n"
    "  set aniso [expr $z-$iso ]\n"
    "  set eta [expr double($y-$x)/$aniso ];\n"
    "  return [list $iso $aniso $eta]\n"
    "}\n"
    "\n"
    "proc csaprinc {iso aniso eta} {\n"
    "  set zz [expr $aniso + $iso]\n"
    "  set xx [expr $iso-$aniso*(1.0+$eta)/2.0]\n"
    "  set yy [expr $xx + $eta*$aniso]\n"
    "  return [list $xx $yy $zz]\n"
    "}\n"
    "\n"
    "\n"
    "proc putmatrix {m {fm \"%9.3g\"}} {\n"
    "   foreach i $m {\n"
    "     foreach j $i {\n"
    "        if {[llength $j] == 2} {\n"
    "          puts -nonewline [format \"($fm,$fm) \" [lindex $j 0] [lindex $j 1]]\n"
    "        } else {\n"
    "          puts -nonewline [format $fm $j]\n"
    "        }\n"
    "     }\n"
    "     puts \"\"\n"
    "   }\n"
    "}\n"
    "\n"
    "\n"
    "\n"
    "proc contourplot {file xlabel ylabel} {\n"
    "  set f [open $file.gnu w]\n"
    "  puts $f \"\n"
    "  set param \n"
    "  set view 0,0,1\n"
    "  set cntrparam bspline\n"
    "  set cntrparam levels 10\n"
    "  set nosurface\n"
    "  set xlabel '$xlabel'\n"
    "  set ylabel '$ylabel'\n"
    "  set contour\n"
    "  set term post\n"
    "  set output '$file.ps'\n"
    "  splot '$file' w l\n"
    "  \"\n"
    "  close $f\n"
    "  exec gnuplot $file.gnu\n"
    "  puts \"Generated: $file.ps\"\n"
    "}\n"
    "\n"
    "proc 2dplot {file xlabel ylabel {title {}}} {\n"
    "\n"
    "  set f [open $file.gnu w]\n"
    "  puts $f \"\n"
    "  set term post\n"
    "  set param\n"
    "  set view 75,20,1\n"
    "  set contour\n"
    "  set title  '$title'\n"
    "  set xlabel '$xlabel'\n"
    "  set ylabel '$ylabel'\n"
    "  set output '${file}.ps'\n"
    "  plot '$file' w l\n"
    "  \"\n"
    "  close $f\n"
    "  exec gnuplot $file.gnu\n"
    "  puts \"Generated: ${file}.ps\"\n"
    "}\n"
    "\n"
    "\n"
    "proc 3dplot {file xlabel ylabel {zrange {}}} {\n"
    "\n"
    "  set zrng {}\n"
    "  if {[llength $zrange] == 2} {\n"
    "     set zrng \"set zrange \\[[join $zrange :]\\]\"\n"
    "  }\n"
    "  set f [open $file.gnu w]\n"
    "  puts $f \"\n"
    "  set term post\n"
    "  set param\n"
    "  $zrng\n"
    "  set view 75,20,1\n"
    "  set contour\n"
    "  set xlabel '$xlabel'\n"
    "  set ylabel '$ylabel'\n"
    "  set output '${file}-3d.ps'\n"
    "  splot '$file' w l\n"
    "  \"\n"
    "  close $f\n"
    "  exec gnuplot $file.gnu\n"
    "  puts \"Generated: ${file}-3d.ps\"\n"
    "}\n"
    ;
/* EOF */
/* This file has been generated from the following Tcl source file(s):
 *    rfshapes.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char rfshapes_tcl[] = 
    "\n"
    "proc shape2fid {{shape -1} {arg1 \"\"} {arg2 \"\"}} {\n"
    "  if {$shape == -1} {\n"
    "    puts stderr \"Usage: <desc> shape2fid <shape> <name> ?<duration>? ?(-ap|-apw|-xy)?\\n\"\n"
    "    puts stderr \"where <desc>     is the SIMPSON file descriptor returned by the command.\"\n"
    "    puts stderr \"      <shape>    is the shape descriptor.\"\n"
    "    puts stderr \"      <duration> is the duration of the shape (optional). If not specified\"\n"
    "    puts stderr \"                 it defaults to the number of points in the shape.\"\n"
    "    puts stderr \"      -ap  creates the fid as amplitude/phase in the real/imaginary part.\"\n"
    "    puts stderr \"      -apw as above but with the phase wrapped to modulo 360 degrees.\"\n"
    "    puts stderr \"      -xy  creates the fid with x-y components (default).\"\n"
    "    exit\n"
    "  }\n"
    "  set type 1\n"
    "  if {$arg1 == \"-ap\" || $arg2 == \"-ap\"} {set type 2}\n"
    "  if {$arg1 == \"-apw\" || $arg2 == \"-apw\"} {set type 3}\n"
    "  set sw 1\n"
    "  if {$arg1 != \"\" && [string is double $arg1]} {\n"
    "    set sw [expr 1.0/$arg1]\n"
    "  } elseif {$arg2 != \"\" && [string is double $arg2]} {\n"
    "    set sw [expr 1.0/$arg2]\n"
    "  }\n"
    "  if { [catch {shape_len $shape} np] } {\n"
    "      puts stderr \"shape2fid error: shape does not exist\"\n"
    "      exit\n"
    "   }\n"
    "  set f [fcreate -type fid -sw $sw -np $np -ref 0]\n"
    "  set i 1\n"
    "  switch $type {\n"
    "    1 {\n"
    "      set d2r [expr 3.1415927/180]\n"
    "      foreach l [shape2list $shape] {\n"
    "        fsetindex $f $i [expr cos($d2r*[lindex $l 1])*[lindex $l 0]] \\\n"
    "                        [expr sin($d2r*[lindex $l 1])*[lindex $l 0]]\n"
    "        incr i\n"
    "      }\n"
    "    }\n"
    "    2 {\n"
    "      foreach l [shape2list $shape] {\n"
    "        fsetindex $f $i [lindex $l 0] [lindex $l 1]\n"
    "        incr i\n"
    "      }\n"
    "    }\n"
    "    3 {\n"
    "      foreach l [shape2list $shape] {\n"
    "        set ph [lindex $l 1]\n"
    "        set ph [expr fmod($ph, 360)]\n"
    "        if {$ph < 0} {set ph [expr $ph + 360]}\n"
    "        fsetindex $f $i [lindex $l 0] $ph\n"
    "        incr i\n"
    "      }\n"
    "    }\n"
    "  }\n"
    "  return $f\n"
    "}\n"
    "\n"
    "proc shape2varian {{shape -1} {name -1} {maxrf -1}} {\n"
    "  global par\n"
    "  if {$name == -1} {\n"
    "    puts stderr \"Usage: shape2varian <shape> <name> ?<maxrf>?\"\n"
    "    exit\n"
    "  }\n"
    "  if { [catch {shape_len $shape} np] } {\n"
    "      puts stderr \"shape2varian error: shape does not exist\"\n"
    "      exit 1\n"
    "   }\n"
    "  if {$maxrf == -1} {\n"
    "    set maxrf [shape_ampl $shape -max]\n"
    "  }\n"
    "  set fp [open $name w]\n"
    "  puts $fp \"# Shape generated by SIMPSON $::simpson_version\"\n"
    "  puts $fp \"# From file: $par(name)\"\n"
    "  puts $fp \"# Date: [clock format [clock seconds]]\"\n"
    "  puts $fp \"# Points in shape: $np\"\n"
    "  foreach l [shape2list $shape] {\n"
    "    set ph [lindex $l 1]\n"
    "    set ph [expr fmod($ph, 360)]\n"
    "    if {$ph < 0} {set ph [expr $ph + 360]}\n"
    "    puts $fp [format \"%10.1f %10.1f\\t 1.0\" $ph \\\n"
    "          [expr round([lindex $l 0]/double($maxrf)*1023.0)]]\n"
    "  }\n"
    "  \n"
    "}\n"
    "\n"
    "proc shape2bruker {{shape -1} {name -1} {maxrf -1}} {\n"
    "  global par\n"
    "  if {$name == -1} {\n"
    "    puts stderr \"Usage: shape2bruker <shape> <name> ?<maxrf>?\"\n"
    "    exit\n"
    "  }\n"
    "  if { [catch {shape_len $shape} np] } {\n"
    "      puts stderr \"shape2bruker error: shape does not exist\"\n"
    "      exit 1\n"
    "  }\n"
    "  if {$maxrf == -1} {\n"
    "    set maxrf [shape_ampl $shape -max]\n"
    "  }\n"
    "  \n"
    "  set fp [open $name w]\n"
    "  puts $fp \"##TITLE= Shape genereted from file: $par(name)\"\n"
    "  puts $fp \"##USAGE= Use for OC-optimized transfer/building block\"\n"
    "  puts $fp \"##JCAMP-DX= 5.00 \\$\\$ Bruker JCAMP library\"\n"
    "  puts $fp \"##DATA TYPE= Shape Data\"\n"
    "  puts $fp \"##ORIGIN= SIMPSON $::simpson_version\"\n"
    "  puts $fp \"##DATE= [clock format [clock seconds] -f \"%D\"]\"\n"
    "  puts $fp \"##TIME= [clock format [clock seconds] -f \"%T\"]\"\n"
    "  puts $fp \"##MINX= 0.000000e+00\"\n"
    "  puts $fp \"##MAXX= 1.000000e+02\"\n"
    "  puts $fp \"##MINY= 0.000000e+00\"\n"
    "  puts $fp \"##MAXY= 3.600000e+02\"\n"
    "  puts $fp \"##\\$SHAPE_EXMODE= None\"\n"
    "  puts $fp \"##\\$SHAPE_TOTROT= 0.000000e+00\"\n"
    "  puts $fp \"##\\$SHAPE_BWFAC= 0.000000e+00\"\n"
    "  puts $fp \"##\\$SHAPE_INTEGFAC= 0.000000e+00\"\n"
    "  puts $fp \"##\\$SHAPE_MODE= 1\"\n"
    "  puts $fp \"##NPOINTS= $np\"\n"
    "  puts $fp \"##XYPOINTS= (XY..XY)\"\n"
    "  foreach l [shape2list $shape] {\n"
    "    set ph [lindex $l 1]\n"
    "    set ph [expr fmod($ph, 360)]\n"
    "    if {$ph < 0} {set ph [expr $ph + 360]}\n"
    "    puts $fp [format \"%.6e, %.6e\" \\\n"
    "          [expr [lindex $l 0]/double($maxrf)*100.0] $ph]\n"
    "  }\n"
    "  puts $fp \"##END= \"\n"
    "  close $fp\n"
    "}\n"
    ;
/* EOF */
/* This file has been generated from the following Tcl source file(s):
 *    relax.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char relax_tcl[] = 
    "\n"
    "proc relax { data } {\n"
    "   global relax par\n"
    "   \n"
    "   set gm_ok { none spherical_top symmetric_top asymmetric_top }\n"
    "   set lm_ok { rigid model_free model_free_ext diffusion_on_a_cone diffusion_in_a_cone 3_sites_jump }\n"
    "   set cc_ok { shift dipole quadrupole }\n"
    "   set auto 0\n"
    "   set cross 0\n"
    "   \n"
    "   set data [split $data \"\\n\"]\n"
    "   foreach lst $data {\n"
    "     if [string match #* [string trimleft $lst]] continue\n"
    "     set nam [lindex $lst 0]\n"
    "     if ![string length $nam] continue\n"
    "\n"
    "     switch -exact $nam {\n"
    "         global_motion {\n"
    "	                set dum [lindex $lst 1]\n"
    "                        if {[lsearch -exact $gm_ok $dum]==-1} {\n"
    "			   puts stderr \"error: unknown global motion '$dum' in relax section, must be one\" \n"
    "			   set f [join $gm_ok {, }]\n"
    "			   puts stderr \"       of $f\"\n"
    "			   exit\n"
    "			}\n"
    "			if [info exists relax($nam)] {\n"
    "                          puts stderr \"error: '$nam' already exists in array relax\"\n"
    "                          exit\n"
    "                        }\n"
    "	                set relax($nam) [lrange $lst 1 9999]\n"
    "		       }\n"
    "	 shift {\n"
    "	        set dum \"$nam\\_[lindex $lst 1]\"\n"
    "		set dum2 [lindex $lst 2]\n"
    "		if {[lsearch -exact $lm_ok $dum2]==-1} {\n"
    "		   set dum1 [regsub -all {_} $dum { }]\n"
    "	           puts stderr \"error: unknown local motion '$dum2' for '$dum1' in relax section, \" \n"
    "	           puts stderr \"       must be one of [join [lrange $lm_ok 0 3] {, }],\"\n"
    "		   puts stderr \"       [join [lrange $lm_ok 4 10] {, }]\"\n"
    "		   exit\n"
    "		}\n"
    "		if [info exists relax($dum)] {\n"
    "                  puts stderr \"error: '$dum' already exists in array relax\"\n"
    "                  exit\n"
    "                }\n"
    "		set relax($dum) [lrange $lst 2 9999]\n"
    "		incr auto\n"
    "	       }\n"
    "	 dipole {\n"
    "	         set dum \"$nam\\_[lindex $lst 1]_[lindex $lst 2]\"\n"
    "		 set dum2 [lindex $lst 3]\n"
    "		 if {[lsearch -exact $lm_ok $dum2]==-1} {\n"
    "		    set dum1 [regsub -all {_} $dum { }]\n"
    "	            puts stderr \"error: unknown local motion '$dum2' for '$dum1' in relax section, \" \n"
    "		    puts stderr \"       must be one of [join [lrange $lm_ok 0 3] {, }],\"\n"
    "		    puts stderr \"       [join [lrange $lm_ok 4 10] {, }]\"\n"
    "		    exit\n"
    "		 }\n"
    "		 if [info exists relax($dum)] {\n"
    "                    puts stderr \"error: '$dum' already exists in array relax\"\n"
    "                    exit\n"
    "                 }\n"
    "		 set relax($dum) [lrange $lst 3 9999]\n"
    "		 incr auto\n"
    "	        }\n"
    "         quadrupole {\n"
    "	             set dum \"$nam\\_[lindex $lst 1]\"\n"
    "		     set dum2 [lindex $lst 2]\n"
    "		     if {[lsearch -exact $lm_ok $dum2]==-1} {\n"
    "		        set dum1 [regsub -all {_} $dum { }]\n"
    "	                puts stderr \"error: unknown local motion '$dum2' for '$dum1' in relax section, \" \n"
    "	                puts stderr \"       must be one of [join [lrange $lm_ok 0 3] {, }],\"\n"
    "		        puts stderr \"       [join [lrange $lm_ok 4 10] {, }]\"\n"
    "		        exit\n"
    "		     }\n"
    "		     if [info exists relax($dum)] {\n"
    "                       puts stderr \"error: '$dum' already exists in array relax\"\n"
    "                       exit\n"
    "                     }\n"
    "		     set relax($dum) [lrange $lst 2 9999]\n"
    "		     incr auto\n"
    "	            }\n"
    "         random_field {\n"
    "	               set dum \"$nam\\_[lindex $lst 1]\"\n"
    "		       if [info exists relax($dum)] {\n"
    "                          puts stderr \"error: '$dum' already exists in array relax\"\n"
    "                          exit\n"
    "                       }\n"
    "		       set relax($dum) [lrange $lst 2 9999]\n"
    "		       incr auto\n"
    "		      }\n"
    "         cross_correlation {\n"
    "                            puts stderr \"error: cross correlation not implemented yet!!!\"\n"
    "			    exit\n"
    "	                    if [info exists relax($nam)] {\n"
    "                               puts stderr \"error: '$nam' already exists in array relax\"\n"
    "                               exit\n"
    "                            }\n"
    "	                    set relax($nam) [lrange $lst 1 9999]\n"
    "			    incr cross\n"
    "			   }\n"
    "         default {\n"
    "	          puts stderr \"error: unknown name '$nam' in relax section, must be one\"\n"
    "		  puts \"       of global_motion, shift, dipole, quadrupole, random_field,\"\n"
    "		  puts \"       cross_correlation\"\n"
    "		  exit\n"
    "		 }\n"
    "     \n"
    "     }\n"
    "   }\n"
    "   \n"
    "   if {[lsearch -exact [array names relax] global_motion] == -1} {\n"
    "     set relax(global_motion) none\n"
    "   }\n"
    "   set relax(auto)  $auto\n"
    "   set relax(cross) $cross\n"
    "\n"
    "}\n"
    ;
/* EOF */
/* This file has been generated from the following Tcl source file(s):
 *    slave.tcl
 * on Thu Sep 23 13:41:47 CEST 2010 by zdenek
 */

char slave_tcl[] = 
	"\n"
	"proc putmatrix {m {fm \"%9.3g\"}} {\n"
	"   foreach i $m {\n"
	"     foreach j $i {\n"
	"        if {[llength $j] == 2} {\n"
	"          puts -nonewline [format \"($fm,$fm) \" [lindex $j 0] [lindex $j 1]]\n"
	"        } else {\n"
	"          puts -nonewline [format $fm $j]\n"
	"        }\n"
	"     }\n"
	"     puts \"\"\n"
	"   }\n"
	"}\n"
	"\n"
	"\n"
	"proc par {data} {  \n"
	"\n"
	"}\n"
	"\n"
	"proc spinsys {data} {  \n"
	"\n"
	"}\n"
	;
/* EOF */
TCLCODE tclcode_pointers[] = {

{"simpson", simpson_tcl},
{"ftools", ftools_tcl},
{"misc", misc_tcl},
{"rfshapes", rfshapes_tcl},
{"relax", relax_tcl},
{"slave", slave_tcl},
{"", (char*)0}
};
