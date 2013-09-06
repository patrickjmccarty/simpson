# External (user) commands and parameters
# opt::xx: External value of parameter xx
# opt::newpar: define a new parameter
# opt::fix: fix parameter = omit from optimization
# opt::release: release parameter = include in optimization
# opt::minimize: run simplex minimization
# opt::scan: scan parameter
# opt::function: specify function
# opt::setpar: set parameter to a specific value

# Internal commands and parameters
# opt::_internal::par(xx): Internal value of parameter xx
# opt::_internal::limits(xx): Bool value specifying if limits are active
# opt::_internal::limit_min(xx): Minimum allowed value of parameter xx
# opt::_internal::limit_max(xx): Maximum allowed value of parameter xx
# opt::_internal::step(xx): Initial step size for parameter xx
# opt::_internal::allparameters: List of all parameters
# opt::_internal::parameters: List of active parameters (not fixed)
# opt::_internal::function: function to minimize
# opt::_internal::wrapper: Wrapper function for the "minimize" command
# opt::_internal::int2ext: Convert from internal to external parameter
# opt::_internal::ext2int: Convert from external to internal parameter

package require optimization
package require math::linearalgebra

package provide opt 1.0

namespace eval ::opt {}
namespace eval ::opt::_internal {}
proc vecadd {a b} {
	set n [llength $a]
	for {set i 0} {$i < $n} {incr i} {
		lappend res [expr [lindex $a $i]+[lindex $b $i]]
	}
	return $res
}

proc vecsub {a b} {
	set n [llength $a]
	for {set i 0} {$i < $n} {incr i} {
		lappend res [expr [lindex $a $i]-[lindex $b $i]]
	}
	return $res
}

proc vecscale {a b} {
	if {[llength $a] == 1} {
		set n [llength $b]
		for {set i 0} {$i < $n} {incr i} {
			lappend res [expr $a*[lindex $b $i]]
		}
	} else {
		set n [llength $a]
		for {set i 0} {$i < $n} {incr i} {
			lappend res [expr $b*[lindex $a $i]]
		}
	}
	return $res
}


proc ::opt::minimize {{tol 1.0e-3}} {
	set n [llength $::opt::_internal::parameters]
	for {set j -1} {$j < $n} {incr j} {
		if [info exists row] {unset row}
		for {set i 0} {$i < $n} {incr i} {
			set nam [lindex $::opt::_internal::parameters $i]
			::opt::_internal::ext2int $nam
			set p $::opt::_internal::par($nam)
			if {$i == $j} {
				set step $::opt::_internal::step($nam)
				set p [expr $p + $step]
				if {$::opt::_internal::limits($nam)} {
					if {$p > $::opt::_internal::limit_max($nam)} {
						set p [expr $p - 2*$step]
					}
				}
			}
			lappend row $p
		}
		lappend simplex $row
	}
	
	set optimizer [optimization -simplex $simplex -function ::opt::_internal::wrapper]
	$optimizer configure -tol $tol
	set result [$optimizer start]
	if {[llength $result] == 2} {
		for {set i 0} {$i < $n} {incr i} {
			set nam [lindex $::opt::_internal::parameters $i]
			set ::opt::_internal::par($nam) [lindex $result 0 $i]
			::opt::_internal::int2ext $nam
		}
		return [lindex $result 1]
	}
	return nan
}

proc ::opt::confidence {parname} {
	set i [lsearch $::opt::_internal::parameters $parname]
	if {$i < 0} {
		puts stderr "Cannot calculate confidence interval for parameter $parname"
		return
	}
	foreach p $::opt::_internal::parameters {
		set store($p) [set ::opt::$p]
	}
	
	set origparameters $::opt::_internal::parameters
	set ::opt::_internal::parameters [lreplace $origparameters $i $i]

	set n 11
	if {$::opt::_internal::limits($parname)} {
		set delta [expr ($::opt::_internal::limit_max($parname) - $::opt::_internal::limit_min($parname))/double($n-1)]
		set p0 $::opt::_internal::limit_min($parname)
	} else {
		set delta $::opt::_internal::step($parname)
		set p0 [expr [set ::opt::$parname] - ($n-1)/2.0*$delta]
	}
	
	set sx4  0.0
	set sx3  0.0
	set sx2  0.0
	set sx   0.0
	set sx2y 0.0
	set sxy  0.0
	set sy   0.0
	if [info exists ::opt::_internal::c_data] {unset ::opt::_internal::c_data}
	for {set i 0} {$i < $n} {incr i} {
		set p [expr $p0 + $i*$delta]
		set ::opt::$parname $p
		set y [::opt::minimize]
		set sx   [expr $sx + $p]
		set sx2  [expr $sx2 + $p*$p]
		set sx3  [expr $sx3 + $p*$p*$p]
		set sx4  [expr $sx4 + $p*$p*$p*$p]
		set sy   [expr $sy + $y]
		set sxy  [expr $sxy + $p*$y]
		set sx2y [expr $sx2y + $p*$p*$y]
		lappend ::opt::_internal::c_data [list $p $y]
	}
	set ::opt::_internal::parameters $origparameters
	
	set mat [::math::linearalgebra::mkMatrix 3 3 0]
	set vec [::math::linearalgebra::mkVector 3 0]
	::math::linearalgebra::setelem mat 0 0 $sx4
	::math::linearalgebra::setelem mat 0 1 $sx3
	::math::linearalgebra::setelem mat 0 2 $sx2
	::math::linearalgebra::setelem mat 1 0 $sx3
	::math::linearalgebra::setelem mat 1 1 $sx2
	::math::linearalgebra::setelem mat 1 2 $sx
	::math::linearalgebra::setelem mat 2 0 $sx2
	::math::linearalgebra::setelem mat 2 1 $sx
	::math::linearalgebra::setelem mat 2 2 $n
	
	::math::linearalgebra::setelem vec 0 $sx2y
	::math::linearalgebra::setelem vec 1 $sxy
	::math::linearalgebra::setelem vec 2 $sy
	
	set res [::math::linearalgebra::solveGauss $mat $vec]

	foreach p $::opt::_internal::parameters {
		set ::opt::$p $store($p)
	}
	
	set a [lindex $res 0]
	if {$a > 0} {
		return [expr 2.0/sqrt($a)]
	}
	return 0
}

proc ::opt::c_data {} {
	if [info exists ::opt::_internal::c_data] {return $::opt::_internal::c_data}
}

proc ::opt::c_calculate {data} {
	set n [llength $data]
	set sx4  0.0
	set sx3  0.0
	set sx2  0.0
	set sx   0.0
	set sx2y 0.0
	set sxy  0.0
	set sy   0.0
	for {set i 0} {$i < $n} {incr i} {
		set p [lindex $data $i 0]
		set y [lindex $data $i 1]
		set sx   [expr $sx + $p]
		set sx2  [expr $sx2 + $p*$p]
		set sx3  [expr $sx3 + $p*$p*$p]
		set sx4  [expr $sx4 + $p*$p*$p*$p]
		set sy   [expr $sy + $y]
		set sxy  [expr $sxy + $p*$y]
		set sx2y [expr $sx2y + $p*$p*$y]
	}

	set mat [::math::linearalgebra::mkMatrix 3 3 0]
	set vec [::math::linearalgebra::mkVector 3 0]
	::math::linearalgebra::setelem mat 0 0 $sx4
	::math::linearalgebra::setelem mat 0 1 $sx3
	::math::linearalgebra::setelem mat 0 2 $sx2
	::math::linearalgebra::setelem mat 1 0 $sx3
	::math::linearalgebra::setelem mat 1 1 $sx2
	::math::linearalgebra::setelem mat 1 2 $sx
	::math::linearalgebra::setelem mat 2 0 $sx2
	::math::linearalgebra::setelem mat 2 1 $sx
	::math::linearalgebra::setelem mat 2 2 $n

	::math::linearalgebra::setelem vec 0 $sx2y
	::math::linearalgebra::setelem vec 1 $sxy
	::math::linearalgebra::setelem vec 2 $sy

	set res [::math::linearalgebra::solveGauss $mat $vec]

	set a [lindex $res 0]
	if {$a > 0} {
		set a [expr 2.0/sqrt($a)]
	}
	return [list $a $res]
}

proc ::opt::fix {parname} {
	set i [lsearch $::opt::_internal::parameters $parname]
	if {$i >= 0} {
		set ::opt::_internal::parameters [lreplace $::opt::_internal::parameters $i $i]
	}
}

proc ::opt::release {parname} {
	if {[lsearch $::opt::_internal::allparameters $parname] >= 0} {
		if {[lsearch $::opt::_internal::parameters $parname] < 0} {
			lappend ::opt::_internal::parameters $parname
		}
	}
}

proc ::opt::newpar {parname value step {min 1e9} {max -1e9}} {
	set ::opt::$parname $value
	set ::opt::_internal::step($parname) $step

	if {$min < $max} {
		set ::opt::_internal::limits($parname) 1
		set ::opt::_internal::limit_min($parname) $min
		set ::opt::_internal::limit_max($parname) $max
	} else {
		set ::opt::_internal::limits($parname) 0
	}

	lappend ::opt::_internal::allparameters $parname
	lappend ::opt::_internal::parameters $parname
}

proc ::opt::setpar {parname value} {
	if {[lsearch ::opt::internal::allparameters $parname] >= 0} {
		set opt::$parname $value
	}
}

proc ::opt::scan {parname} {
	set n 21
	if {$::opt::_internal::limits($parname)} {
		set delta [expr ($::opt::_internal::limit_max($parname) - $::opt::_internal::limit_min($parname))/double($n-1)]
		set p0 $::opt::_internal::limit_min($parname)
	} else {
		set delta $::opt::_internal::step($parname)
		set p0 [expr [set ::opt::$parname] - ($n-1)/2.0*$delta]
	}
	for {set i 0} {$i < $n} {incr i} {
		set p [expr $p0 + $i*$delta]
		set ::opt::$parname $p
		set f [$::opt::_internal::function]
		if {$i == 0} {
			set fmin $f
			set pmin $p
		} else {
			if {$f < $fmin} {
				set fmin $f
				set pmin $p
			}
		}
	}
	set ::opt::$parname $pmin
	return $fmin
}

proc ::opt::function {funcname} {
	set ::opt::_internal::function $funcname
}

proc ::opt::_internal::ext2int {name} {
	if {[set ::opt::_internal::limits($name)]} {
		set pext [set ::opt::$name]
		set a $::opt::_internal::limit_min($name)
		set b $::opt::_internal::limit_max($name)
		set ::opt::_internal::par($name) [expr asin(2.0*($pext-$a)/double($b-$a) - 1)]
	} else {
		set ::opt::_internal::par($name) [set ::opt::$name]
	}
}

proc ::opt::_internal::int2ext {name} {
	if {[set ::opt::_internal::limits($name)]} {
		set pint $::opt::_internal::par($name)
		set a $::opt::_internal::limit_min($name)
		set b $::opt::_internal::limit_max($name)
		set ::opt::$name [expr $a + ($b-$a)/2.0*(sin($pint)+1)]
	} else {
		set ::opt::$name $::opt::_internal::par($name)
	}
}

proc ::opt::_internal::wrapper {vertex} {
	set n [llength $::opt::_internal::parameters]
	for {set i 0} {$i < $n} {incr i} {
		set nam [lindex $::opt::_internal::parameters $i]
		set ::opt::_internal::par($nam) [lindex $vertex $i]
		::opt::_internal::int2ext $nam
	}
	return [$::opt::_internal::function]
}

