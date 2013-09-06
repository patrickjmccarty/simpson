lappend ::auto_path ../ 
package require opt

proc y {} {
	set y [expr pow($opt::a,2) + pow($opt::b-1,2)+1]
	puts -nonewline [format "%10.5f %10.5f %10.5f\015" $opt::a $opt::b $y]
	flush stdout
	after 100
	return $y
}

opt::newpar a 5 0.1 -1 10
opt::newpar b -5.0 1

opt::function y

puts "Scan a"
set min [opt::scan a]
puts -nonewline "\nBest values: "
y
puts {}

puts "Scan b"
set min [opt::scan b]
puts -nonewline "\nBest values: "
y
puts {}

puts "Minimize while fixing a"
opt::fix a
set min [opt::minimize]
puts {}

puts "Minimizing"
opt::release a
set min [opt::minimize 1e-5]

puts "\nFinal parameters"

puts "a: $opt::a"
puts "b: $opt::b"
puts "Y: $min"

