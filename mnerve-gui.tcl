#!/usr/bin/env wish



#---- Auxilary procedures --------------------------------


proc framed-label-and-entry {parent name text} {
# packs a label and a text entry into a frame
    set fullname .$parent.$name
    pack [frame $fullname] -side top -fill x
    label $fullname.l -width 15 -text $text -padx 0 -anchor w
    entry $fullname.e -width 15 -relief sunken
    pack $fullname.l .$parent.$name.e -side left

}

proc make-keys {} {
# make command arguments for mrunner.lua
    set y {}
    global entry_options
    foreach opt $entry_options {
	global [lindex $opt 0]
	set optname [lindex $opt 0]
	if {[lindex $opt 3] != 0 && \
		  [set $optname] != {}} {
	    set optname [lindex $opt 0]
	    lappend y [lindex $opt 3] [set $optname]
	}
    }
    lappend y -e
    lappend y [stims->lua]
    return $y
}

# Run simulation
proc Run {} {
    global input runbut out-file modeltime output
    global command xdata compartments yvectors

    upvar layout lo

    set dictfile $lo-dict.gp

    xdata delete :
    foreach name [array names compartments -glob "*.*"] {
	global y_${name}
	y_$name delete :
    }


    set keys [make-keys]
    puts "$command mrunner.lua $keys > ${out-file}"
    
    #set input [open "| $command mrunner.lua $keys "]
    if [catch {open "| $command mrunner.lua $keys "} input] { 
	bgerror "Couldn't start lua process"
    } else {
	set output [open ${out-file} w]

	fconfigure $input -blocking 0; # needed for win32
    
	$runbut config -text "Stop" -command Stop \
	    -activebackground red \
	    -background pink
	fileevent $input readable Log

	if { [catch {open $dictfile r} dict]} {

	    after 15000 set dict fail
	
	    while 1 {
		after 20
		if [file exists $dictfile] {
		    set dict [open $dictfile]
		    break
		}
	    }

	    vwait dict 
	}
	
	if {$dict != "fail"} {
	    read-gp-dict $dict compartments
	    close $dict
	    .controls.addplot config -command \
		{add-plot-window compartments plot_windows}
	    setup-yvectors yvectors compartments
	    .controls.addplot config -state active
	}

    }
}

# log model time from the output file
proc Log {} {
    global compartments
    global input modeltime output
    global xdata yvectors
    global max-points
    if [eof $input] { Stop 
    } else {
	gets $input line
	if {[string length $line] > 0} {
	    puts $output $line
	    set modeltime [format "%4.3f" [lindex $line  0]]
	    set dlength [xdata length]
	    if {$dlength>${max-points}} {
		    xdata delete 0
	    }
	    set xdata(++end) $modeltime
	    foreach name [array names compartments -glob "*.*"] {
		# this must be slow: random acces on a list 
		# todo: think of better line parsing
		global y_${name}
		if {$dlength>${max-points}} {
		    y_${name} delete 0
		}
		set y_${name}(++end) [lindex $line [set yvectors(${name})]]
	    }
	}
    }
}

# Stop the running lua/luajit process
proc Stop {} {
    global runbut input output
    global xdata ydata
    global tcl_platform
    if [regexp -nocase {win} [set tcl_platform(platform)]] {
	    catch {exec tskill [pid $input]} 
    }
    catch {close $input}
    catch {close $output}
    

    $runbut config -text "Run" -command Run \
	-background lightgreen -activebackground green
}

# Delete stimulation row
proc remove-stim {id} {
    global stimuli stimulation_params
    destroy $stimuli($id,frame)
    foreach p $stimulation_params {
	unset stimuli($id,$p)
    }
    unset stimuli($id,frame)
}

# Add new stimulation row
proc add-stim {} {
    global stimuli stimulation_params
    
    incr stimuli(id)
    set id $stimuli(id)
    pack [set stimuli($id,frame) [frame .stimuli.body.$id]]
    
    foreach p $stimulation_params {
	pack [entry .stimuli.body.$id.$p \
		  -textvariable stimuli($id,$p) \
		  -width 5] -side left
    }
    
    pack [button .stimuli.body.$id.del -text "x"\
	      -command "remove-stim $id"] -side left
    
}

proc stims->lua {} {
    global stimuli stimulation_params
    set res "return "
    append res "{"
    for {set id 0} {$id <= $stimuli(id)} {incr id} {
	set stim [array get stimuli $id*]
	if {$stim != {} } {
	    set st {}
	    append st "{"
	    foreach p $stimulation_params {
		append st "$p=$stimuli($id,$p), "		
	    }
	    append st "}, "
	    append res $st
	}
    }
    append res "} "
}

proc read-gp-dict {fobj compartments} {
# reads gnuplot dictionary
    upvar $compartments comps
    set comp_prev ""
    set comps(names) [list]
    while {[gets $fobj line] >= 0} {
	if [regexp {([a-z]+[0-9]+)_([a-z]+[0-9]+)_(.+) = ([0-9]+)} \
		$line all tag comp var id] {
	    if {$comp != $comp_prev} {
		lappend comps(names) $comp
	    }
	    set comp_prev $comp
	    set comps($comp.$var) \
		[list $var [expr $id - 1]]
	}
    }
}


proc make-var-menu {id compartments} {
    upvar $compartments comps
    
    set win .plot-$id

    menubutton $win.menubar.pl -text {Variable} \
	-menu $win.menubar.pl.m

    pack $win.menubar.pl -side left

    set m [menu $win.menubar.pl.m]
    set i 0
    foreach key $comps(names) {
	$m add cascade -label $key -menu $m.sub$i
	set sub [menu $m.sub$i]
	incr i 
	foreach name [lsort [array names comps $key.*]] {
	    $sub add checkbutton -label [lindex [set comps($name)] 0] \
		-variable varflags($id-$name) \
		-command "add/remove-element $id $name ng"
	}
    }
}

proc make-save-graph-menu {id graph} {
    set win .plot-$id
    menubutton $win.menubar.sv -text {Save} \
	-menu $win.menubar.sv.m
    pack $win.menubar.sv -side left
    set m [menu $win.menubar.sv.m]
    $m add command -label {Eps} -command "save-eps $graph"
    $m add command -label {Bitmap} -command "save-bitmap $graph"
}

proc save-eps {graph} {
    set savename [tk_getSaveFile \
		      -defaultextension .eps \
		      -filetypes {
			  {{Encapsulated postscript} {.eps}}
			  {{All files} * } 
		      }]
    if {$savename != ""} {
	$graph postscript configure \
	    -landscape yes -center yes -maxpect yes \
	    -decorations no
	$graph postscript output $savename
    }
}

proc save-bitmap {graph} {
	$graph snap -format EMF CLIPBOARD
}

proc UnZoom {graph} {
    $graph axis configure x y -min {} -max {}
}

proc Zoom {graph x1 y1 x2 y2 } {
    if {$x1 > $x2} {
	$graph axis configure x -min $x2 -max $x1
    } elseif {$x2 > $x1} {
	$graph axis configure x -min $x1 -max $x2
    }
    if {$y1 > $y2} {
	$graph axis configure y -min $y2 -max $y1
    } elseif {$y2 > $y1} {
	$graph axis configure y -min $y1 -max $y2
    }
    
}

proc RegionStart {graph x y} {
    global x0 y0
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    $graph marker create line \
	-coords ""  \
	-name regionline \
	-dashes dash -xor yes 
    set x0 $x 
    set y0 $y
}

proc slide-start {graph x y} {
    global x0 y0
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    set x0 $x
    set y0 $y
}

proc slide-move {graph x  y} {
    global x0 y0
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    

    set xmin [$graph axis cget x -min]
    set xmax [$graph axis cget x -max]
    set ymin [$graph axis cget y -min]
    set ymax [$graph axis cget y -max]

    set xdelta [expr $x0 - $x]
    set ydelta [expr $y0 - $y]
    
    if {$xmin != ""} {
	
	set xmin [expr $xmin + $xdelta]
	set xmax [expr $xmax + $xdelta]
	set ymin [expr $ymin + $ydelta]
	set ymax [expr $ymax + $ydelta]

	$graph axis configure x -min $xmin -max $xmax
	$graph axis configure y -min $ymin -max $ymax
    }
}

proc RegionMotion {graph x y} {
    global x0 y0 
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    $graph marker configure regionline \
	-coords "$x0 $y0 $x0 $y $x $y $x $y0 $x0 $y0" \
}

proc RegionEnd {graph x y} {
    global x0 y0 
    $graph marker delete regionline
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    Zoom $graph $x0 $y0 $x $y
}



proc update-graph-xy0 {graph x y} {
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    upvar graphxy gxy
    set gxy [format "%4.3f,%4.3f" $x $y]
}

proc update-graph-xy {graph x y label} {
    set x [$graph axis invtransform x $x]
    set y [$graph axis invtransform y $y]
    set gxy [format "%4.3f,%4.3f" $x $y]
    $label configure -text $gxy
}


set colors {black blue green red magenta yellow}

proc last-dir {} {
    lindex [split [pwd] / ] end
}

proc add-plot-window {compartments plot_windows} {
    global varflags graphs_count
    upvar $plot_windows plots
    upvar $compartments comps
    blt::vector create ydata
    incr plots(id)
    set id $plots(id)
    set win [toplevel .plot-$id -width 600 -height 600]

    wm title $win "[last-dir]: plot-$id"    

    set graphs_count($id) 0

    pack [frame $win.menubar] -fill x

    make-var-menu $id comps

    set graph [blt::graph .plot-$id.g -width 400 -height 260]
    pack $graph -fill both
    
    $graph grid on

    bind $graph <ButtonPress-1> {slide-start %W %x %y}
    bind $graph <B1-Motion> {slide-move %W %x %y}

    bind $graph <ButtonPress-3> {RegionStart %W %x %y}
    bind $graph <B3-Motion> {RegionMotion %W %x %y}
    bind $graph <ButtonRelease-3> {RegionEnd %W %x %y}
    bind $graph <ButtonRelease-2> {UnZoom %W}


    make-save-graph-menu $id $graph

    bind $win <g> "$graph grid toggle"
    #bind $win <l> "$graph legend toggle"
    bind $win <a> "UnZoom $graph"

    set graphxy [label $win.menubar.graphxy \
		     -relief sunken \
		     -text "undefined"] 
    pack $graphxy -side right

    bind $graph <Motion> \
	"update-graph-xy %W %x %y $graphxy"

    return $graph

}



proc setup-yvectors {yvectors compartments} {
#ODO: colt blt vectors into an array
    upvar $yvectors yvs $compartments comps
    foreach name [array names comps -glob "*.*"] {
	set yvs($name) [lindex [set comps($name)] 1]
	global y_$name
	if {![info exists y_$name]} {
	    blt::vector create y_${name}}
    }
}

proc add/remove-element {plotid name ngraphs} {
    global xdata varflags yvectors graphs_count
    global colors

    set ng [set graphs_count($plotid)]
    
    foreach n [array names comps -glob "*.*"] {
	global y_${n}
    }
    if $varflags($plotid-$name) {
	.plot-$plotid.g element create $name -symbol "" \
	    -color [lindex $colors [expr $ng % [llength $colors]]] \
	    -xdata xdata -ydata y_$name
	incr graphs_count($plotid)
    } else {
	#incr graphs_count($plotid) -1
	.plot-$plotid.g element delete $name
    }
}

#---------------------------------------------------------

### Main program starts here ###

#package require platform
package require BLT

# Set window title

set model_suffix  [lindex $argv 0]

wm title . "MGUI-$model_suffix: [last-dir]"

# Container frame for run parameters
labelframe .options -text "Parameters"

# Container frame for run  controls
frame .controls -borderwidth 5

# Container frame for stimulation protocols
labelframe .stimuli -text "Stimulation"

# Container frame for time logging
labelframe .log -text "Log"

pack .options .stimuli .log .controls -side top -fill x

# Control buttons
button .controls.quit -text "Quit" -command exit

set plot_windows(id) -1

button .controls.addplot -text "Add plot"

global xdata
blt::vector create xdata; # this one will be shared

#if {![catch {open } dictfile]} {
#    read-gp-dict $dictfile compartments
#    .controls.addplot config -command {add-plot-window \
#					   compartments plot_windows}
#    setup-yvectors yvectors compartments##
#
#} else {
#    .controls.addplot config -state disabled
#}

    .controls.addplot config -state disabled
    
set runbut [button .controls.run -text "Run" \
		-background lightgreen -activebackground green \
		-command Run]

pack .controls.run  .controls.addplot .controls.quit  -side left

## Option entries ##

set model_suffix  [lindex $argv 0]

set entry_options [list \
		       {start-time  "Start time:" 0 -t} \
		       {stop-time  "Stop time:" 1000 -T} \
		       {stepper "Stepper func:" rkc_a -i}\
		       {init-h "Init. time step:" 1e-3 -h}\
		       {rkc-s "RKC dimension:" 200 -s} \
		       {config "Config file:" dc-conf -c} \
		       {model "Model file:" dc-mod -m} \
		       "layout {Nerve layout} $model_suffix -L"\
		       {load-file "Load state from:" "" -l}\
		       "{save-file} {Save state to:} $model_suffix-state.lua -S"\
		       {save-period "Save period:" 1e4 -w}\
		       {max-points "Max pts/graph:" 5e4 0}\
		       {out-file "Output file:" out.0 0}\
		       {command "Run with:" "lua" 0}
		      ]

# Create entries for the available options
foreach option $entry_options {
    set optname [lindex $option 0]
    set $optname [lindex $option 2]
    framed-label-and-entry options $optname [lindex $option 1]
    .options.$optname.e config -textvariable $optname
}

if [regexp -nocase {win} [set tcl_platform(platform)]] {
	set command lua
} else {set command "lua"}


## Stimulation ##

set stimulation_params {amp start stop int width}

pack [frame .stimuli.head] [frame .stimuli.body] -side top -anchor w

foreach p $stimulation_params {
    pack [label .stimuli.head.$p -text $p -width 5] -side left
}

set stimuli(id) -1

set addstim [button .stimuli.head.add -text "+" -command add-stim]
pack $addstim -side left

## A bit of logging ##
set modeltime 0

label .log.lab -text "Model time:"
label .log.time -width 15 -relief sunken -textvariable modeltime

pack .log.lab .log.time -side left


