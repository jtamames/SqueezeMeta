# File:  xfig.pm
# Author: Mihai Pop
#
# Routines that print out xfig output.
#

package AMOS::xfig;

use strict;

BEGIN {
    use Exporter ();
    use vars     qw(@EXPORT @EXPORT_OK @ISA %EXPORT_TAGS);

    @ISA         = qw(Exporter);
    @EXPORT      = qw(&print_xfig_header 
		      &print_solid_line
		      &print_box
		      &print_diamond
		      &print_tick
		      &print_x
		      &print_text
		      &print_horiz_text
		      &print_low_text
		      &print_circle
		      $BLACK   $BLUE   $GREEN $CYAN $RED 
		      $MAGENTA $YELLOW $WHITE $FILL $NOFILL);
    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();
}

use vars @EXPORT;
use vars @EXPORT_OK;

# define some global variables
$FILL    = 20;
$NOFILL  = -1;

$BLACK   =  0;
$BLUE    =  1;
$GREEN   =  2;
$CYAN    =  3;
$RED     =  4;
$MAGENTA =  5;
$YELLOW  =  6;
$WHITE   =  7;

# print_xfig_header
#
# prints an xfig header to standard output
sub print_xfig_header {
    my($ppi) = @_;

    print "\#FIG 3.2\n";
    print "Landscape\n";
    print "Center\n";
    print "Inches\n";
    print "Letter\n";
    print "100.00\n";
    print "Multiple\n";
    print "-2\n";
    print "$ppi 2\n";
}

# print_solid_line 
#
# prints a line between the given x1, y1, x2, y2 coordinates
sub print_solid_line {
    my($x1, $y1, $x2, $y2) = @_;

# Fields are: 
#    2 - polyline
#    1 - polyline
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    0 - join style
#    0 - cap style
#   -1 - radius
#    0 - forward arrow
#    0 - backward arrow
#    2 - number of points
    print "2 1 0 1 0 0 100 0 -1 0.000 0 0 -1 0 0 2\n";
# point coords
    print "\t$x1 $y1 $x2 $y2\n";
}

# print_box
#
# prints a box centered at the given coordinates
sub print_box {
    my($x, $y, $color, $fill) = @_;

# Fields are: 
#    2 - polyline
#    2 - box
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    0 - join style
#    0 - cap style
#   -1 - radius
#    0 - forward arrow
#    0 - backward arrow
#    2 - number of points
    print "2 2 0 1 $color $color 100 0 $fill 0.000 0 0 -1 0 0 5\n";
# point coords
    printf ("\t%d %d %d %d %d %d %d %d %d %d\n", $x - 50, $y - 50,
	    $x - 50, $y + 50, $x + 50, $y + 50, $x + 50, $y - 50,
	    $x - 50, $y - 50);
}

# print_diamond 
#
# prints a diamond centered at the given coordinates
sub print_diamond {
    my($x, $y, $color, $fill) = @_;

# Fields are: 
#    2 - polyline
#    2 - box
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    0 - join style
#    0 - cap style
#   -1 - radius
#    0 - forward arrow
#    0 - backward arrow
#    2 - number of points
    print "2 2 0 1 $color $color 100 0 $fill 0.000 0 0 -1 0 0 5\n";
# point coords
    printf ("\t%d %d %d %d %d %d %d %d %d %d\n", $x - 50, $y,
	    $x, $y + 50, $x + 50, $y, $x, $y - 50,
	    $x - 50, $y);
}

# print_tick
#
# prints tick mark at coordinate and with specified color and width
sub print_tick {
    my($x, $y, $color, $size) = @_;

# Fields are: 
#    2 - polyline
#    1 - line
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    0 - join style
#    0 - cap style
#   -1 - radius
#    0 - forward arrow
#    0 - backward arrow
#    2 - number of points
    print "2 1 0 1 $color 0 100 0 -1 0.000 0 0 -1 0 0 2\n";
# point coords
    printf ("\t%d %d %d %d\n", $x, $y - $size, $x, $y + $size);
}

# print_x
# 
# prints an x mark at coordinate and with specified color
sub print_x {
    my($x, $y, $color) = @_;

# Fields are: 
#    2 - polyline
#    1 - line
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    0 - join style
#    0 - cap style
#   -1 - radius
#    0 - forward arrow
#    0 - backward arrow
#    2 - number of points
    print "2 1 0 1 $color 0 100 0 -1 0.000 0 0 -1 0 0 2\n";
# point coords
    printf ("\t%d %d %d %d\n", $x - 50, $y - 50, $x + 50, $y + 50);

# next line
    print "2 1 0 1 $color 0 100 0 -1 0.000 0 0 -1 0 0 2\n";
    printf ("\t%d %d %d %d\n", $x + 50, $y - 50, $x - 50, $y + 50);
}

# print_text
#
# prints a string vertically at the given coordinates
sub print_text {
    my($x, $y, $string, $color) = @_;

# Fields are: 
#    4 - polyline
#    0 - left justified
#    0 - color
#  100 - depth    
#    0 - pen style
#    0 - font
#   12 - font size
# 1.5708 - angle (radians)
#    4 - font flags
#    h - height
#    l - length
#    x - x
#    y - y
#    string
    print "4 0 $color 100 0 0 12 1.5708 4 135 ", 89 * length($string), 
    " ", $x + 68, " ", $y - 68, " $string\\001\n";
}

# print_text
#
# prints a string vertically at the given coordinates
sub print_low_text {
    my($x, $y, $string, $color) = @_;

# Fields are: 
#    4 - polyline
#    2 - right justified
#    0 - color
#  100 - depth    
#    0 - pen style
#    0 - font
#   12 - font size
# 1.5708 - angle (radians)
#    4 - font flags
#    h - height
#    l - length
#    x - x
#    y - y
#    string
    print "4 2 $color 100 0 0 12 1.5708 4 135 ", 89 * length($string), 
    " ", $x + 68, " ", $y + 68, " $string\\001\n";
}

# print_horiz_text
#
# prints a string vertically at the given coordinates
sub print_horiz_text {
    my($x, $y, $string, $color) = @_;

# Fields are: 
#    4 - polyline
#    1 - center justified
#    0 - color
#  100 - depth    
#    0 - pen style
#    0 - font
#   12 - font size
# 1.5708 - angle (radians)
#    4 - font flags
#    h - height
#    l - length
#    x - x
#    y - y
#    string
    print "4 1 $color 100 0 0 12 0 4 135 ", 89 * length($string), 
    " ", $x, " ", $y, " $string\\001\n";
}

# print_circle 
#
# draws a circle centered at the given coordinates in the specified color
sub print_circle {
    my($x, $y, $color, $fill) = @_;

# Fields are: 
#    1 - ellipse
#    3 - circle
#    0 - solid
#    1 - thickness
#    0 - black (line)
#    0 - black (fill)
#  100 - depth    
#    0 - pen style
#   -1 - area fill
# 0.000 - style val ??
#    1 - direction
#    0 - angle
#    cx,cy - center coords
#    rx,ry - radius
#    sx,sy - first point
#    ex,ey - last point
    print "1 3 0 1 $color $color 100 0 $fill 0.000 1 0.000 $x $y ",
    "50 50 $x $y ", $x + 50, " $y\n";
}

END {}
