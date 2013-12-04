#!/usr/bin/env perl

# Packages
use LWP::Simple;
use URI::URL;

# Defaults
$TRUE=1;
$FALSE=0;
$CLOBBER=$TRUE;
$DTOR=3.14159265359/180;
%MAGKEYS=("SDSS","UGRIZ");

## Not doing anything with DEBUG right now
#$DEBUG=$TRUE;
#$DEBUG=$FALSE;

$addnums=$FALSE;
$addmags=$FALSE;
$radius=12;
$catalog="SDSS";
$color="green";
$regfile="";
$fullcat="";
$p60=0;
$keck=0;
@now=gmtime(time());

# Following formula is correct only to ~1day
$epoch=(1900+$now[5]+$now[7]/365.0);

($exec)=($0=~/\/?([^\/]+)$/);

################
# Command Line #
################

&usage if (@ARGV<3);

while(@ARGV) {
    $_=shift(@ARGV);
    if (/^[-\+]\D/) {
	# Options
	if (/^[-\+]c/i) {
	    # Color for regions
	    $color=shift(@ARGV);
	    print STDERR "# Using color $color for regions\n";
	}elsif(/^[-\+]r/i) {
	    # Search radius
	    $radius=shift(@ARGV);
	    print STDERR "# Setting search radius to $radius arcmin\n";
	}elsif (/^[-\+]m/i) {
	    # Output magnitudes
	    $addmags=$TRUE;
	}elsif (/^[-\+]f/i) {
	    # Write a region file
	    $regfile=shift(@ARGV);
	    print STDERR "# Writing region file $regfile\n";
	}elsif (/^[-\+]s/i) {
	    # Save full catalog
	    $fullcat=shift(@ARGV);
	    print STDERR "# Saving full catalog to $fullcat\n";
        }elsif (/^[-\+]p/i) {
            # P60 Query
            $p60=1;
            print STDERR "# P60 Query: Rejecting Galaxies\n";
        }elsif (/^[-\+]k/i) {
            # Keck Query
            $keck=1;
            print STDERR "# Keck Query: Rejecting Galaxies\n";
        }elsif (/^[-\+]u/i) {
            # Subtraction Query w/o HPM Stars
            $subtnopm=1;
            print STDERR "# Subtraction Query: Rejecting Galaxies and HPM Stars\n";
	}elsif (/^[-\+]b/i) {
            # Subtraction Query w/o HPM Stars
            $subtall=1;
            print STDERR "# Subtraction Query: Rejecting Galaxies\n";
	}else{
	    # Unrecognized options
	    print STDERR "# Unrecognized option $_... ignoring.\n";
	}
    }else{
	# Arguments
	push(@args,$_);
    }
}

&usage if (@args<3);
($ra,$dec,$out)=@args;

######################
# Parse Input Coords #
######################

# Parse RA
if ($ra=~/[:\s]/) {
    # Sexagesimal RA
    ($rah,$ram,$ras)=($ra=~/(\d+)[\s:]+(\d+)[\s:]+([\d\.]+)/);
    $rahd=ten($rah,$ram,$ras);
    $rad=$rahd*15;
}else{
    $rad=$ra;
    $rahd=$ra/15;
    ($rah,$ram,$ras)=&sixty($rahd);
}

# Parse Dec
if ($dec=~/[:\s]/) {
    # Sexagesimal Dec
    ($dcd,$dcm,$dcs)=($dec=~/([\+\-]?\d+)[\s:]+(\d+)[\s:]+([\d\.]+)/);
    $decd=ten($dcd,$dcm,$dcs);
}else{
    $decd=$dec;
    ($dcd,$dcm,$dcs)=&sixty($decd);
}


#####################
# Execute Retrieval #    
#####################

$astrom=&getsdss;

exit($astrom);

######################################################################

sub getsdss {
    # gets astrometry catalog from SDSS
    # puts data into file $out, which contains:
    #     index ra DEC distance [mag magu]
    # where distance = output distance in arcsec
    #            mag = catalog magnitude if requested
    #           magu = catalog magnitude uncertainty

    # Size in kb of HTML doc's returned for failed requests
    my $failsize=120;
    my $failsafe;

    my $cmd, $url;
    my $RA,$DEC,$dra,$ddec,@text,$i,$r,$content;

    # Calculate min/max boundaries for Sloan
    if ($catalog eq "SDSS") {
	($ramin,$dcmin)=&cooshift($rad,$decd,-$radius,-$radius);
	($ramax,$dcmax)=&cooshift($rad,$decd,+$radius,+$radius);
    }

    # Sloan Catalog Retrieval
    $query="select p.objID,p.ra,p.dec,p.psfMag_u,p.psfMagErr_u," .
	          "p.psfMag_g,p.psfMagErr_g,p.psfMag_r,p.psfMagErr_r," . 
	          "p.psfMag_i,p.psfMagErr_i,p.psfMag_z,p.psfMagErr_z" .
		  " FROM fGetObjFromRect($ramin,$ramax,$dcmin,$dcmax) n,";

    if ($p60 eq 1) {
       $query=$query . "Star p WHERE n.objID=p.objID AND p.psfMag_r < 22.0";
    }elsif ($keck eq 1) {
       $query=$query . "Star p WHERE n.objID=p.objID AND p.psfMag_r < 24.5";
    }elsif ($subtnopm eq 1) {
       $query=$query . "Star p, ProperMotions a WHERE n.objID=p.objID AND " .
              "a.objID=n.objID AND (p.psfMag_r BETWEEN 12.0 AND 18.0) AND " .
              "SQRT(POWER(a.pmRa,2) + POWER(a.pmDec,2)) < 30 AND " .
	      "(p.flags & 0x10) = 0 AND (p.flags & 0x40000) = 0 AND " .
              "(p.flags & 0x40) = 0 AND (p.flags & 0x1000) " .
              " = 0 AND (p.nChild = 0 OR (p.flags & 0x08) = 0)";
    }elsif ($subtall eq 1) {
	$query=$query . "Star p WHERE n.objID=p.objID AND (p.psfMag_r " .
               "BETWEEN 12.0 AND 18.0) AND (p.nChild = 0 OR " .
               "(p.flags & 0x08) = 0)";
              #"(p.flags & 0x10) = 0 AND (p.flags & 0x08) = 0 AND (p.flags & " .
              #"0x40) = 0 AND p.nChild = 0"; 
    }else {
       $query=$query . "PhotoObj p WHERE n.objID=p.objID";
    }

    $cmd="sdss_sql -q \"$query\"";

    # Report on progress
    if ($catalog eq "SDSS") {
	print "Retrieving SDSS catalog with command \'$cmd\'\n";

	unless (open(CMD,"$cmd |")) {

	    killoff("Failed to spawn SDSS catalog retrieval command \'$cmd\'\n");
	}
	@catlines=<CMD>;
	close(CMD);
	$content=join("\n",@catlines);
    }else{
	# Report on progress
	unless ($catalog eq "SDSS") {
	    printout("\# Fetching $url\n");
	    unless (defined ($content=get($url))) {
		killoff("Could not get $url $!\n");
	    }
	}
    }

    # Save all catalog data to file if requested
    if ($fullcat) {
	open(FCAT,"> $fullcat");
	print FCAT $content;
	close(FCAT);
    }

    # Write a region file if requested
    if ($regfile) {
	if(!open(REG,"> $regfile")){
	    print STDERR "# Failed to open region file $regfile for writing\n";
	    $regfile="";
	}else{
	    print REG "\# Region file format: DS9 version 3.0\n";
	    print REG "global color=green font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n";
	}
    }

    # Reformat output
    @text=split(/\n/,$content);
    
    if (-e $out && !($CLOBBER)) {
	killoff("Clobber set to no and file $out exists\n");
    }

    if (!open(OUT,"> $out")) {
	killoff("Couldn't open $out for writing: $!",__LINE__);
    }

    my $i=0;
    foreach $_ (@text) {

	if ($catalog eq "SDSS") {
	    if (/^\d{12}/) {
		# Data line
		s/^\s+//;
		@data=split(/,/,$_);

		($RA,$DEC,$UMAG,$UMAGU,$GMAG,$GMAGU,$RMAG,$RMAGU,
		 $IMAG,$IMAGU,$ZMAG,$ZMAGU) = @data[1..12];

		++$i;
		$rdist{$i}=&separation($RA,$DEC,$rad,$decd)*3600;
		$mags{$i}=join(",",@data[3..12]);
		$dataline{$i}="$RA $DEC 0.2 0.2";
	    }
	}
    }

    my @lines=(1..$i);
    @lines = sort { $rdist{$a} <=> $rdist{$b} } @lines;
    $ndig=1+int(log(0.1+@lines)/log(10));
    my $n=0;
    @magkeys=split(//,$MAGKEYS{$catalog});
    foreach $i (@lines) {
	++$n;
	$lab=sprintf("%s-%0${ndig}d",$catalog,$n);
	@mags=split(/,/,$mags{$i});
	# Catalog line
	$catline="";
	$catline=sprintf("%${ndig}d ",$n) if $addnums;
	$catline.=sprintf("%s-%0${ndig}d %s %.2f",
			  $catalog,$n,$dataline{$i},$rdist{$i});
	$catline.=" @mags" if $addmags;
	print OUT "$catline\n";
	# Region file line
	if ($regfile) {
	    @d=split(/\s+/,$dataline{$i});
	    $regline="fk5;ellipse($d[0],$d[1],$d[2]\",$d[3]\",0) \# text=\{$lab\} color=$color";
	    for ($ii=0;$ii<@magkeys;++$ii) {
		$jj=2*$ii;
		$regline.=" ${magkeys[$ii]}MAG=\{$mags[$jj]\}";
		++$jj;
		$regline.=" ${magkeys[$ii]}MAGU=\{$mags[$jj]\}";
	    }
	    print REG "$regline\n";
	}
    }
    close(OUT);
    close(REG) if ($regfile);

    return (-s $out > 0);
    }

######################################################################

sub usage {
    my $diestr=(@_ ? $_[0] : "");
    $diestr .= "Usage:  $exec [-d <catalog>] [-r <radius/arcmin>] [-e <epoch>] [-s <savefile>] [-f <regfile>] [-m(ags)] [-n(umber)] [-c <color>] <ra> <dec> <outfile>\n";
    die $diestr;
}

sub sixty {
    my $dcml=$_[0];
    my $sign = ($dcml<0 ? -1 : 1);
    $dcml*=$sign;
    my $hexi=int($dcml);
    my $hexmd=60*($dcml-$hexi);
    my $hexm=int($hexmd);
    my $hexs=sprintf("%6.3f",60*($hexmd-$hexm));
    return ($sign*$hexi,$hexm,$hexs);
}

sub ten {
    my ($din,$min,$sec)=@_[0..2];
    my $sign = ($din=~/^-/ ? -1 : 1);
    $din *= $sign;
    my $deg=$sign*($din+($min+$sec/60)/60);
    return $deg;
}

sub max {
    my @a=@_;
    my $max=shift(@a);
    while (@a) {
	$el=shift(@a);
	$max=$el if ($el>$max);
    }
    return $max;
}

sub acos2 {
    # Range:  [0, pi]
    my ( $inp ) = @_;
    while( $inp > 1.0 ) { $inp -= 1.0; }
    while( $inp < -1.0 ) { $inp += 1.0; }
    my $outp = atan2( sqrt( 1.0 - ( $inp * $inp ) ), $inp );
}

sub cooshift {

    # Shift coordinates (RA, Dec) by (dRA, dDec) arcmin
    my ($raz, $dcz, $dra, $ddc) = @_[0..3];
    my ($newra, $newdc) = (0,0);
    my $DEGREES = 57.2957795130823; # Degrees in a radian.
    my $ARCMIN=60.0; # arcmin per degree

    # RA shift
    if (abs($dcz)<90) {
	$newra=$raz + ($dra/($ARCMIN*cos($dcz/$DEGREES)));
    }else{
	# Degenerate case
	$newra=$raz;
    }

    # Check RA limits
    while ($newra<0) {
	$newra += 360.0;
    }
    while ($newra>=360.0) {
	$newra -= 360.0;
    }

    # Dec shift
    $newdc=$dcz + $ddc/$ARCMIN;

    # Check Dec limits
    if (abs($newdc)>90) {
	# Oops, over the pole
	$flipdc=180.0 - abs($newdc);
	$newdc = ($newdc>0 ? $flipdc : -$flipdc);
    }

    return ($newra,$newdc);
}

sub separation { 
    # Input RA and DEC in decimal degrees.  Return angular separation, degrees.
    my ( $ra1, $dec1, $ra2, $dec2 ) = @_;
    my $degrees = 57.2957795130823; # Degrees in a radian.
    my $x1 = cos( $ra1 / $degrees ) * cos( $dec1 / $degrees );
    my $y1 = sin( $ra1 / $degrees ) * cos( $dec1 / $degrees );
    my $z1 = sin( $dec1 / $degrees );
    my $x2 = cos( $ra2 / $degrees ) * cos( $dec2 / $degrees );
    my $y2 = sin( $ra2 / $degrees ) * cos( $dec2 / $degrees );
    my $z2 = sin( $dec2 / $degrees );
    my $theta = acos2( $x1*$x2 + $y1*$y2 + $z1*$z2 ) * $degrees;
    return ($theta);
}   

sub printout {
    my $inp=shift;
    print STDERR $inp;
}

sub killoff {
    my $message=shift;
    my $line=0;
    $line=shift if ($#_ >= 0);
    my $s;

    $s=$message;
    $s.= " at " . __FILE__ . " line " . $line . ".\n" if ($message !~ /\n$/ && $line > 0);    
    die $s;
}

######################################################################
######################################################################

__END__

