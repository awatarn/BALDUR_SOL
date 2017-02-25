#!/usr/bin/perl
#
# V1.1 POST Fortran Pre-Processor
#
# Parse fortran file and 
#
# 1) substitute ${pre}."QUOTE".${post} by '
# 2) substitute ${pre}."DOUBLEQUOTE".${post} by "
# 3) substitute ${pre}."SLASH".${post} by /
# 4) recover fortran comments 
# 5) map double spaces to single
# 6) warn about lines exceeding 72 characters
#
# A. Pletzer April 14 1999
#
# V1.0 accepts -free option: APLET August 26 1999
# V1.1 turn back ^ in first 6 fields into spaces. 
#
# USAGE:
# 
# cat file | postfpp [-132] [-free] > newfile
#
# options:
#
# -132 Generate warning if line width exceeds 132 columns (72 is default)
# -free Source is f90 free-form. This overrides the -132 option.
#
####################################################################
$MAX_COLUMNS = 72; # max number of columns

foreach $opt (@ARGV){
    if( $opt =~ /\-132/ || $opt =~ /\-free/){ $MAX_COLUMNS = 132; }
}

$FREE_FLAG = 0;
foreach $opt (@ARGV){
        if( $opt =~ /-free/ ) {$FREE_FLAG = 1;}
}


#FILE=STDIN;
#if($#ARGV > 0){
#    $file = $ARGV[$#ARGV];
#    open(FILE, $file) ||  die "Can't open $file\n\n";
#}

#
# substitution grammar
#
$pre = "\@";
$post = "\@";

$quote = "\'";
$doublequote = "\"";
$slash = "\/";
$exclam = "!";
$leftparen = "\(";
$rightparen = "\)";
$rightbracket = "\}";

@meta = ('\(', '\)', '\$', '\<', '\>', 
	 '\{', '\}', '\^', '\[', '\]', 
	 '\+', '\.', '\*', '\-', '\@', '\|', '\?');

@special = ($quote, $doublequote, $slash, $rightbracket);

$substitute{$special[0]} = ${pre}."QUOTE".${post};
$substitute{$special[1]} = ${pre}."DOUBLEQUOTE".${post};
$substitute{$special[2]} = ${pre}."SLASH".${post};
$substitute{$special[3]} = ${pre}."RIGHTBRACKET".${post};

#
#
while(<STDIN>){
#
# don't touch preprocessor command lines (if any left)
#
  if(!/^ *$/){                      # CAL: rm lines with blanks
                                    #      Perl 4
  if(!/^\s*#.*/){
#
     s/^\s*\n//g;                   # rm empty lines
     s/^\s*${pre}_0${post}\n//g;    # rm empty lines
#
# remove extra spaces -- only in non-verbatim fields
#
     ($beyond6 = $_) =~ s/^(.{6})//;
     $first6 = $1;
     @fields = split( /${pre}verbatim\{[^\}]*\}/, $beyond6 );
     foreach $f (@fields){
	 ($newf = $f) =~ s/[ ]+/ /g;
	 $newf =~ s/(\w)\s+\(/$1(/g; # remove gcc's space between sub/fct names and (
	#
        # protect meta-characters
	#
	foreach $m (@meta){
	    $f =~ s/($m)/\\$1/g;
	}
				 #print "\n\n   f=$f\n";
				 #print "newf=$newf\n\n";
	 #s/${first6}${f}/${first6}${newf}/;
	 #$_ = $first6 . $beyond6;
     }
#
# remove verbatim wrapper
#
    s/${pre}verbatim\{([^\}]*)\}/$1/g;
#
# map ${pre}EXCLAM${post} back to !
#
    s/${pre}EXCLAM${post}/${exclam}/g;
#
# special characters -- which may interfere with preprocessor --
# are mapped back
#
    foreach $c (@special){
	s/$substitute{$c}/$c/g;
    }
#
# warn about line overflow
#
  s/${pre}_(\d+)${post}//;
  $ncols = $1;
  @chars = split(//);
  if($#chars > $ncols && $#chars > $MAX_COLUMNS ){
      print STDERR  "*WARNING* ncols=$ncols MAX_COLUMNS=$MAX_COLUMNS \n";
      print STDERR  "*WARNING* inflated line->$_"; 
      print STDERR  "*WARNING* The script will now attempt to squeeze out blanks. ".
                    "This can potentially yield wrong code!\n";
      $oldwidth = $#chars;
      if(/\"/){
	  print STDERR "*ERROR* Double quote(s) indicating start/end of comment\n";
	  print STDERR "*ERROR* will exit with error code 1";
	  exit(1);
      } elsif (/\'/){
	  print STDERR "*ERROR* Single quote(s) indicating start/end of comment\n";
	  print STDERR "*ERROR* will exit with error code 2";
	  exit(2);
      } else {
	  #
          # compression is safe?
	  #
	  s/(.{6})([^ ]*)[ ]+/${1}${2} /g;
	  s/(\w)[ ]+\(/${1}(/g;
	  @chars = split(//);
	  if( $#chars >= $oldwidth || $#chars > $MAX_COLUMNS ){
	      print STDERR "*ERROR* Compression failed with error code 3\n";
	      exit(3);
	  }
      }
	  
	  
	  
  }
  
  }
    # turn back  ^ in first 6 columns into blanks

    if( $FREE_FLAG == 0 ) {
       if( $_ !~ /^[Cc\*!]/ ){
     # 
     # translate spaces located in first 6 columns into ^
     # (some gcc precompilers remove them)
     #
     /^(.{6})(.*)$/;
     $first_part = $1;
     $second_part= $2;
     $first_part =~ s/\^/ /g;
     $_ = $first_part . $second_part . "\n";
    }
   }
    print;
}

}










