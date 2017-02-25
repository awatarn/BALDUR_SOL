#!/usr/bin/perl
#
# V1.1 PRE Fortran Pre-Processor
#
# Parse fortran file and 
#
# 1) substitute ' by ${pre}."QUOTE".${post}
# 2) substitute " by ${pre}."DOUBLEQUOTE".${post}
# 3) substitute / by ${pre}."SLASH".${post}
# 4) identify preserve fortran comments 
# 5) capitalize fortran command lines
#
# A. Pletzer April 13 1999
#
# V1.0 added -free option: APLET August 26 1999
# V1.1 turn spaces in first 6 fields into ^ to correct for new gcc 
#      space compression behaviour.
#
# USAGE:
# 
# prefpp [-free] file > newfile
#
# OPTIONS:
#
# -free Source file is f90 free-form.
#
####################################################################
$MAX_EXCLAM = 10; # max number of exclamation marks on a line

$FREE_FLAG = 0;
foreach $opt (@ARGV){
        if( $opt =~ /-free/ ) {$FREE_FLAG = 1;}
}
$file = $ARGV[$#ARGV];

open(FILE, $file) ||  die "Can't open $file\n\n";

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


@special = ($quote, $doublequote, $slash);

@meta = ('\(', '\)', '\$', '\<', '\>', 
	 '\{', '\}', '\^', '\[', '\]', 
	 '\+', '\.', '\*', '\-', '\@', '\|', '\?');

$substitute{$special[0]} = ${pre}."QUOTE".${post};
$substitute{$special[1]} = ${pre}."DOUBLEQUOTE".${post};
$substitute{$special[2]} = ${pre}."SLASH".${post};
#$substitute{$special[3]} = ${pre}."EXCLAM".${post};
#
# verbatim delimiting fields
#
@verbatim = (${quote}, ${doublequote});
#
# exclamation mark. treat differently when in verbatim delimiters
#
$exclamation = "(\{[^\}]*)${exclam}([^\}]*\})";
#
# comments
#
@comments = ("^(!)(.*)\$");
if( $FREE_FLAG == 0 ){
    @comments = ("^([Cc\*])(.*)\$"); # f77 comments are skipped
                                     # ! type comments are treated 
                                     # separately because they can
                                     # be anywhere.
	     #"^([^!]*!)(.*)\$");
	     #"^([^!|\'|\"]*!)(.*)\$");
}

$quote_terminated = 1;
$doublequote_terminated = 1;
while(<FILE>){

if (! /^ *$/) {        # 05/16/01 CAL: for Perl 4

if( $FREE_FLAG == 0) {
   if( $_ !~ /^[Cc\*!]/ ) {
      if( $_ !~ /^\s*#/) {
     # 
     # translate spaces located in first 6 columns into ^
     # (some gcc precompilers remove them)
     #
     s/^\t/       /o;
     /^(.{6})(.*)$/;
     $first_part = $1;
     $second_part= $2;
     $first_part =~ s/[ ]/^/g;
     $_ = $first_part . $second_part . "\n";
     }
  }
}


    
#
# skip preprocessor directives
#
  if(!/^\s*#.*/){
     #
     # store line width
     # 
     @chars = split(//);
     $ncols = $#chars;
#
# } will have a special meaning. In order not to have it entangled
# with } in strings, we substitute..
#
     s/$rightbracket/${pre}RIGHTBRACKET${post}/g;
#
# detect comments. First ! indicates a comment except in between '' or ""
#
    foreach $c (@comments){
        s/$c/${1}${pre}verbatim{$2}/;
    }
     $is_a_comment = 0;
     if( $_ =~ /^[Cc\*]/ && $FREE_FLAG == 0 ){ $is_a_comment = 1;}
     if( $_ =~ /^\s*!/ && $FREE_FLAG == 1 ) { $is_a_comment = 1;}
     if( $is_a_comment == 0) {
#
# the following are logical flags
#
    $in_comment = 0;
    $in_quote = 0;
    $in_doublequote = 0;

    @chars = split(//,$_);
    $line = "";
    $k = 0;
    foreach $c (@chars){
       if( $k == 0 && $quote_terminated == 0 ){ 
	   $in_quote = 1;
	   $c = ${pre}."verbatim{".${c};
       }
       if( $k == 0 && $doublequote_terminated == 0 ){ 
	   $in_doublequote = 1; 
	   $c = ${pre}."verbatim{".${c};
       }
       if( $c eq "'"  
	  && $in_comment + $in_quote + $in_doublequote == 0){
           #
           # a ' string starts
           #
           $c = ${pre}."verbatim{".${c};
           $in_quote = 1;
       }
       if( $c eq '"'  
	  && $in_comment + $in_quote + $in_doublequote == 0){
           #
           # a " string starts
           #
           $c = ${pre}."verbatim{".${c};
           $in_doublequote = 1;
       }
       if( $c eq "'"  && $in_quote == 1) {
           #
           # a ' string ends
           #
           $c .= "}";
           $in_quote = 0;
	   $quote_terminated = 1;
       }
       if( $c eq '"'  && $in_doublequote == 1) {
           #
           # a " string ends
           #
           $c .= "}";
           $in_doublequote = 0;
	   $doublequote_terminated = 1;
       }
          if( $c eq "!" && $in_comment + $in_doublequote + $in_quote == 0 ){
           #
           # start comment
           #
           $c = "!".${pre}."verbatim{";
           $in_comment = 1;
       }
       if( $k == $#chars  && $in_comment == 1){
	   #
	   # finish up comment
	   #
	   $c = "}\n";
       }
       if( $k == $#chars  && $in_quote == 1 ){
	   #
           # string flows over to next line
           #
           $quote_terminated = 0;
	   $c = "}\n";
       }
       if( $k == $#chars  && $in_doublequote == 1 ){
	   #
           # string flows over to next line
           #
           $doublequote_terminated = 0;
	   $c = "}\n";
       }
       $line .= $c;
       $k++;
       #print "+++$line\n"; 
    }
    $_ = $line;
     }
#
# special characters -- which interfere with cpp -- are mapped
#
    foreach $c (@special){
	s/$c/$substitute{$c}/g;
    }
#
# capitalize non-verbatim text
#
    @cap = split(/${pre}verbatim\{[^\}]*\}/);
    foreach $c (@cap){
        ($C = $c) =~ tr/a-z/A-Z/;
	#
        # protect meta-characters
	#
	foreach $m (@meta){
	    $c =~ s/($m)/\\$1/g;
	}
	#
        # capitalize
        #
        s/$c/$C/;
    }
    # 
    # add line width for info
    #
    s/^(.*)$/${1}${pre}_${ncols}${post}/;
  }
    print;
}
}
