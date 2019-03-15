use strict;
use warnings;
use Data::Dumper;
use IO::File;

 

my $docompocount=1;
my %h;
my $inpf=$ARGV[0];
my $cfilepath=$ARGV[1];
my $extid=9249;
my $doextract=0;
my $maxnei=150;

my @remspe=(
 
);

sub keepId($){
 my($idi)=@_;
 my $fnd=1;
 for(my $i=0;$i<=$#remspe;$i++){
	my $speid=$remspe[$i];
        if($idi=~/^$speid/){
		$fnd=0;
		 last;
        }
 }
return $fnd;
 
}
 


 
my $alllinks = {};
open FILE, "<", "$inpf" or die $!;
while (<FILE>){
	 chomp;
	 my $l=$_;
	 my @r=split/\t/,$l;  
	 my $node=$r[0];
	 my $node2=$r[1]; 

	 if(!defined $node || length($node)<3 ){
		        die($node ."too short");
	 }
	 if(!defined $node2 ||length($node2)<3 ){
		        die($node2 ."too short");
	 }
	  if(keepId($node)==1 && keepId($node2)==1){
	   $alllinks->{$node}->{$node2} = 1;
	   $alllinks->{$node2}->{$node} = 1;
	 }
	 my @ar1 = keys %{ $alllinks->{$node} };
	 my @ar2 = keys %{ $alllinks->{$node2} };
	 if($#ar1>$maxnei){ die ("too mush element ar1 ");}
	 if($#ar2>$maxnei){ my $st=$node2."=".Dumper \$alllinks->{$node2}; die ("too mush element ar2 ".$st);}
}
#print Dumper \$alllinks;
#exit();
my $df=IO::File->new(">subset.txt");
#print "----------------------\n";
my $currentccid = 0;
my %visited = ();
 
foreach my $node (keys %$alllinks) {
    	 #print "***********KEY: $node\n";
    if (defined  $visited{$node}  ){}
    else{	
            #print "\tnot visited : $node\n";
	    # start the next segment
	    
	    $currentccid++;

	    #my @to_visit;
	    my @tvt;
	    push @tvt,$node;

	    my @nta;	
	    my %hta;
	    push @nta,$node;
	    $hta{$node}=1;
	    my %vst;
	    while ($#nta>=0) {
	      my $nodeta = shift @nta;
	      if(! exists $vst{$nodeta} ){
		      $vst{$nodeta}=1;
		      my @ngb = keys %{ $alllinks->{$nodeta} };
		      for(my $i=0;$i<=$#ngb;$i++){
			 my $n=$ngb[$i];
		         if(! exists $vst{$n} ){
				if(! exists $hta{$n}){
					push @nta,$n;
					$hta{$n}=1;
				}
			 }
	      	     }
		}
	    }	
             foreach my $g(keys %vst){
			$visited{$g}=$currentccid;
	   } 
 
   }
 
  #if($ccid==$extid){
	#last;
  #}
}

$df->close();
my %cc;
my %hh;
foreach my $k(keys %visited){
	my $v=$visited{$k};
	
       if($doextract==1){
	 if($v==$extid){
                $hh{$k}= 1;

	 }
	}
	if($docompocount==1){
		my $ct=0;
		if(exists $cc{$v}){
			$ct=$cc{$v};
		}
		$ct++;
		$cc{$v}= $ct;
	}
     	print $k."\t".$v."\n";
};
if($docompocount==1){
	my $cf=IO::File->new(">".$cfilepath);
	my $maxmbc=0;
	my $largestcc=0;
	print  $cf "----cc members count-----\n";
	print  $cf "cc_id\tmember_count\n";
	foreach my $cn(sort keys %cc){
		my $countmb=$cc{$cn};
		print  $cf $cn."\t".$countmb."\n";
		if($maxmbc<$countmb){$maxmbc=$countmb;$largestcc=$cn};
	}
	print  $cf "----max members count:-----\n";
	print  $cf "----member_count:".$maxmbc."(cc_id ".$largestcc.")";
	$cf->close();	


};

if($doextract==1){
 my $cf2=IO::File->new(">".$cfilepath.".ext");
 my $dum=Dumper \%hh;
 print  $cf2 $dum."\n------------------\n";
 foreach my $el(sort keys %hh){
	 my @ar = keys %{ $alllinks->{$el} };
	 for(my $j=0;$j<=$#ar;$j++){
		my $el2=$ar[$j];
		 print  $cf2 $el."\t".$el2."\n";
	 }
 }
 $cf2->close();
}

my $eFile=IO::File->new(">End.txt");
$eFile->close();

