#!/user/bin/perl -w

use Bio::SeqIO;

$input_pattern = "TCCAAATGAAGTCATTATCAAA";

#sampled:
#$input_pattern = "ATGGCTGCCGGCGGCCGAGGTC";

$fasta_file = "/Users/kdimatteo/bin/ben/gene/data/CCDS_nucleotide.current.fna";
#$fasta_file = "/Users/kdimatteo/bin/ben/gene/data/ringer.fasta";
$start = (times)[0];


#$FILE = "IAV_p1_fixed.txt";
#data/CCDS_nucleotide.current.fna

@p = split(/ */, $input_pattern);
$seqio_obj = Bio::SeqIO->new(-file => $fasta_file, -format=>"fasta");
  
while ($seq_obj = $seqio_obj->next_seq){   
    
	#pattern for verification at end of the run 
	@L = (); 
	
	#string from the sequence for verification
	@Q = (); 
	
	#mismatches, up to 3
	$digit_errors = 0;
	
	#spaced between errors, up to 6
	$space_errors = 0;
	
	$subs = 0;
	$seqIndex = 1;
	$i = 0;
	$j = 0;
	
    # print the sequence   
	#print ">>>>" , $seq_obj->seq,"\n";
	
	@s = split(//, $seq_obj->seq);
	
	#print scalar(@Q);
	#print scalar(@p);
	
	while ($i<scalar(@s)){
    	
    	
    	
    	if (scalar(@Q) == scalar(@p)){
            print "==============================================\n";
        	last;
        }
            
       	#print "s:", $i, ":", $s[$i], "?=", "p:", $j, ":", $p[$j], "(e:", $digit_errors, ")", join("", @L), "(", scalar(@L), ")\n";
		
		if ($s[$i] eq $p[$j]){
        	push(@L, $p[$j]);
            push(@Q, $s[$i]);
            $i++;
            $j++;
                
         } elsif (($s[$i] eq 'T') && ($p[$j] eq 'C')){
         	push(@L, $p[$j]);
            push(@Q, $s[$i]);
            $i++;
            $j++;
            $subs++;
         } else {
         	# no match
            if (length(@L) > 0){
            	if ($digit_errors <= 3){
                	$digit_errors++;
                    $space_errors++;
                    push(@L, $p[$j]);
                    push(@Q, $s[$i]);
                    $i++;
                    $j++;
                 }else{
                 	$i = $seqIndex;
                    $seqIndex++;
                    $j = 0;
                    $digit_errors = 0;
                    @L = ();
                    @Q = ();
                }
            }else{
            	$i++;
                $j = 0;
                # nothin in L anway
            }
         }
    }
    
    #print "checking", $seq_obj['description'], "\n";
    
    if (join("", @L) eq $input_pattern){
    	print join("", @L) , " === ", $input_pattern, "\n";
    	print "[@][SUCCESS]: ", $seq_obj->id, " : ", $seq_obj->description, "\n";
        print "Found this in the sequence:", join("", @Q), "\n\n";
    } else {
    	#print "() fail",  $seq_obj->description, "\n\n";
    }            
}
        

$end = (times)[0];
printf "completed in %.2f seconds", $end - $start;