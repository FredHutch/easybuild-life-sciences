use strict;

my $readinput = shift;
my $fastqoutput = shift;

open(READS, $readinput);
open(FASTQ1, ">$fastqoutput"."_1.fq");
open(FASTQ2, ">$fastqoutput"."_2.fq");

my $readid = '';
my $readseq = '';
my $maxseq = 0;
my $maxid = '';
while(my $line = <READS>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		if($readseq){
			printFQ($readid, $readseq);
			if(length($readseq) > $maxseq){$maxseq = length($readseq);$maxid = $readid;}
		}
		$readid = $1;
		$readseq = '';
	}else{
		$readseq .= $line;
	}
}
printFQ($readid, $readseq);
print $maxseq."\n".$maxid."\n";

sub printFQ
{
	my $readid = shift;
	my $readseq = shift;
	
	my $nbfrag = 2;
	my $maxlen = 400;
	my $minlen = 35;
	my $nbPairs = 1;
	
	my $lenseq = length($readseq);
	if($lenseq < $minlen){return;}
	while($lenseq > $maxlen){
		if($readid eq 'F0DRRT203F7J2D')
		{
			print $lenseq."\n";
		}
		$nbPairs++;
		$nbfrag = $nbfrag * 2;
		$lenseq = int($lenseq/2);
	}
	my $lenFrag = int(length($readseq)/$nbfrag);
	my @readseqs;
	
	my $curseq;
	for (my $i = 0; $i < $nbfrag; $i++)
	{
		if($i < ($nbfrag-1))
		{
			$curseq = substr($readseq, $i*($lenFrag),$lenFrag);
		}else{
			$curseq = substr($readseq, $i*($lenFrag));
		}
		if($readid eq 'F0DRRT203F7J2D')
		{
			print "OK $i $lenFrag\n$curseq\n";
		}
		printFQout($readid,$i,$curseq);
	}

}

sub printFQout
{
	my $readid = shift;
	my $count = shift;
	my $readseq = shift;

	if($count % 2 == 0)
	{
		if($readid eq 'F0DRRT203F7J2D'){print "FQ1\n";}
		print FASTQ1 '@'.$readid."_".$count."\n";
		print FASTQ1 $readseq."\n";
		print FASTQ1 "+\n";
		print FASTQ1 ('I' x length($readseq))."\n";
	}else{
		if($readid eq 'F0DRRT203F7J2D'){print "FQ2\n";}
		print FASTQ2 '@'.$readid."_".$count."\n";
		print FASTQ2 $readseq."\n";
		print FASTQ2 "+\n";
		print FASTQ2 ('I' x length($readseq))."\n";
	}
}