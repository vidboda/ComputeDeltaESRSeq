package sequence;
use strict;
use REST::Client;

sub new {
	my ($class, $chr, $pos, $ref, $alt, $strand) = @_;
        my $self;
	if ($chr =~ /^([\dXY]{1,2})$/o) {$chr = "chr$1"}
	$self = {
                'chr' => $chr,
		'pos' => $pos,
		'ref' => $ref,
		'alt' => $alt,
		'strand' => $strand,
		#'start' => $start,
		#'end' => $end,
		#'seq' => $seq
		};
	bless ($self, $class);
	return $self;
}
sub setChr {my ($self, $chr) = @_;$self->{chr} = $chr}
sub getChr {my $self = shift;return $self->{chr}}
sub setPos {my ($self, $pos) = @_;$self->{pos} = $pos}
sub getPos {my $self = shift;return $self->{pos}}
sub setRef {my ($self, $ref) = @_;$self->{ref} = $ref}
sub getRef {my $self = shift;return $self->{ref}}
sub setAlt {my ($self, $alt) = @_;$self->{alt} = $alt}
sub getAlt {my $self = shift;return $self->{alt}}
sub setStrand {my ($self, $strand) = @_;$self->{strand} = $strand}
sub getStrand {my $self = shift;return $self->{strand}}
sub setStart {my ($self, $start) = @_;$self->{start} = $start}
sub getStart {my $self = shift;return $self->{start}}
sub setEnd {my ($self, $end) = @_;$self->{end} = $end}
sub getEnd {my $self = shift;return $self->{end}}
sub setwtSeq {my ($self, $seq) = @_;$self->{wtseq} = $seq}
sub getwtSeq {my $self = shift;return $self->{wtseq}}
sub setmtSeq {
	my ($self, $seq, $alt) = @_;
	substr($seq,5,1,$alt);
	$self->{mtseq} = $seq;
}
sub getmtSeq {my $self = shift;return $self->{mtseq}}
#sub set {my ($self, ) = @_;$self->{} = }
#sub get {my $self = shift;return $self->{}}

sub getSurroundings {
	my $self = shift;
	my $client = REST::Client->new();
	#define start, end
	$self->setStart($self->getPos()-5);
	$self->setEnd($self->getPos()+5);
	my $response = $client->GET("http://togows.org/api/ucsc/hg19/".$self->getChr().":".$self->getStart()."-".$self->getEnd());
	if ($client->responseCode() == 200) {
		my $data = $client->responseContent();
		if (length($data) == 11) {
			#my $seq = $client->responseContent();
			#print $client->responseContent();
			if ($self->getStrand() eq '-') {
				my $seqrev = reverse $data;
				$seqrev =~ tr/acgtACGT/tgcaTGCA/;
				$data = $seqrev;
			}
			#we must check strandness of nucleotides
			if ($self->getRef() eq substr($data,5,1)) {
				#strandness is coherent
				$self->setwtSeq($data);
				$self->setmtSeq($data, $self->getAlt());
			}
			elsif(&compl($self->getRef()) eq substr($data,5,1)) {
				#strandness is not coherent, we complement alt
				$self->setwtSeq($data);
				$self->setmtSeq($data, &compl($self->getAlt()));
			}
			else {die "ERROR the wt nucleotide provided is not coherent with the genomic position for variant ".$self->getChr()." ".$self->getPos()." ".$self->getRef()." ".$self->getAlt()." : ".$data."\n"}
		}
		else {die "\nERROR: togows.org returned a wrong numebr of nucleotides for variant ".$self->getChr()." ".$self->getPos()." ".$self->getRef()." ".$self->getAlt()." : ".$data."\n"
}
	}
	else {
		die "\nERROR: togows.org badly respond for variant ".$self->getChr()." ".$self->getPos()." ".$self->getRef()." ".$self->getAlt()." http code: ".$client->responseCode()."\n"
	}	
}

sub toPrint {
	my $self = shift;
	print $self->getChr()." ".$self->getPos()." ".$self->getRef()." ".$self->getAlt()." ".$self->getStart()." ".$self->getEnd()." ".$self->getwtSeq()." ".$self->getmtSeq()."\n";
}

sub compl {
	my $nt = shift;
	if ($nt eq 'A') {return 'T'}
	elsif ($nt eq 'T') {return 'A'}
	elsif ($nt eq 'G') {return 'C'}
	elsif ($nt eq 'C') {return 'G'}
}

sub ESR {
	my $self = shift;
	my (@wt_score, @mt_score);
	for (my $i = 0;$i < 6;$i++) {
		my $hexa_wt = substr($self->getwtSeq(),$i,6);
		#my $data = `grep -e '$hexa_wt' data/hexamers.txt`;
		my @out = split(/\t/, `grep -e '$hexa_wt' data/hexamers.txt`);
		push @wt_score, $out[1];
		my $hexa_mt = substr($self->getmtSeq(),$i,6);
		@out = split(/\t/, `grep -e '$hexa_mt' data/hexamers.txt`);
		push @mt_score, $out[1];
		#print "$hexa_wt\t$hexa_data\n";
	}
	my ($wtESR, $mtESR);
	foreach (@wt_score) {$wtESR += $_}
	foreach (@mt_score) {$mtESR += $_}
	my $ESR = $mtESR - $wtESR;
	$ESR = sprintf('%.3f', $ESR);
	$wtESR = sprintf('%.3f', $wtESR);
	$mtESR = sprintf('%.3f', $mtESR);
	return ($wtESR, $mtESR, $ESR);
}


1;