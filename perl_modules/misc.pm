
package My::misc;
use strict;
use warnings;
 
use Exporter qw(import);
our @EXPORT_OK = qw(newHist addToHist igvSession);



# my $IntersObj=My::InterS->new;
sub newHist{

	my ($self,$intList)=@_;
	my $histObj={};
	my $cCat="";
	
	$$histObj{intList}=$intList;

	$cCat="<".$$intList[0];
	push @{$$histObj{catOrder}},$cCat;
	$$histObj{$cCat}=0;
	my $i=0;
	for ($i=1; $i<=$#$intList;$i++){ 
		$cCat=">=".$$intList[$i-1].", <".$$intList[$i];
		push @{$$histObj{catOrder}},$cCat;
		$$histObj{$cCat}=0;
	} 
	$cCat=">=".$$intList[$#$intList];
	push @{$$histObj{catOrder}},$cCat;
	$$histObj{$cCat}=0;

	return $histObj;
}


sub addToHist{
	my ($self,$histObj,$cVal)=@_;
	my $intList=$$histObj{intList};
	my $cCat="";
# 	print "$cVal $$intList[$#$intList]\n";
	if($cVal<$$intList[0]){$cCat="<".$$intList[0]; $$histObj{$cCat}++;}
	elsif($cVal>=$$intList[$#$intList]){$cCat=">=".$$intList[$#$intList]; $$histObj{$cCat}++;}
	else{
		my $i=1; 
		for ($i=1; $i<=$#$intList;$i++){  
			if ($$intList[$i-1] <= $cVal and $cVal < $$intList[$i]){
				$cCat=">=".$$intList[$i-1].", <".$$intList[$i];
				$$histObj{$cCat}++;
# 				print "$cCat\n";
				last;
			} 
		}
	}
	return \$cCat;
}

sub igvSession{ 
	# %IGVS igv session hash $IGVS{panelNumber}{TrackNumber}{option}=""; 
	# minimum e.g. $IGVS{1}{1}{path}="/g/data2/gd7/resources/UCSC/GRCh37/GRCh37GenomicSuperDup.bed.gz";

	use File::Basename;
	my($x,$cIGVS,$cIGVOpt,$sessionPath)=@_;
	my %IGVS=%{$cIGVS};
	my %IGVOpt=%{$cIGVOpt};
	
	my %defaultTrack;
	
	
	my $ENSEMBLEgenes=($sessionPath =~ /^http/)? "http://130.56.244.184/external/andmin/controlTracks/Homo_sapiens.GRCh37.75.gff.gz":"/g/data2/gd7/research/andmin/resources/genes/Homo_sapiens.GRCh37.75.gff.gz";
	
	$defaultTrack{"org.broad.igv.track.DataSourceTrack"}={"altColor"=>"255,0,0","autoScale"=>"true", "clazz"=>"org.broad.igv.track.DataSourceTrack", "color"=>"102,102,102", 
	"displayMode"=>"COLLAPSED", "featureVisibilityWindow"=>"-1", "fontSize"=>"10", "normalize"=>"false", "renderer"=>"BAR_CHART", 
	"sortable"=>"true", "visible"=>"true", "windowFunction"=>"mean"};

	$defaultTrack{"org.broad.igv.track.FeatureTrack"}={"altColor"=>"0,0,178","autoScale"=>"false", "clazz"=>"org.broad.igv.track.FeatureTrack", "color"=>"0,0,178", 
	"displayMode"=>"SQUISHED", "featureVisibilityWindow"=>"-1", "fontSize"=>"10", "renderer"=>"BASIC_FEATURE", 
	"sortable"=>"false", "visible"=>"true", "windowFunction"=>"count"};
	
	$defaultTrack{"DataRange"}={"baseline"=>"0.0", "drawBaseline"=>"true", "flipAxis"=>"false", "maximum"=>"66.8", "minimum"=>"0.0", "type"=>"LINEAR"};
	
	$IGVOpt{dividerFractions}="0.33,0.52,0.72,0.89" if !exists($IGVOpt{dividerFractions});
	
	# fill track options in with default values
	foreach my $cPanel (sort {$a <=> $b} keys %IGVS){
		foreach my $cTrack (sort {$a <=> $b} keys %{$IGVS{$cPanel}}){
# 			print STDERR $IGVS{$cPanel}{$cTrack}{"path"}." $cPanel}{$cTrack\n";
			my $cBaseName=basename($IGVS{$cPanel}{$cTrack}{"path"});
			my $cClazz=($cBaseName =~ /[.]bw$/)? "org.broad.igv.track.DataSourceTrack":"org.broad.igv.track.FeatureTrack";
			$IGVS{$cPanel}{$cTrack}{"clazz"}=$cClazz;
			$IGVS{$cPanel}{$cTrack}{"id"}=$IGVS{$cPanel}{$cTrack}{"path"};
			$IGVS{$cPanel}{$cTrack}{"name"}=$cBaseName if !exists($IGVS{$cPanel}{$cTrack}{"name"});
			
			foreach my $cOption (keys %{$defaultTrack{$cClazz}}){
				if(!exists($IGVS{$cPanel}{$cTrack}{$cOption})){
					$IGVS{$cPanel}{$cTrack}{$cOption}=$defaultTrack{$cClazz}{$cOption};
				}
			}
			if($cClazz eq "org.broad.igv.track.DataSourceTrack"){
				foreach my $cOption (keys %{$defaultTrack{"DataRange"}}){
					if(!exists($IGVS{$cPanel}{$cTrack}{"DataRange"}{$cOption})){
						$IGVS{$cPanel}{$cTrack}{"DataRange"}{$cOption}=$defaultTrack{"DataRange"}{$cOption};
					}
				}
			}
			
			
		}
	}
	

#### output the xml
my $xmlOUT="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
<Session genome=\"b37\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"1:197493130-197514036\" version=\"8\">
\t<Resources>\n";
	foreach my $cPanel (sort {$a <=> $b} keys %IGVS){
		foreach my $cTrack (sort {$a <=> $b} keys %{$IGVS{$cPanel}}){
			$xmlOUT.="\t\t<Resource path=\"".$IGVS{$cPanel}{$cTrack}{"path"}."\"/>\n";
		}
	}

	$xmlOUT.="\t\t<Resource path=\"$ENSEMBLEgenes\"/>\n";
    $xmlOUT.="\t</Resources>\n";
    
foreach my $cPanel (sort {$a <=> $b} keys %IGVS){

	my $cPanelName=($cPanel eq 1)?  "DataPanel":"Panel$cPanel";
	$xmlOUT.="\t<Panel height=\"390\" name=\"$cPanelName\" width=\"1133\">\n";

	foreach my $cTrack (sort {$a <=> $b} keys %{$IGVS{$cPanel}}){
		
		my $cTrackOption=""; foreach my $cOption (sort keys %{$IGVS{$cPanel}{$cTrack}}){ next if $cOption eq "DataRange"; $cTrackOption.="$cOption=\"".$IGVS{$cPanel}{$cTrack}{$cOption}."\" ";	}
		
		if(exists($IGVS{$cPanel}{$cTrack}{"DataRange"})){
		
			my $cDataROption=""; foreach my $cOption (sort keys %{$IGVS{$cPanel}{$cTrack}{"DataRange"}}){ $cDataROption.="$cOption=\"".$IGVS{$cPanel}{$cTrack}{"DataRange"}{$cOption}."\" ";	}
			$xmlOUT.="\t\t<Track $cTrackOption >\n\t\t\t<DataRange $cDataROption/>\n\t\t</Track>\n";
		}else{ $xmlOUT.="\t\t<Track $cTrackOption />\n"; }
		
	}
	$xmlOUT.="\t</Panel>\n";

}

$xmlOUT.="\t<Panel height=\"126\" name=\"FeaturePanel\" width=\"1133\">
\t\t<Track altColor=\"0,0,178\" autoScale=\"false\" color=\"0,0,178\" displayMode=\"EXPANDED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sortable=\"false\" visible=\"true\"/>
\t\t<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" colorScale=\"ContinuousColorScale;0.0;251.0;255,255,255;0,0,178\" displayMode=\"EXPANDED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" height=\"35\" id=\"b37_genes\" name=\"Gene\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">
\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"251.0\" minimum=\"0.0\" type=\"LINEAR\"/>
\t\t</Track>
\t<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" displayMode=\"SQUISHED\" featureVisibilityWindow=\"1000000\" fontSize=\"10\" id=\"$ENSEMBLEgenes\" name=\"ENSEMBL Genes\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\"/>
\t</Panel>
\t<PanelLayout dividerFractions=\"".$IGVOpt{dividerFractions}."\"/>
\t\t<HiddenAttributes>
\t\t\t<Attribute name=\"DATA FILE\"/>
\t\t\t<Attribute name=\"DATA TYPE\"/>
\t\t\t<Attribute name=\"NAME\"/>
\t\t</HiddenAttributes>
</Session>";

return \$xmlOUT;

}







