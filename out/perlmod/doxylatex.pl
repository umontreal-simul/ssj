use DoxyStructure;
use DoxyDocs;

sub latex_quote($) {
	my $text = $_[0];
	$text =~ s/\\/\\textbackslash /g;
	$text =~ s/\|/\\textbar /g;
	$text =~ s/</\\textless /g;
	$text =~ s/>/\\textgreater /g;
	$text =~ s/~/\\textasciitilde /g;
	$text =~ s/\^/\\textasciicircum /g;
	$text =~ s/[\$&%#_{}]/\\$&/g;
	print $text;
}

sub generate_doc($) {
	my $doc = $_[0];
	for my $item (@$doc) {
		my $type = $$item{type};
		if ($type eq "text") {
			latex_quote($$item{content});
		} elsif ($type eq "parbreak") {
			print "\n\n";
		} elsif ($type eq "style") {
			my $style = $$item{style};
			if ($$item{enable} eq "yes") {
				if ($style eq "bold") { print '\bfseries'; }
				if ($style eq "italic") { print '\itshape'; }
				if ($style eq "code") { print '\ttfamily'; }
			} else {
				if ($style eq "bold") { print '\mdseries'; }
				if ($style eq "italic") { print '\upshape'; }
				if ($style eq "code") { print '\rmfamily'; }
			}
			print '{}';
		} elsif ($type eq "symbol") {
			my $symbol = $$item{symbol};
			if ($symbol eq "copyright") { print '\copyright'; }
			elsif ($symbol eq "szlig") { print '\ss'; }
			print '{}';
		} elsif ($type eq "accent") {
			my ($accent) = $$item{accent};
			if ($accent eq "umlaut") { print '\"'; }
			elsif ($accent eq "acute") { print '\\\''; }
			elsif ($accent eq "grave") { print '\`'; }
			elsif ($accent eq "circ") { print '\^'; }
			elsif ($accent eq "tilde") { print '\~'; }
			elsif ($accent eq "cedilla") { print '\c'; }
			elsif ($accent eq "ring") { print '\r'; }
			print "{" . $$item{letter} . "}"; 
		} elsif ($type eq "list") {
			my $env = ($$item{style} eq "ordered") ? "enumerate" : "itemize";
			print "\n\\begin{" . $env ."}";
		  	for my $subitem (@{$$item{content}}) {
				print "\n\\item ";
				generate_doc($subitem);
		  	}
			print "\n\\end{" . $env ."}";
		} elsif ($type eq "url") {
			latex_quote($$item{content});
		}
	}
}

sub generate($$) {
	my ($item, $node) = @_;
	my ($type, $name) = @$node[0, 1];
	if ($type eq "string") {
		print "\\" . $name . "{";
		latex_quote($item);
		print "}";
	} elsif ($type eq "doc") {
		if (@$item) {
			print "\\" . $name . "{";
			generate_doc($item);
			print "}%\n";
		} else {
#			print "\\" . $name . "Empty%\n";
		}
	} elsif ($type eq "hash") {
		my ($key, $value);
		while (($key, $subnode) = each %{$$node[2]}) {
			my $subname = $$subnode[1];
			print "\\Defcs{field" . $subname . "}{";
			if ($$item{$key}) {
				generate($$item{$key}, $subnode);
			} else {
#					print "\\" . $subname . "Empty%\n";
			}
			print "}%\n";
		}
		print "\\" . $name . "%\n";
	} elsif ($type eq "list") {
		my $index = 0;
		if (@$item) {
			print "\\" . $name . "{%\n";
			for my $subitem (@$item) {
				if ($index) {
					print "\\" . $name . "Sep%\n";
				}
				generate($subitem, $$node[2]);
				$index++;
			}
			print "}%\n";
		} else {
#			print "\\" . $name . "Empty%\n";
		}
	}
}

generate($doxydocs, $doxystructure);
