use DoxyStructure;

sub process($) {
	my $node = $_[0];
	my ($type, $name) = @$node[0, 1];
	my $command;
	if ($type eq "string") { $command = "String" }
	elsif ($type eq "doc") { $command = "Doc" }
	elsif ($type eq "hash") {
		$command = "Hash";
		for my $subnode (values %{$$node[2]}) {
			process($subnode);
		}
	}
	elsif ($type eq "list") {
		$command = "List";
		process($$node[2]);
	}
	print "\\" . $command . "Node{" . $name . "}%\n";
}

process($doxystructure);
