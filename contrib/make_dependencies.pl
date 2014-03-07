# Copyright W. Bangerth, University of Heidelberg, 1998, 1999, 2000, 2001, 2002


# Make a dependency file tree
# usage: make_dep -Iinc_path1 -Iinc_path2 ... -Sobject_suffix

# This program makes for each of the given files a makefile dependency
# list, also considering nested includes. It only considers included
# files which are located in the given include pathes (you can give any
# number of pathes). The output looks like this:




# list of include pathes (from command line)
@include_path = ();
# list of files to be checked (from command line)
@input_files  = ();
# associate list of files together with the files they include
%include_files = ();


# fill include paths
while ($ARGV[0] =~ /^-I/) {
    $_ = shift;
    /^-I(.*)/;
    $_ = $1;
    if (m![^/]$!) {
	@include_path = (@include_path, $_ . "/");
    } else {
	@include_path = (@include_path, $_);
    }
}

# get object suffix
$_ = shift;
/^-S(.*)/;
$object_suffix = $1;

#Now i'll make the C++ source extension a variable
$cplusplusext=".cpp";

#fill list of files to be processed
while ($ARGV[0]) {
    @input_files = (@input_files, shift);
}


foreach $file (@input_files) {
    make_include_tree ($file);
};


foreach $file (keys %include_files) {
    # complete list of included files by nested include files
    foreach $include (split(' ',$include_files{$file})) {
	complete ($file, $include);
    }
}

# print dependency list
foreach $file (sort @input_files) {
    $file =~ /(.*)$cplusplusext/;

    # replace the .C or .c or .cc with .$object_suffix
    $rulename = $file;
    $rulename =~ s/$cplusplusext?\s*$/.$object_suffix/g;


    
    @include_file_list = sort (split (' ', $include_files{$file}));

    # write rule for the .o file
    print "$rulename:";
    print "\\\n    $file";
    foreach $f (@include_file_list) {
	print "\\\n    $f";
    }
    print "\n";
}








# complete the list of included files by the files
# included by a file included by the original one
sub complete {
    local ($file, $include) = ($_[0], $_[1]);
    foreach $second_include (split(' ',$include_files{$include})) {
        # check whether $second_include is in the list of filenames
	# $include_file{$file}. in order to avoid that characters
	# in the filename and/or path of $second_include are
	# interpreted as special characters in the regexp matches,
	# we have to escape all such characters (well, there
	# may be more, but I hope that no one uses them in filenames).
	my $pattern = $second_include;
	$pattern =~ s/\+/\\+/g;
	$pattern =~ s/\*/\\*/g;
	$pattern =~ s/\?/\\?/g;
	$pattern =~ s/\(/\\(/g;
	$pattern =~ s/\)/\\)/g;
	$pattern =~ s/\[/\\[/g;
	$pattern =~ s/\]/\\]/g;
	$pattern =~ s/\{/\\{/g;
	$pattern =~ s/\}/\\}/g;
	if (! ($include_files{$file} =~ $pattern)) {
	    #second_include not yet in list of included files
	    $include_files{$file} =
		join(' ', $second_include, $include_files{$file});
	    complete ($file, $second_include);
	}
    }
}




# make the include tree for a file
sub make_include_tree {
    local ($filename) = $_[0];

    open (FILE, $filename);
    while (<FILE>) {
	# look out for include statements
	if (/^#\s*include\s*(["])([^"]*)["]/) {
	    local($include) = $2;
	    # included by "..." Try to find real path
	    
	    # first find out the path of the present file
	    # and then join the path to the included one
	    for $include_dir (@include_path) {
		$if = $include_dir . $include;
		if (-r $if) {
		    # file found; delete ./ at beginning
		    if ($if =~ m!^\./(.*)!) {
			$if = $1;
		    }
		    $include_files{$filename} =
			    join (' ', $if, $include_files{$filename});
		}
	    }
	}
    }

    # for each file included here: make up include tree itself
    for $if (split(/ /,$include_files{$filename})) {
	# if include file list not yet made up
	if (! defined ($include_files{$if})) {
	    make_include_tree($if);
	}
    }
}
