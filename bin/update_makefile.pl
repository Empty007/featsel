my $date_string = localtime ();

print "Loading the lists of algorithms and cost functions...";
my $list_of_cost_functions_file = "ListOfCostFunctions.txt";
my $list_of_algorithms_file = "ListOfAlgorithms.txt";
my $algorithm_is_in_the_list = 0;
open (LIST, $list_of_cost_functions_file) or 
  die "Error: could not open $list_of_cost_functions_file for reading!\n";
my @list_of_function_codes;
my @list_of_function_class_names;

while (<LIST>)
{
  chomp $_;
  if ($_ =~ /^\s*(\w+)\s+(\w+)\s+/) # $1 == code, $2 == class name
  {
    push @list_of_function_codes, $1;
    push @list_of_function_class_names, $2; 
  }
  elsif ($_ =~ /^\s*$/)
  {
    # Blank line; do nothing.
  }
  else
  {
    die "Error: could not parse line of $list_of_cost_functions_file!\n";
  }
}
close (LIST);

open (LIST, $list_of_algorithms_file) or 
  die "Error: could not open $list_of_algorithms_file for reading!\n";
my @list_of_algorithm_codes;
my @list_of_algorithm_class_names;
while (<LIST>)
{
  chomp $_;
  if ($_ =~ /^\s*(\w+)\s+(\w+)\s*$/) # $1 == code, $2 == class name
  {  
    push @list_of_algorithm_codes, $1;
    push @list_of_algorithm_class_names, $2; 
  }
  elsif ($_ =~ /^\s*$/)
  {
    # Blank line; do nothing.
  }
  else
  {
    die "Error: could not parse line of $list_of_algorithms_file!\n";
  }
}
close (LIST);
print " [done]\n";


print "Updating Makefile file...";
open (TEMPLATE, "Makefile.tpl") or
  die "Error: could not open Makefile template file!\n";
open (OUT, ">Makefile") or
  die "Error: could not open Makefile output file!\n";

printf OUT "# Makefile automatically generated by $0\n" .
           "# in $date_string.\n\n";

while (<TEMPLATE>)
{
  if ($_ =~ /\%template_class\%/)
  {
    foreach my $index (0..$#list_of_function_codes)
    {
      printf OUT "      src/functions/%s.o \\\n",
      $list_of_function_class_names[$index];
    }
    foreach my $index (0..$#list_of_algorithm_codes)
    {
      printf OUT "      src/algorithms/%s.o \\\n",
      $list_of_algorithm_class_names[$index];
    }
  }
  elsif ($_ =~ /\%template_test\%/)
  {
    foreach my $index (0..$#list_of_function_codes)
    {
      printf OUT "      test/functions/%s" . "Test.o \\\n",
        $list_of_function_class_names[$index];
    }
    foreach my $index (0..$#list_of_algorithm_codes)
    {
      printf OUT "      test/algorithms/%s" . "Test.o \\\n",
        $list_of_algorithm_class_names[$index];
    }
  }
  else
  {
    print OUT $_;
  }
}
close (TEMPLATE);
close (OUT);
print " [done]\n";