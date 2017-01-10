
#
# run_benchmarking.pl : program to create a benchmarking table among
#     different algorithms and cost functions. In this table, for a given
#     cost function, the rows represent different instance sizes, while
#     the columns are the performance of different algorithms, both in terms
#     of required computational time and number of calls of the cost function.
#
#    This file is part of the featsel program
#    Copyright (C) 2016  Marcelo S. Reis
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use POSIX qw(ceil floor);
use strict;
use Time::HiRes qw (gettimeofday tv_interval);

use lib './lib';

# %template_use%


# Constant that works as an upper bound limit for cost function values.
#
my $INFINITY = 1000000;


# Output file prefix
#
my $output_file_prefix = "foo";


# Flag that verifies whether the output graphs should be generated or not.
#
my $PRINT_OUTPUT_GRAPHS = "ON"; # "ON" to print output graphs, "OFF" otherwise


# Directory where the input instances are stored.
#
my $INPUT_DIR = "input/tmp";


# Directory where the output files are stored.
#
my $OUTPUT_DIR = "output";


# Name of the log file used when the featsel binary is executed.
#
my $LOG_FILE = $OUTPUT_DIR . "/result.log";


# featsel binary file.
#
my $FEATSEL_BIN = "./bin/featsel";


#------------------------------------------------------------------------------#
#
# Process the input parameters.
#

# Instance mode:
#
# 0: for each instance size in [1,n], creates $k random instances of the 
# specified cost function in the default directory "input/tmp"
#
# 1: read input files from the directory established in $input_directory
#
my $instance_mode = 0;      


# Employed cost function code.
#
my $cost_function = "bar";                


# Number of instances to be processed for each instance size (integer >= 1)
#
my $number_of_instances_per_size = 1;


# Maximum instance size (integer >= 1)
#
my $maximum_instance_size = 1;


# Search mode:
#
# 0: complete search (may be optimal or not, depending on the algorithm) 
#
# 1: search constrained by a given number of cost function calls
#
my $search_mode = 0;


# Cost function calls limit for the search mode 1
#
my $COST_FUNCTION_LIMIT = 2 ** $maximum_instance_size;


if ((@ARGV == 6) || (@ARGV == 7))
{
  $output_file_prefix = $ARGV[0];

  $instance_mode = $ARGV[1];

  $cost_function = $ARGV[2];
  
  $number_of_instances_per_size = $ARGV[3];

  $maximum_instance_size = $ARGV[4];

  $search_mode = $ARGV[5];

  if (!(($search_mode == 0) || ($search_mode == 1)))
  {
    die "Error: search mode must be either 0 or 1!\n";
  }
  elsif (($search_mode == 1) && (@ARGV < 6))
  {
    die "Error: search mode $search_mode requires a threshold value!\n";
  }
  (@ARGV == 7) and $COST_FUNCTION_LIMIT = $ARGV[6];
}
else
{
  die "\nSyntax: $0 " . 
      " OUTPUT_FILE_PREFIX\n        INSTANCE_MODE  COST_FUNCTION_CODE k  n   " . 
      "SEARCH_MODE  [max_number_of_calls]\n\n" .
      
      "Where:\n\n" .

      "    OUTPUT_FILE_PREFIX: prefix of the output file name.\n\n" .

      "    COST_FUNCTION_CODE: code of the employed cost function.\n\n" .
      
      "    INSTANCE_MODE: must be 0 or 1; if the value is 0, then for each\n" .
      "    size of instance in [1,n], it creates k random instances of the\n" .
      "    chosen cost function. Otherwise, it reads input files already\n".
      "    generated, which must have k files per instance size, which in\n" .
      "    turn must range from 1 to n.\n\n" .

      "    SEARCH_MODE: must be 0 or 1; if the value is 0, then it performs\n" .
      "    a complete search, which may be optimal or not depending on the\n" .
      "    algorithm. Otherwise, it runs a search constrained by\n" .
      "    'max_number_of_calls' calls of the cost function.\n\n"
      ;
}


#------------------------------------------------------------------------------#
#
# Loading the list of available cost functions, to verify whether the chosen
# cost function is available.
#
print "Starting $0 program.\nExecution on instances from size 1 to ". 
      "$maximum_instance_size.\nFor each size, it calculates the average of " . 
      "$number_of_instances_per_size executions.\n";

(defined $COST_FUNCTION_LIMIT) 
  and print "Maximum number of calls of the cost function was set as " . 
      $COST_FUNCTION_LIMIT . "\n";

print "\nVerifying the availability of cost function '$cost_function'...";

my $list_of_cost_functions_file = "ListOfCostFunctions.txt";

open (LIST, $list_of_cost_functions_file) or 
  die "Error: could not open $list_of_cost_functions_file for reading!\n";

my $cost_function_name;
my $cost_function_file_type;

while (<LIST>)
{
  chomp $_;
  # $1 == code, $2 == class name, $3 == file type
  #
  if ($_ =~ /^\s*(\w+)\s+(\w+)\s+(\w+)\s*$/) 
  {
    $cost_function eq $1
      and $cost_function_name = $2
      and $cost_function_file_type = $3;
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

defined $cost_function_name 
  or die "Error: cost function '$cost_function' is not available!\n\n";

print " [done]\n";

#------------------------------------------------------------------------------#
#
# Loading the list of available algorithms.
#

print "\nLoading the list of available algorithms...";

my $list_of_algorithms_file = "ListOfAlgorithms.txt";

open (LIST, $list_of_algorithms_file) or 
  die "Error: could not open $list_of_algorithms_file for reading!\n";

my @list_of_algorithm_codes;
my %list_of_algorithm_class_names;

while (<LIST>)
{
  chomp $_;
  if ($_ =~ /^\s*(\w+)\s+(\w+)\s*$/) # $1 == code, $2 == class name
  {  
    push @list_of_algorithm_codes, uc $1;
    $list_of_algorithm_class_names{uc $1} = $2; 
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


#------------------------------------------------------------------------------#
#
# At this point the user can choose a subset of algorithms, through removal 
# of codes from the @algorithms array.
#
my @algorithms = @list_of_algorithm_codes;

# Number of algorithms.
#
my $number_of_algorithms = @algorithms;


open (OUT, ">output/" . $output_file_prefix . "_table.html") 
  or die "Error: could not open HTML output file!\n";

#------------------------------------------------------------------------------#
#
# Printing HTML output file header
#
print OUT "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";

print OUT "<HTML><HEAD>" .
"<meta http-equiv=\"content-type\" content=\"text/html; charset=iso-8859-1\">" .
  "<TITLE>Feature selection algorithms and cost functions " . 
  "benchmarking</TITLE></HEAD>\n";

print OUT "<BODY><H1>Feature selection algorithms and cost functions " .
  "benchmarking</H1><HR>\n";

print OUT "<H3>Table description:</H3>\n";

print OUT "<EM>2^|S|:</EM> size of the search space of an instance of S " .
  "elements.<BR>\n";

print OUT "<EM>Total Time (sec):</EM> required time to run each algorithm " .
  "(average of up to $number_of_instances_per_size runnings).<BR>\n";

print OUT "<EM>Cost Function Time (sec):</EM> required time for the " .
  "computation of all calls the cost function, during the running of each " .
  "algorithm (average of $number_of_instances_per_size runnings).<BR>\n";

print OUT "<EM>\# Computed nodes:</EM> number of times the chosen " . 
  "cost function ($cost_function_name) " .
  "is computed by each algorithm (average of $number_of_instances_per_size " . 
  "runnings).<BR>\n";

print OUT "<EM>\# The best solution:</EM> number of times " .
  "(out of $number_of_instances_per_size) that each " .
  "algorithm had the best solution among all algorithms.<BR><HR><BR>\n";

print OUT "<TABLE border='5'>\n";

print OUT "<TR bgcolor='yellow'>\n";

print OUT "<TD colspan=2><CENTER>Size of instance</CENTER></TD>\n";

print OUT "<TD><CENTER>&nbsp;</CENTER></TD>\n" .
  "<TD colspan=$number_of_algorithms><CENTER>Total Time (sec)</CENTER></TD>\n" .
  "<TD><CENTER>&nbsp;</CENTER></TD>\n" .
  "<TD colspan=$number_of_algorithms><CENTER>Cost Function Time (sec)</CENTER>" .
  "<TD><CENTER>&nbsp;</CENTER></TD>\n" .
  "<TD colspan=$number_of_algorithms><CENTER>\# Computed nodes</CENTER></TD>\n" .
  "<TD><CENTER>&nbsp;</CENTER></TD>\n" .
  "<TD colspan=$number_of_algorithms><CENTER>\# The best solution</CENTER>" .
  "</TD>\n</TR>\n";

print OUT "<TR bgcolor='yellow'>\n";

# Instance size
#
print OUT "<TD><CENTER>&nbsp;|S|&nbsp;</CENTER></TD><TD><CENTER>&nbsp;2^|S|" . 
  "&nbsp;</CENTER></TD>";

# Total time
#
print OUT  "<TD><CENTER>&nbsp;</CENTER></TD>";
for (my $i = 0; $i < $number_of_algorithms; $i++)
  {
    printf OUT "<TD><CENTER>&nbsp; %s &nbsp;</CENTER></TD>", $algorithms[$i];
  }

# Cost function time
#
print OUT  "<TD><CENTER>&nbsp;</CENTER></TD>";
for (my $i = 0; $i < $number_of_algorithms; $i++)
  {
    printf OUT "<TD><CENTER>&nbsp; %s &nbsp;</CENTER></TD>", $algorithms[$i];
  }

# Number of computed nodes
#
print OUT  "<TD><CENTER>&nbsp;</CENTER></TD>";
for (my $i = 0; $i < $number_of_algorithms; $i++)
  {
    printf OUT "<TD><CENTER>&nbsp; %s &nbsp;</CENTER></TD>", $algorithms[$i];
  }

# Number of best solutions found
#
print OUT  "<TD><CENTER>&nbsp;</CENTER></TD>";
for (my $i = 0; $i < $number_of_algorithms; $i++)
  {
    printf OUT "<TD><CENTER>&nbsp; %s &nbsp;</CENTER></TD>", $algorithms[$i];
  }
printf OUT "</TR>\n";


#------------------------------------------------------------------------------#

# Create random instances (if is the case) and load the input files.
#

my @instance_file;             # Store the list of instance file names.

my $MAX_ELEM_VALUE = 100000;   # Maximum value of an element of S.

if ($instance_mode == 0)
{
  foreach my $i (1..$maximum_instance_size)    
  {
    # Create files (either .dat or .xml) containing random instances.
    #
    foreach (1..$number_of_instances_per_size)
    {
      if (0)
      {
        # NOP.
      }
# %template_if%

      elsif ($cost_function eq "mce")
      {
        MeanConditionalEntropy::random_mce_instance
        ($i, $MAX_ELEM_VALUE, $INPUT_DIR, $_);
      }

      else
      {
        die "The '$cost_function' instance generator is not implemented " .
            "yet! To this end, a subroutine must be included into this code!\n";
      }
    }
  }
}

# Load file names.
#
opendir (my $dh, $INPUT_DIR) or die "Cannot open input directory: $!\n";
@instance_file = grep { /.*\.$cost_function_file_type$/ && -f "$INPUT_DIR/$_" }
  readdir ($dh);
closedir $dh;

my %experiment;


(scalar(@instance_file) >= $number_of_instances_per_size*$maximum_instance_size)
  or die "Insufficient number of instance files stored in '$INPUT_DIR'!\n\n";

foreach my $file (sort @instance_file)
{
  if ($file =~ /Test_0*(\d+)_\w+/)
  {
    push @{$experiment{$1}}, $file;
  }
  else
  {
    die "\nError: '$file' file name must follow this format: Test_\\d+_\\w+." .
        $cost_function_file_type . ".\n\n";
  }
}


#------------------------------------------------------------------------------#

# Run experiments
#
my %total_time;
my %cost_function_time;
my %cost_function_calls;
my $max_number_of_calls = " ";

# This variable stores the maximum required time among all algoritms; this is
# useful in order to generate the graphs later.
#
my $MAX_TIME_VALUE = 0;

print "\nRunning benchmarking experiments with $number_of_algorithms " .
      "algorithms and instances of size up to $maximum_instance_size.\n\n";

foreach my $i (1..$maximum_instance_size)    
{
  print "Starting iteration $i... ";

  # Current instance size (|S| = $i, implying on a 2^$i search space)
  #
  printf OUT "<TR><TD><CENTER>&nbsp;%2d&nbsp;</CENTER></TD>", $i;
  printf OUT  "<TD><CENTER>&nbsp;%7d&nbsp;</CENTER></TD>", 2 ** $i;

  # This informaton is extracted from all algorithms
  #
  my @average_calls_of_cost_function;
  my @average_time_of_cost_function;
  my @average_time_of_algorithm;
  my @number_of_times_that_has_a_best_solution;

  for (my $j = 0; $j < $number_of_algorithms; $j++)
  { 
    $average_calls_of_cost_function[$j] = 0;
    $number_of_times_that_has_a_best_solution[$j] = 0;
  }

  # If the search mode is one, then the procedure has an upper limit
  # for the number of calls of the cost function.
  #
  if ($search_mode == 1)
  {
    $max_number_of_calls = " -t $COST_FUNCTION_LIMIT ";
  } 

  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    $average_calls_of_cost_function[$j] = 0;
    $average_time_of_cost_function[$j] = 0;
    $average_time_of_algorithm[$j] = 0;
    $number_of_times_that_has_a_best_solution[$j] = 0;
  }
    
  foreach my $k (1..$number_of_instances_per_size)
  {
    my $best_solution = $INFINITY;
    my @minimum_of_algorithms;

    for (my $j = 0; $j < $number_of_algorithms; $j++)
    {
      my ($t0, $t1);
      my $current_algorithm = lc $algorithms[$j];

      $t0 = [gettimeofday];
      system ("$FEATSEL_BIN -n $i -a $current_algorithm " . 
              "-c $cost_function -f $INPUT_DIR/" . $experiment{$i}->[$k-1]  . 
              $max_number_of_calls . " > $LOG_FILE");
      $t1 = [gettimeofday];

      $average_time_of_algorithm[$j] += tv_interval ($t0, $t1);

      open (LOG, $LOG_FILE);
      while (<LOG>)
      {
        if ($_ =~ /(\<\d+\>)\s+\:\s+(\S+)/)
        {
          $minimum_of_algorithms[$j] = $2;
        }
        elsif ($_ =~ /^Number\s+of\s+visited\s+subsets\:\s+(\S+)/)
        { 
          $average_calls_of_cost_function[$j] += $1;
        }
        elsif ($_ =~ /subsets\:\s+(\d+)\s+microseconds/)
        {
          $average_time_of_cost_function[$j] += $1;
        }
      }
      close(LOG);

      if ($best_solution > $minimum_of_algorithms[$j])
      {
        $best_solution = $minimum_of_algorithms[$j];
      }

    } # for (my $j = 0; $j < $number_of_algorithms; $j++)
      
    for (my $j = 0; $j < $number_of_algorithms; $j++)
    {
      if ($minimum_of_algorithms[$j] == $best_solution)
      {
 	      $number_of_times_that_has_a_best_solution[$j]++;
      }
    }
      
  } # foreach my $k (1..$number_of_instances_per_size)
  
  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    $average_calls_of_cost_function[$j] /= $number_of_instances_per_size;	 

    $cost_function_calls{$algorithms[$j]}->[$i] = 
      $average_calls_of_cost_function[$j];      

    $average_time_of_cost_function[$j] /= 
      ($number_of_instances_per_size * 1000000); # convert to seconds
      
    $average_time_of_algorithm[$j] /= $number_of_instances_per_size;

    $total_time{$algorithms[$j]}->[$i] = $average_time_of_algorithm[$j];
    
    if ($average_time_of_algorithm[$j] > $MAX_TIME_VALUE){    
	     $MAX_TIME_VALUE = $average_time_of_algorithm[$j];
      }
    
    $cost_function_time{$algorithms[$j]}->[$i] = 
      $average_time_of_cost_function[$j];      
  }
  
  # Total time output
  #
  printf OUT "<TD><CENTER>&nbsp;</CENTER></TD>\n";
  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    printf OUT "<TD><CENTER>&nbsp;%4.2f&nbsp;</CENTER></TD>",
      $average_time_of_algorithm[$j];
  }
    
  # Cost function time output
  #
  printf OUT "<TD><CENTER>&nbsp;</CENTER></TD>\n";
  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    printf OUT "<TD><CENTER>&nbsp;%4.2f&nbsp;</CENTER></TD>",
      $average_time_of_cost_function[$j];
  }
   
  # Number of computed nodes output
  #
  printf OUT "<TD><CENTER>&nbsp;</CENTER></TD>\n";
  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    printf OUT "<TD><CENTER>&nbsp;%5.2f&nbsp;</CENTER></TD>",
      $average_calls_of_cost_function[$j];
  }

  # Number of best solutions found output
  #
  printf OUT "<TD><CENTER>&nbsp;</CENTER></TD>\n";
  for (my $j = 0; $j < $number_of_algorithms; $j++)
  {
    printf OUT "<TD><CENTER>&nbsp;%d&nbsp;</CENTER></TD>",
      $number_of_times_that_has_a_best_solution[$j];
  }
  print OUT "</TR>\n";

  print "[done]\n";

} # foreach my $i (1..$maximum_instance_size)


# Create the table and HTML footers.
#

print OUT "</TABLE><BR><HR>\n";
print OUT "
<p>
    <a href=\"http://validator.w3.org/check?uri=referer\"><img
        src=\"http://www.w3.org/Icons/valid-html401-blue\"
        alt=\"Valid HTML 4.01 Strict\" height=\"31\" width=\"88\"></a>
</p>\n";
print OUT "</BODY></HTML>\n";

close (OUT);


# Print the output graphs.
#
if ($PRINT_OUTPUT_GRAPHS eq "ON")
{
  print_time_graphs ();
  
  # TODO: the user can include here new subs for printing graphs. 
}


# Remove temporary files.
#
system ("rm $LOG_FILE");


# End of main program.
#
print "\nEnd of execution.\n\n";

exit 0;


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# TODO: the user can add here new subs for printing other type of graphs.
#

#------------------------------------------------------------------------------#

# This sub prints, for each algorithm, a graphic depicting how the required 
# computational time evolves according to the instance size.
#
sub print_time_graphs
{
  my $GNUPLOT_DATA_FILE = $OUTPUT_DIR . "/dat.temp";

  my $GNUPLOT_PLOT_FILE = $OUTPUT_DIR . "/gnuplot.temp";

  print "\nPrinting output graphs... ";

  foreach my $algo (@algorithms)
  {
    open (DATA, ">$GNUPLOT_DATA_FILE");

    for (my $i = 1; $i <= $maximum_instance_size; $i++)
    {
      printf DATA "$i %.4f %.4f\n", $cost_function_calls{$algo}->[$i],
                                    $total_time{$algo}->[$i];
    }

    close (DATA);

    my $Xaxis = 75 * $maximum_instance_size;

    my $grid_points = $maximum_instance_size * 3;

    open (PLOT, ">$GNUPLOT_PLOT_FILE");

    printf PLOT "set terminal png enhanced crop size $Xaxis, $Xaxis\n";
    printf PLOT "set output '$OUTPUT_DIR/$output_file_prefix\_$algo.png'\n";
    printf PLOT "unset key\n";
    printf PLOT "set xlabel \"|S|\" rotate parallel\n";
    printf PLOT "set ylabel \"Number of computed nodes\" rotate parallel\n";
    printf PLOT "set zlabel \"Total time (seconds)\" rotate parallel\n";
    printf PLOT "unset colorbox\n";
    printf PLOT "set dgrid3d $grid_points, $maximum_instance_size\n";
    printf PLOT "set hidden3d\n";  
    printf PLOT "splot \"$GNUPLOT_DATA_FILE\" using 1:2:3 with lines\n";

    # Execute Gnuplot.
    #
    system ("gnuplot $GNUPLOT_PLOT_FILE");

    # Remove temporary files.
    #
    system ("rm $GNUPLOT_DATA_FILE");
    system ("rm $GNUPLOT_PLOT_FILE");
  }
  print "[done]\n";
}


#------------------------------------------------------------------------------#

#
# End of file.

