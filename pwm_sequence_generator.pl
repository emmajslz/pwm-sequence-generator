#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min), qw(max);
use Data::Dumper;

# ------------------------------------------------------------------------------
# MAIN PROGRAM
# ------------------------------------------------------------------------------

# We get the PWM from the filepath we had on the command line
my @pwmfile = get_pwm_from_file();

# We get rid of the lines we don't need and create an array with ONLY the PWM
my @pwm = extract_pwm(\@pwmfile);

my @p0 = @{$pwm[0]};

# We create a second PWM reordered by the score of the most frequent nucleotide
# per position (without the P0 row)
my @pwm2 = reorder_by_frequency(\@pwm);

# We use our recursive function to obtain all of the sequences (1000 at most)
my @sequences_in_arrays = sequence_generator(\@pwm2);

# We sort our sequences by frequency
my @sorted_sequences_in_arrays = sort_sequences_by_frequency(\@sequences_in_arrays);

# We calculate the total number of sequence we generated
my $total_nb_sequences = scalar @sorted_sequences_in_arrays;

# Our sequences are now stored in arrays,
# We convert each sequence to a string
my @sequences = sequence_arrays_to_strings(\@sorted_sequences_in_arrays, $total_nb_sequences);

# We print our results (either to an output file or the standard output)
if ((scalar @ARGV) > 1){
  # If we indicated an output file in the command line, we deposit the results
  # in said file
  output_to_file(\@sequences, $total_nb_sequences);

} else {
  #If we didn't specify an output file, we print the results to the STDOUT
  print Dumper \@sequences;
  print "\nTotal number of sequences: ";
  print scalar @sequences . "\n";

}

# ------------------------------------------------------------------------------
# PROGRAM FUNCTIONS
# ------------------------------------------------------------------------------

sub get_pwm_from_file{

  # Function to save the contents of the PWM TRANSMAC file to an array
  # The filepath was passed as a command line argument

  # We get the filepath from the command line
  my $pwmfilepath = $ARGV[0];

  # We open the file
  # We print an error message to the STDERR and exit in case of error in opening
  open(PWM_FILE, "< $pwmfilepath") || do {
    print STDERR "ERROR: Cannot open file $pwmfilepath";
    exit(1);
  };

  # We save the contents of the file to an array called pwmfile using chomp
  chomp(my @pwmfile = <PWM_FILE>);

  close(PWM_FILE);

  return @pwmfile;

}
sub extract_pwm{

  # Function to create an array that only has the pwm matrix data (from line P0)
  # We pass as argument the array that has the whole of the file content

  # We declare the argument of the function (we pass the array as a reference)
  my @pwmfile = @{$_[0]};

  my $length_pwmfile = scalar @pwmfile;
  my @pwm;
  my $first_pwm_line = 0;
  my $last_pwm_line = $length_pwmfile - 2;

  # We loop through the pwmfile array until we find the P0 line and we know
  # in which position it is. Once we reach that position, we then split every
  # line and add it to the new array
  for (my $i = 0; $i < $last_pwm_line; $i++) {
    if (substr($pwmfile[$i], 0, 2) eq "P0"){
      $first_pwm_line = 1;
    };
    if ($first_pwm_line){
      my @split = split(/\s+/, $pwmfile[$i]);
      @split = @split[0..4];
      push @pwm, \@split;
    }
  };

  return @pwm;
}

sub sorting_method{

  # Function to define the sorting method for the sort function.
  # This function will be used when we sort our PWM by nucleotide frequency

  return (max($b->[1], $b->[2], $b->[3], $b->[4]) <=>
          max($a->[1], $a->[2], $a->[3], $a->[4]));

}
sub reorder_by_frequency{

  # Function to create a second pwm that where we reorder by nucleotide
  # frequency.
  # We pass an argument with the pwm.

  # We declare the argument of the function (we pass the pwm as a reference)
  my @pwm = @{$_[0]};

  shift @pwm;

  my @pwm2 = sort sorting_method @pwm;

  return @pwm2;
}

sub get_max_frequency{

  # Function to get the highest frequency in a row, and the index it's stored in
  # The function requires the current row as argument
  my @row = @{$_[0]};

  my $max_frequency = 0;
  my $max_frequency_i;

  # We loop through the row to determine the maximum frequency, and store the
  # index it's found on.
  for (my $i = 1; $i < 5; $i++){
    if ($row[$i] > $max_frequency){
      $max_frequency = $row[$i];
      $max_frequency_i = $i;
    }
  }

  return $max_frequency, $max_frequency_i;
}
sub get_nucleotide{

  # Function that returns the nucleotide in the position $index (on the P0 row)
  my $index = $_[0];
  my $nucleotide = "";

  $nucleotide = $p0[$index];

  return $nucleotide;

}
sub row_only_frequency{

  # Function that will return a boolean indicating if the current row only has
  # one nucleotide with the maximum frequency (only one nucleotide possible)
  my @current_row = @{$_[0]};
  my $current_frequency = $_[1];

  my $maximum_frequency = 1;

  for (my $i = 1; $i < 5; $i++){
    if ($current_row[$i] != $current_frequency and $current_row[$i] != 0){
      $maximum_frequency = 0;
    }
  }

  return $maximum_frequency;
}

sub count_frequencies_in_row{

  # Function that will count how many nucleotides have a frequency in a row
  # (how many positions in the row are > 0)

  # Function used in get_sequence_count
  my @current_row = @{$_[0]};
  my $count = 0;

  for (my $i = 1; $i < 5; $i++){
    if (($current_row[$i] - 0) > 0){
      $count++;
    }
  }

  return $count;
}
sub get_sequence_count{

  # Function that will determine how many sequences we will have after the current
  # row on the pwm2 array.

  # This function will be used to determine wether proceeding with the following
  # row in the recursive function will make us have more sequences than allowed
  # by our cut-off value

  my @pwm2 = @{$_[0]};
  my @current_row = @{$_[1]};

  my @first_row = @{$pwm2[0]};
  my $sequence_count = count_frequencies_in_row(\@first_row);

  my $frequencies_in_row;
  my @row;

  for (my $i = 1; $i < (scalar @pwm2); $i++){

    @row = @{$pwm2[$i]};

    $frequencies_in_row = count_frequencies_in_row(\@row);
    $sequence_count = $sequence_count * $frequencies_in_row;

    if ((join ", ", @row) eq (join ", ", @current_row)){
      return $sequence_count;
    }
  }
  return $sequence_count;

}

sub sequence_generator{

  # Function that will initialize our array with the resulting sequences and
  # will call recursive_generator, the recursive function that will get us the
  # sequences
  my @pwm2 = @{$_[0]};

  my @current_row = @{$pwm2[0]};

  my $current_frequency;
  my $max_frequency_i;
  ($current_frequency, $max_frequency_i) = get_max_frequency(\@current_row);

  my @sequences;
  $sequences[0][0][$current_row[0]-1] = get_nucleotide($max_frequency_i);
  $sequences[0][1] = $current_frequency;

  my $max_frequency = $current_frequency;

  if (!row_only_frequency(\@current_row, $current_frequency)){
    for (my $i = 1; $i < 5; $i++){
      if ($i != $max_frequency_i and $current_row[$i] != 0){

        $current_frequency = $current_row[$i];

        my @new_sequence;

        @{$new_sequence[0][0]} = @{$sequences[0][0]};

        $new_sequence[0][0][$current_row[0]-1] = get_nucleotide($i);
        $new_sequence[0][1] = $current_frequency;

        push @sequences, recursive_generator(\@pwm2, \@new_sequence, \@pwm2, 0);
      }
    }
  }

  return recursive_generator(\@pwm2, \@sequences, \@pwm2, 0);

}
sub recursive_generator{

  # Recursive function that will generate all possible sequences from the PWM
  # (or only 1000 if there are more possible sequenes)

  my @pwm2 = @{$_[0]};
  my @sequences = @{$_[1]};
  my @og_pwm2 = @{$_[2]};
  my $limit_hit = $_[3];

  shift @pwm2;

  if(@pwm2){

    my @current_row = @{$pwm2[0]};

    my $current_frequency;
    my $max_frequency_i;
    ($current_frequency, $max_frequency_i) = get_max_frequency(\@current_row);

    my $old_frequency = $sequences[0][1];

    $sequences[0][0][$current_row[0]-1] = get_nucleotide($max_frequency_i);
    $sequences[0][1] = $old_frequency + $current_frequency;

    if (row_only_frequency(\@current_row, $current_frequency)){

      return recursive_generator(\@pwm2, \@sequences, \@og_pwm2, $limit_hit);

    } else {

      if (!$limit_hit and get_sequence_count(\@og_pwm2, \@current_row) < 1000){

        my $max_frequency = $current_frequency;

        for (my $i = 1; $i < 5; $i++){
          if ($i != $max_frequency_i and $current_row[$i] != 0){

            $current_frequency = $current_row[$i];

            my @new_sequence;

            @{$new_sequence[0][0]} = @{$sequences[0][0]};

            $new_sequence[0][0][$current_row[0]-1] = get_nucleotide($i);
            $new_sequence[0][1] = $old_frequency + $current_frequency;

            push @sequences, recursive_generator(\@pwm2, \@new_sequence, \@og_pwm2, $limit_hit);
          }
        }
      } else {

        # If the next row we'll do will make us surpass the cut-off value, we'll
        # go through the recursive function ONE more time, to surpass the
        # cut-off value only once, and then be able to cut at 1000 sequences.

        # If we don't do it this way, we'll possibly not hit 1000 exact sequences

        if (!$limit_hit){
          my $max_frequency = $current_frequency;

          for (my $i = 1; $i < 5; $i++){
            if ($i != $max_frequency_i and $current_row[$i] != 0){

              $current_frequency = $current_row[$i];

              my @new_sequence;

              @{$new_sequence[0][0]} = @{$sequences[0][0]};

              $new_sequence[0][0][$current_row[0]-1] = get_nucleotide($i);
              $new_sequence[0][1] = $old_frequency + $current_frequency;

              push @sequences, recursive_generator(\@pwm2, \@new_sequence, \@og_pwm2, 1);
            }
          }
        }

        $limit_hit = 1;

      }

      return recursive_generator(\@pwm2, \@sequences, \@og_pwm2, $limit_hit);
    }

  } else{

    return @sequences;
  }

}

sub sorting_method_sequences_by_frequency{

  return ($b->[1] <=> $a->[1]);
}
sub sort_sequences_by_frequency{

  # Function to sort our array with resulting sequences by frequency (highest to lowest)
  my @sequences = @{$_[0]};

  my @sorted_sequences = sort sorting_method_sequences_by_frequency @sequences;

  if ((scalar @sorted_sequences) > 1000){
    @sorted_sequences = @sorted_sequences[0..999];
  }

  return @sorted_sequences;

}

sub sequence_arrays_to_strings{

  my @sequences_in_arrays = @{$_[0]};
  my $total_nb_sequences = $_[1];

  my @sequences;

  for (my $i = 0; $i < $total_nb_sequences; $i++){
    $sequences[$i][0] = join("", @{$sequences_in_arrays[$i][0]});
    $sequences[$i][1] = $sequences_in_arrays[$i][1];
  }

  return @sequences;

}

sub output_to_file{

  my @sequences = @{$_[0]};
  my $total_nb_sequences = $_[1];

  # We get the filepath from the command line
  my $outputfilepath = $ARGV[1];

  # We open the file in write mode
  # We print an error message to the STDERR and exit in case of error in opening
  open(OUTPUT_FILE, "> $outputfilepath") || do {
    print STDERR "ERROR: Cannot open file $outputfilepath";
    exit(1);
  };

  my $pwmfilepath = $ARGV[0];

  print OUTPUT_FILE "Sequences generated from PWM matrix in $pwmfilepath\n\n";

  print OUTPUT_FILE ("Total number of sequences: $total_nb_sequences\n\n");
  print OUTPUT_FILE ("Sequence\t\tFrequency\n");

  for (my $i = 0; $i < $total_nb_sequences; $i++){
    print OUTPUT_FILE ("$sequences[$i][0]\t$sequences[$i][1]\n");
  }

}
