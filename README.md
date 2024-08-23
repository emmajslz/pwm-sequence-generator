# PWM Sequence Generator - PERL Script

This report explains the goal and functioning of `pwm_sequence_generator.pl`

> *“Position Weight Matrices are a simple way to model signals appearing on DNA and protein sequences. They summaryze the frequencies of every letter of the nucleotide or amino acid alphabets at a given position of the signal, for instance a splice site or a transcription factor binding site. This means that we can obtain a score, either from the absolute frequencies or by transforming them into log-likelihoods, for a given sequence to determine if it contains the signal pattern or not. Yet, they can also be used as sequence generators.”*
> 

Our goal is to create PERL script that will take as input a Position Weight Matrix and will generate all of the possible sequences from this PWM.

Some guidelines we need to take into account:

- We will only work with nucleotides
- We will set up a cut-off at 1000 sequences
    - Depending on the way our PWM is set up, and being the complexity of this program exponential, we could end up with an execution time too long to handle. Thus the reason we need a cut-off value.
- Our program will use a recursive function

---

# Calling the script (command line)

To execute our script, we have to call it using the following specified syntax:

`perl pwm_sequence_generator.pl pwm_filepath [output_filepath]`

**The arguments of our script are the following:**

- `pwm_filepath` → the filepath where out PWM in TRANSFAC format is stored.
- `[output_filepath]` (*optional)*
    - If given → The resulting sequences will be written into the file at `output_filepath`. If the file doesn’t exist, it’ll be created
    - If NOT given → The resulting sequences will be printed to the STDOUT.

# Script description

## Main pipeline

The main steps the script follows are the following:

1. We read a file containing the PWM in TRANSFAC format.
2. We transfer all of the information into a two dimensional array → `@pwm`
3. We create a new array, sorting this array (`@pwm`) by the score of the most frequent nucleotide per position → `@pwm2`
4. We generate all of the sequences recursively, until we have looped through all of `@pwm2` or until we hit the cut-off value.
    1. We now have a three-dimensional array with the sequences (as arrays) and the frequencies.
5. We convert the sequences into strings
6. *If* we specified an output file, we write the results to said file. If not, we print them to the STDOUT.

## Detailed description of the code

The code is structured in **functions**. Because of this, we have firstly, shebang, pragmas and library declarations, we then have the main program, in which we call all of our functions, and then we have the functions.

### Shebang, pragmas and libraries

```perl
#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min), qw(max);
use Data::Dumper;
```

We declare the shebang, the strict and warnings pragmas to restrict our code, and we include the following libraries:

- `List::Until` → Methods `qw(min)` and `qw(max)`
    - Used to indentify the minimum and maximum values in a list of values.
- `Data::Dumper`
    - Useful for data visualization (very useful in this case to correctly visualize multi-dimensional arrays.

### Main program

The main program follows the following pipeline:

1. We get the pwm file
2. We extract the information from the file and deposit it into an array
3. We save the P0 row, containing the nucleotides (used later to know which nucleotide is in each position)
4. We reorder the pwm array by the score of the most frequent nucleotide per position (we don’t include the P0 row)
5. We use the `sequence_generator` function to create our array with the sequences (this function will initialize the array and then call the recursive function that will build our sequence array)
6. We obtain the total number of sequence we generated
7. We convert our sequences (saved in arrays) into strings
8. *If* we specified an output file, we write the results to said file. If not, we print them to the STDOUT.

### Functions

A detailed description of each function:

<aside>
⚠️ When we refer to “current row” it means the row we are currently on as the recursive function runs and goes through the pwm2 array

</aside>

- **`get_pwm_from_file`**
    - Function to save the contents of the PWM TRANSMAC file to an array
    - The filepath was passed as a command line argument
    
    ---
    
    1. We obtain the filepath from the command line arguments
    2. We open the file and print an error message to the STDERR and exit in case of error in opening
    3. We save the contents of the file to an array called pwmfile using chomp
    4. We close the file
    
    ---
    
    We return `@pwmfile`
    

- **`extract_pwm`**
    - Function to create an array that only has the pwm matrix data (from line P0)
    - **Arguments:**
        - `my @pwmfile = @{$_[0]};`
            - The PWM file contents in array format
    
    ---
    
    1. We loop through the pwmfile array (using a for loop) until we find the P0 line and we know in which position it is. Once we reach that position, we then split every line and add it to the new array
    
    ---
    
    We return `@pwm`
    

- **`sorting_method`**
    - Function to define the sorting method for the sort function.
    - This function will be used when we sort our PWM by nucleotide frequency

- **`reorder_by_frequency`**
    - Function to create a second pwm that where we reorder by nucleotide frequency
    - **Arguments:**
        - `my @pwm = @{$_[0]};`
            - The PWM matrix in our two dimensional array
    
    ---
    
    1. We erase the first row using shift (the P0 row).
    2. We use the sort function (with our predefined sorting method) to reorder our matrix
    
    ---
    
    We return `@pwm2`
    

- **`get_max_frequency`**
    - Function to get the highest frequency in a row, and the index it's stored in
    - **Arguments:**
        - `my @row = @{$_[0]};`
            - The row of the PWM matrix from wich we want to obtain the maximum frequency
    
    ---
    
    1. We initialize our variables
    2. We loop through the row to determine the maximum frequency, and store the index it's found on.
    
    ---
    
    We return `@max_frequency` and `$max_frequency_i` (the index where max frequency is stored).
    

- **`get_nucleotide`**
    - Function that returns the nucleotide in the position $index (on the P0 row)
    - **Arguments:**
        - `my $index = $_[0];`
            - The index of the nucleotide we want to obtain
    
    ---
    
    1. We simply access the `$p0` array (P0 row) at position `$index` to get the character corresponding to our desired nucleotide.
    
    ---
    
    We return `@nucleotide`
    

- **`row_only_frequency`**
    - Function that will return a boolean indicating if the current row only has one nucleotide with the maximum frequency (only one nucleotide possible)
    - **Arguments:**
        - `my @current_row = @{$_[0]};`
            - The current row
        - `my $current_frequency = $_[1];`
            - The current frequency, meaning the maximum frequency on the row that was previously obtained (with the `get_max_frequency function`)
    
    ---
    
    1. We initialize the “boolean” variable `$maximum_frequency` at True (1).
    2. We loop through the row
        1. If we encounter a frequency that isn’t the current frequency, we change `$maximum_frequency` to False (0).
    
    ---
    
    We return `@maximum_frequency`
    

- **`count_frequencies_in_row`**
    - Function that will count how many nucleotides have a frequency in a row (how many positions in the row are > 0)
    - Function used in `get_sequence_count`
    - **Arguments:**
        - `my @current_row = @{$_[0]};`
            - The current row
    
    ---
    
    1. We initialize `$count` at 0
    2. We loop through the row
        1. If we encounter a frequency that is greater than 0, we increment `$count`.
    
    ---
    
    We return `@count`
    

- **`get_sequence_count`**
    - Function that will determine how many sequences we will have after the current row on the pwm2 array.
    - This function will be used to determine wether proceeding with the following row in the recursive function will make us have more sequences than allowed by our cut-off value
    - **Arguments:**
        - `my @pwm2 = @{$_[0]};`
            - The pwm2 array
        - `my @current_row = @{$_[1]};`
            - The current row
    
    ---
    
    1. We create an array with the first row of the pwm2 array.
    2. We initialize `$sequence_count` using the function `count_frequencies_in_row`. We pass as argument the first row to obtain the number of possible nucleotides in the first row (so the number of sequences we will start with).
    3. We initialize the variables `$frequencies_in_row` and `@row` that we’ll use in every iteration of the for loop.
    4. We loop through the pwm2 array. In each iteration:
        1. We get the current row
        2. We get the number of possible nucleotides in said row
        3. We multiply this number by the sequence count we already have.
        4. If we reach the current row, we return the sequence count
    
    ---
    
    We return `@sequence_count`
    

- **`sequence_generator`**
    - Function that will initialize our array with the resulting sequences and will call recursive_generator, the recursive function that will get us the sequences
    - **Arguments:**
        - `my @pwm2 = @{$_[0]};`
            - The pwm2 array
    
    ---
    
    1. We set `@current_row` to the first row in the pwm2 array
    2. We get `$current_frequency` and `$max_frequency_i` (the maximum frequency in the current row and its position / index)
    3. We initialize the `@sequences` array, where we’ll store all of our resulting sequences and their frequencies
        1. This array is multi-dimensional. It’s should look like:
            1. `[[[nucleotide 1-1, nucleotide 1-2, ... nucleotide 1-n], frequency1], 
              [[nucleotide 2-1, nucleotide 2-2, ... nucleotide 2-n], frequency2], 
              ...                                                                 
              [[nucleotide m-1, nucleotide m-2, ... nucleotide m-n], frequency m]]`
    4. We initialize `$sequences[0][0][$current_row[0]-1]` with its corresponding nucleotide. We use the `get_nucleotide` function.
        1. `$current_row[0] - 1` will give us the sequence position where our nucleotide has to go. 
        2. The sequence position is saved as a string (`“01”` or `“02”` for example).
        3. We do `$current_row[0] - 1` to obtain its index (because indexes start with 0) as a numerical variable.
    5. We initialize `$sequences[0][1]` to the current frequency.
    6. If the first row has more than one frequency (we check with the `row_only_frequency` function), ee loop through the row. When we reach a position that isn’t the one we already initialized (the `$i` is different than `$max_frequency_i`):
        1. We change $current_frequency to the current frequency.
        2. We initialize an array `@new_sequence` where our new sequence will go.
        3. We initialize `$new_sequence[0][0][$current_row[0]-1]` with its corresponding nucleotide. We use the `get_nucleotide` function.
        4. We initialize `$new_sequence[0][1]` to the current frequency.
        5. We use push to add this new sequence to our `@sequences` array. We can’t simply push the array, we have to push the array that will result of calling the recursive function with this new sequence, to be able to complete it.
    7. We call the recursive function (`@recursive_generator`) when we return the sequence, thus starting the recursion and completing the `@sequences` array.

- **`recursive_generator`**
    - Recursive function that will generate all possible sequences from the PWM (or only 1000 if there are more possible sequenes)
    - **Arguments:**
        - `my @pwm2 = @{$_[0]};`
            - The pwm2 array
        - `my @sequences = @{$_[1]};`
            - The sequences array (where our resulting sequences go)
        - `my @og_pwm2 = @{$_[2]};`
            - The original pwm2. As progress (recursively), we will erase rows from the pwm2 array. We need to keep the original pwm2 array to calculate the number of sequences we’ll get each iteration (to respect the cut-off value).
        - `my $limit_hit = $_[3];`
            - Boolean that indicated wether or not we have hit the cut-off value
    
    ---
    
    1. We use shift to erase the first row in the pwm2 array. This way, each iteration, we’ll have the row we need to work on as the first row.
    2. If the pwm2 array is not empty:
        1. We set `@current_row` to the first row in the pwm2 array
        2. We get `$current_frequency` and `$max_frequency_i` (the maximum frequency in the current row and its position / index)
        3. We initialize the variable `$old_frequency` where we’ll store the current frequency (which we’ll need later to be able to add the current frequency and have the total frequency for each new sequence.)
        4. We assign to `$sequences[0][0][$current_row[0]-1]`  its corresponding nucleotide. We use the `get_nucleotide` function.
            1. `$current_row[0] - 1` will give us the sequence position where our nucleotide has to go. 
            2. The sequence position is saved as a string (`“01”` or `“02”` for example).
            3. We do `$current_row[0] - 1` to obtain its index (because indexes start with 0) as a numerical variable.
        5. We assign to `$sequences[0][1]` the current frequency: We add `$current_frequency` to `$old_frequency` to get the total frequency at this moment.
        6. If the row only has one frequency (we check with the `row_only_frequency` function), we simply call the function (`recursive_generator`) to start the recursion and continue the sequence.
        7. If the row has more than one frequency (we check with the `row_only_frequency` function), 
            1. If we **haven’t** reached the cut-off value, we loop through the row. When we reach a position that isn’t the one we already initialized (the `$i` is different than `$max_frequency_i`):
                1. We change $current_frequency to the current frequency.
                2. We initialize an array `@new_sequence` where our new sequence will go.
                3. We initialize `$new_sequence[0][0][$current_row[0]-1]` with its corresponding nucleotide. We use the `get_nucleotide` function.
                4. We initialize `$new_sequence[0][1]` to the current frequency.
                5. We use push to add this new sequence to our `@sequences` array. We can’t simply push the array, we have to push the array that will result of calling the recursive function with this new sequence, to be able to complete it.
            2. If we **have** reached the cut-off value, we do the process once again to ensure we will surpass the cut-off value, but we set the boolean `$limit_hit` to true (1), to not go deeper than that. This will allow for us to have 1000 sequences (if we don’t do another level of recursion, we probably won’t reach 100’ exact sequences).
        8. We call the recursive function (`@recursive_generator`) when we return the sequence, thus starting the recursion and completing the `@sequences` array.
    3. If pwm2 is empty, there’s no more sequence positions to go through, so we simply return the `@sequences` array and we finish the recursion

- **`sorting_method_sequences_by_frequency`**
    - Function to define the sorting method for the sort function.
    - This function will be used when we sort our resulting sequences by frequency in the `sort_sequences_by_frequency` function.

- **`sort_sequences_by_frequency`**
    - Function to sort our array with resulting sequences by frequency (highest to lowest)
    - **Arguments:**
        - `my @sequences = @{$_[0]};`
            - The sequences array with all of the resulting sequences.
    
    ---
    
    1. We create a new arrray with the sorted sequences using the sort function with the previously defined sorting method.
    2. If we have more than 1000 sequences, we only take the first 1000 (to keep only the highest sequences).
    
    ---
    
    We return `@sorted_sequences`
    

- **`sequence_arrays_to_strings`**
    - Function to convert our sequences (saved as arrays) to strings
    - **Arguments:**
        - `my @sequences_in_arrays = @{$_[0]};`
            - The sequences array with all of the resulting sequences (in arrays)
        - `my $total_nb_sequences = $_[1];`
            - The total number of sequences we have
    
    ---
    
    1. We initialize the array where we’ll store our sequence as strings
    2. We loop through the `@sequences_in_arrays` array and we use the join function to turn the arrays to strings.
    
    ---
    
    We return `@sequences`
    

- **`output_to_file`**
    - Function to write our results to the output file if it was given
    - **Arguments:**
        - `my @sequences = @{$_[0]};`
            - The sequences array with all of the resulting sequences
        - `my $total_nb_sequences = $_[1];`
            - The total number of sequences we have
    
    ---
    
    1. We obtain the filepath from the command line arguments
    2. We open the file (in write mode) and print an error message to the STDERR and exit in case of error in opening
    3. We print the header into the output file
    4. We print each line in our array to the output file
    5. We close the file

# Results

Our resulting file can look something like this:

---

Sequences generated from PWM matrix in /Users/emmajuansalazar/Arxius/MBI/T1/PER/Emma_JUAN_PERL_project/pwm1.txt

Total number of sequences: 384

Sequence		Frequency
GGACATGCCCGGGCATGTCC	307
GGACATGCCCGGGCATGTCT	307
GGACATGCCCGGGCATGTCG	302
GAACATGCCCGGGCATGTCC	300
GAACATGCCCGGGCATGTCT	300
AGACATGCCCGGGCATGTCC	298
AGACATGCCCGGGCATGTCT	298
GGACATGTCCGGGCATGTCC	298
GGACATGTCCGGGCATGTCT	298
GGACATGCCCGGGCATGTTC	298
GGACATGCCCGGGCATGTTT	298
GAACATGCCCGGGCATGTCG	295
GGGCATGCCCGGGCATGTCC	294
GGGCATGCCCGGGCATGTCT	294
GGACATGCCCGGACATGTCC	294
GGACATGCCCGGACATGTCT	294
GGACATGCCCGGGCATGCCC	294

etc.

The resulting files from both matrices are included in the zip file.