# PWM Sequence Generator

**“Position Weight Matrices are a simple way to model signals appearing on DNA and protein sequences. They summaryze the frequencies of every letter of the nucleotide or amino acid alphabets at a given position of the signal, for instance a splice site or a transcription factor binding site. This means that we can obtain a score, either from the absolute frequencies or by transforming them into log-likelihoods, for a given sequence to determine if it contains the signal pattern or not. Yet, they can also be used as sequence generators.”**

This PERL script will take as input a Position Weight Matrix and will generate all of the possible sequences from this PWM.

Some guidelines we need to take into account:
- We will only work with nucleotides
- We will set up a cut-off at 1000 sequences: Depending on the way our PWM is set up, and being the complexity of this program exponential, we could end up with an execution time too long to handle. Thus the reason we need a cut-off value.
- Our program will use a recursive function

## Execute the script:

`perl pwm_sequence_generator.pl pwm_filepath [output_filepath]`

Arguments:
- `pwm_filepath` → the filepath where our PWM is stored, in TRANSFAC format.
- `[output_filepath]` (optional)
    - If given → The resulting sequences will be written into the file at output_filepath. If the file doesn’t exist, it’ll be created
    - If **NOT** given → The resulting sequences will be printed to the STDOUT.
