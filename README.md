# PAS2RS

PAS2RS is a simple tool for protein absolute structure transformation to relative structure.

## Running the script

The tool can be run two different ways. Either with input file passed as parameter or input passing through pipe.

```bash
# reading from STDIN
> cat input_file | ./PAS2RS.sh
```


```bash
# reading from file
> ./PAS2RS.sh input_file 
```

## Input

Protein absolute structure is a sequence of 3D coordinates of central carbon of every residue in protein.

PAS2RS expect relative structure to be passed through STDIN or regular text file.
Either way, every line of input has to consist of 3 floats representing 3 coordinates of a given central carbon.
See the [example input file](https://github.com/michalfilippi/PAS2RS/blob/master/test/coordinates.txt).

## Output

Relative structure of a protein consist of 5 features per every residue in protein.
That is central carbon angles (CCA), forward central carbon plane torsion angles (CCPTA_F) and directions (CCPTD_F) and forward central carbon plane torsion angles (CCPTA_B) and directions (CCPTD_B). 

Output is printed to the STDOUT. 

### Output example

```bash
> cat test/coordinates.txt | ./PAS2RS.sh

CCA CCPTA_F CCPTD_F CCPTA_B CCPTD_B
0 0 0 0 0
0.333333333333 0 0 0.608172322575 -1.0
0.333333333333 0.608172322575 -1.0 0.695913171072 1.0
0.5 0.695913171072 1.0 0.695913171072 -1.0
0.333333333333 0.695913171072 -1.0 0.0 0.0
0.333333333333 0.0 0.0 0.391827677425 1.0
6.70787927625e-09 0 0 0 0
0.333333333333 0.391827677425 1.0 0.391827677425 1.0
0.333333333333 0.391827677425 1.0 0.608172322575 -1.0
0.333333333333 0.608172322575 -1.0 0.391827677425 1.0
0.333333333333 0.391827677425 1.0 0.0 0.0
0.666666666667 0.0 0.0 1.0 0.0
6.70787927625e-09 0 0 0 0
0.666666666667 1.0 0.0 0.0 0.0
0.333333333333 0.0 0.0 0.391827677425 -1.0
0.666666666667 0.391827677425 -1.0 1.0 0.0
0.666666666667 1.0 0.0 0.0 0.0
0.333333333333 0.0 0.0 0.608172322575 1.0
6.70787927625e-09 0 0 0 0
0.666666666667 0.608172322575 1.0 0 0
0 0 0 0 0

```

