# Project Title

This is my personal project to try to culminate my understanding of C++. It is Matrix annd Vector class that allows for the creation of a number container in which Matrix and Vector operations can be performed on. This project was meant to be portable in graphical libraries to allow for certain algorithms, but, as of writing this, has not yet been implemented.

## Getting Started

Simply clone or download this entire repositry to have full access to Matrices, Vectors and their associated functions

### Sample
Start by #include "mvmath.h"

- Initialize a Vector as
    Vector v{1,2,3} or 
    Vector w(1,2,3) or 
    Vector u = v

- Run vector functions included in mvmath.h under VMath namespace (To see implementation, simply take a look at mvmath.cpp) --
    VMath::proj(v,w)
- Run vector functions include in vector.h --
    v.setAll(1) or v.magnitude()

- Initialize a Matrix as
    Matrix A(2,3) and 
    A = { 1, 2, 3,
          4, 5, 6 }
    (This initialization starts with your explicit declaration of Matrix size, followed by the input of values); or 
          
    Matrix B(Vector{1,2,3}, Vector{1,2,3}, Vector{1,2,3}) or
    Matrix C = A;
    
- Run matrix functions included in mvmath.h under MMath namespace --
  MMath::inv(B)
- Run matrix functions include in matrix.h --
  B.rowReduce() or B.setAll(1) or B.augRank()

## Authors

* **Matthew Ong** - [MTau314](https://github.com/MTau314)

## License

Not licensed. Okay for personal use, acknowledgements if used publicly.

## Acknowledgments

* javidx9 for showing the possibilites of C++ *
* wikipedia and rosettacode for explanations and [pseudocode](https://rosettacode.org/wiki/Reduced_row_echelon_form#C.2B.2B) on row reducing algortihms *
