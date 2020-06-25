# Project Title

This is my personal project to try to culminate my understanding of C++. It is Matrix annd Vector class that allows for the creation of a number container in which Matrix and Vector operations can be performed on. This project was meant to be portable in graphical libraries to allow for certain algorithms, but, as of writing this, has not yet been implemented.

June update: I simplified a couple of implementations of functions and what files to clone. I plan to refactor code to make it more efficient and add more functionality in future, such as reading from a file.

## Getting Started

Simply clone or download this entire repositry to have full access to Matrices, Vectors and their associated functions and just #include "matrix.h" in your project.

### Declaring Matrices and Vectors
To declare a vector simply code: Vector a{1,2,3,4};<br/>
To declare a matrix you can code it in different ways: <br/><br/>

(i) List initialization:<br/>
Matrix m(3,3); //declare dimension first<br/>
m = {1,2,3, 3,4,5};<br/><br/>

(ii) Vector initialization:<br/>
Matrix c{a,b}; //a and b are vectors<br/><br/>

(iii) Augmented matrix initialization:<br/>
Matrix i{c,a}; //c is a matrix and a is vector<br/><br/>

Copy and assignment initialization are fully implemented
## Authors

* **Matthew Ong** - [MTau314](https://github.com/MTau314)

## License

Not licensed. Okay for personal use, acknowledgements if used publicly.

## Acknowledgments

* javidx9 for showing the possibilites of C++ *
* wikipedia and rosettacode for explanations and [pseudocode](https://rosettacode.org/wiki/Reduced_row_echelon_form#C.2B.2B) on row reducing algortihms *
