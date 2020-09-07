# Project Title

This is my personal project to try to culminate my understanding of C++. It is Matrix and Vector class that allows for the creation of a number container in which Matrix and Vector operations can be performed on. This project was meant to be portable in graphical libraries to allow for certain algorithms, but, as of writing this, has not yet been implemented.

September update: I have added a simple Matrix initialization from a file, as well as refactored parts of the code (with Move semantics implementation to reduce memory usage). 

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

(iv) File initialization:<br/>
Matrix f{3,3}; //must declare size<br/>
std::ifstream fileIn("Filesrc.txt"); //check for file open success<br/>
readMatrix(fileIn,f);<br/>

In Filesrc.txt:<br/>
1 0 2<br/>
3 12 4<br/>
9 -3 4<br/>
<br/><br/>

Copy and assignment initialization are fully implemented

### Matrices and Vectors Functions
**Both Matrices and Vectors can be outputted with <<**<br/>
**September update note: applicable operators now handle commutativity**<br/>

Vector functional functions:
<ul>
  <li>.setAll(int value) will set the vector values to value paramter</li>
  <li>.size() will return the size of vector</li>
  <li>.magnitude() will return vector magnitude</li>
  <li>dot(Vector a, Vector b) will return the dot product between a and b</li>
  <li>cross(Vector a, Vector b) will return the vector of a cross b</li>
  <li>proj(Vector a, Vector b) will return the vector from a projected onto b</li>
  <li>Vector class has a .at(int index) and operator[] for element retrieval</li>
</ul>

Matrix functional functions:
<ul>
  <li>.setAll(int value) will set the matrix values to value paramter<li>.eye() is a special setAll for an identity matrix</li></li>
  <li>.cofRank(int m, int n) will return the mxn matrix rank and .augRank() will return augmented matrix rank</li>
  <li>.det() will return determinant</li>
  <li>**.eApply(function) will apply the specified function to each element**</li>
  <li>.rowReduce() will row reduce matrix</li>
  <li>.transpose() will return a transposed matrix</li>
  <li>.dim() will matrix mxn dimension and .size() returns a vector with m and n values</li>
  <li>adj() and inv() will return matrix adjective and inverse</li>
  <li>multiply(Matrix matrix, Vector vector) and multiply(Matrix matrix1, Matrix matrix2) performs matrix multiplication</li>
  <li>Matrix class has a .at(int index), operator[], and operator() for element retrieval</li>
</ul>

## Authors

* **Matthew Ong** - [MTau314](https://github.com/MTau314)

## License

Not licensed. Okay for personal use, acknowledgements if used publicly.

## Acknowledgments

* javidx9 for showing the possibilites of C++ *
* wikipedia and rosettacode for explanations and [pseudocode](https://rosettacode.org/wiki/Reduced_row_echelon_form#C.2B.2B) on row reducing algortihms *
