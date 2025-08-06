
This is *single header file* that implements a set of double-double primitive operations.

** IMPORTANT**:  any relaxation of floating-point rules (via command line options) will break a
reasonable number of functions. This is includes enabling so-call floating point contractions
(automatic usage of FMAs) which is typically enabled by default. Exceptions to this are 
*math-errno* which should have stop being a thing decades ago and *trapping math*.

Spitball overview:

* modern versions of classic methods: FMA variants and doesn't replace divides by a large
number of products and additions
* Lange and Rump's relaxed (CPair) methods
* 

