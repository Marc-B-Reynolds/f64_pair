
This is *single header file* that implements a set of double-double primitive operations.

**IMPORTANT**:  any relaxation of floating-point rules (via command line options) will break a
reasonable number of functions. This is includes enabling so-call floating point contractions
(automatic usage of FMAs) which is typically enabled by default. Exceptions to this are 
*math-errno* which should have stopped being a thing decades ago and *trapping math*.

Spitball overview:
* REFERENCES! (probably the most useful bit)
* More of copy/paste/modify set of code snippets than library
* modern versions of classic methods: FMA variants and doesn't replace divides by a large
number of products and additions
* Lange and Rump's relaxed (CPair) methods
* *Graillat & Muller* 2025 methods reworked to fast-path/slow-path and bit manipulation for the test. The *fast-path* is expected to execute ~100% of inputs.
* some special casing of standard methods
* some (not to be trusted too much) homegrown routines

