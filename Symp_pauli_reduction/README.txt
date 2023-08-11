README for Symplectic Pauli Reductions

Code related to "Occam's razor for quantum simulations"
By Lane Gunderman, Andrew Jena, and Luca Dellantonio

Project description:
This method takes in a collection of Pauli operators transforms them via a Clifford operation into minimal representations with separable registers (isotropic subspace) at the end (those only using I/Z), while the other registers provide the non-separable registers (symplectic subspace). Coefficients are discarded in this method and must be tracked separately, however, the i-th original Pauli operator is transformed into the i-th minimal representation Pauli operator.

Dependencies: numpy, (and random, math, os)

Use: input is a textfile, or folder of textfiles, with each row being a single Pauli operator expressed with I,X,Y,Z terms. Beyond this formatting is not important as these characters are picked out. Capitalization is required.
Output (and options): a priori the code states the number of Pauli operators in the file, the length (number of qubits) of the original Pauli operators, and the minimal number of qubits needed to simulate those operations (half the dimension of the symplectic subspace).
The code can also output the equivalent Pauli operators that obtain this minimum (along with the separable registers), "pauli_forms".

Credit: Lane Gunderman

Included for testing are the following Pauli files:
hydrogen (H2) in STO3G: n=4, c=3, half-symp rank=1
hydrogen (H2) in 6-31G: n=8, c=3, half-symp rank=5
methane (CH4) in STO3G: n=18, c=3, half-symp rank=15
ethyne (C2H2) in STO3G: n=24, c=4, half-symp rank=20

Improvements for the future:
1) Turn numpy arrays into numpy Boolean arrays
2) Add prime qudit functionality, which will require integer arrays
3) Add non-prime qudit functionality, which is ongoing research involving Smith normal forms