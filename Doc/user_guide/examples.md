# Examples

## Estimation of Disconnected Contributions

### Conventions

We choose the hermitian basis of gamma matrices given in Tab. 1. Each element of the basis is referred by an index in \[0,15\] shown in the following table

| No | Matrix			  					|
|:---|:---------------------------------------|
| 0  | $\gamma_5$ 		  					|
| 1  | $\gamma_1$		  					|
| 2  | $\gamma_2$		  					|  
| 3  | $\gamma_3$		  					|
| 4  | $-\mathrm{i}\gamma_0\gamma_5$ 			|
| 5  | $-\mathrm{i}\gamma_0\gamma_1$			|
| 6  | $-\mathrm{i}\gamma_0\gamma_2$			|
| 7  | $-\mathrm{i}\gamma_0\gamma_3$	    		|
| 8  | 1										|
| 9  | $-\mathrm{i}\gamma_5\gamma_1$			|
| 10 | $-\mathrm{i}\gamma_5\gamma_2$			|
| 11 | $-\mathrm{i}\gamma_5\gamma_3$			|
| 12 | $\gamma_0$							|
| 13 | $-\mathrm{i}\gamma_5\gamma_0\gamma_1$ 	|
| 14 | $-\mathrm{i}\gamma_5\gamma_0\gamma_2$  |
| 15 | $-\mathrm{i}\gamma_5\gamma_0\gamma_3$  |

### Singlet Two-Point Functions

Consider a gauge theory on a group G coupled to $N_f$ fermions in an arbitrary representation $R$. Let us denote:

$$C(t, x_0) = \dfrac{1}{N_f}\sum_{\vec{x}}\langle \bar{q}\Gamma q(x)\bar{q}\Gamma q(x_0)\rangle$$

where $q$,$\bar{q}$ are the $N_f$ quark fields and $\Gamma$ denotes as arbitrary Dirac structure. The $1/N_f$ factor is only there for convenience. The Wick contractions read:

$$C(t, x_0) = \sum_{\vec{x}} \langle -\mathrm{tr}\left(\Gamma S(x,x_0)\Gamma S(x_0,x)\right) + N_f\,\mathrm{tr}\left(\Gamma S(x,x)\right)\mathrm{tr}\left(\Gamma S(x_0,x_0)\right)\rangle$$

### Stochastic Evaluation of Disconnected Loops

The simple one consist to evaluate stochastically the disconnected contribution without any variance reduction techniques. Considering a general volume source $\xi$, we define $\phi$ using the Dirac operator $D$:

$$\phi = D^{-1}\xi$$

For a given element X of the basis defined in the previous section, we then have

$$\sum \left(\xi^{*}X\phi\right)_{R} = \sum XM^{-1} + \text{noise}$$

where the symbol $(\ldots)_{R}$ refers to the average over R samples of the stochastic source. 

It should be observed that in evaluating the disconnected contributions to the neutral meson correlators each one of the two quark loops arising from Wick contractions must be averaged over completely _independent_ smaples of stochastic sources for the purpose of avoiding unwanted biases.

#### Implemented Source Types

We use XX noise sources. The user can switch between the following different source types

* type 0: Pure volume source
* type 1: Gauge fixed wall source
* type 2: Volume sources with time and spin dilution
* type 3: Volume sources with time, spin and color dilution
* type 4: Volume source with time, spin, color and even-odd dilution
* type 6: Volume source with spin, color and even-odd dilution

### Output

The code does not perform any average on the stochastic noise or on the dilution indices. This allows to keep as much information as possible and to vary the number of stochastic sources at the analysis level.

The indices are always

`#T #iGamma #iSrc #\[color and/or e/o \] #Re #Im`

where iGamma refers to the index of the Gamma matrix defined in Table 1.

### Debugging Options

If the code is executed with the following additional arguments

```{bash}
-p  <propagator_name> -s <source_name>
```

the code will read the two files and perform the contraction accordingly computing $\chi^{\dagger}\Gamma \psi$

