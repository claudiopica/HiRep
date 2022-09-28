# Clover Term

The clover term can be written as

$$ D_{sw} = -\frac{c_{sw}}{4}\sum_x\sum_{\mu,\nu}\sigma_{\mu\nu}F_{\mu\nu}(x), $$(eq:clover)

with the (unconventional) definition of $\sigma_{\mu\nu}$ given by

$$ \sigma_{\mu\nu} = \frac{1}{2}[\gamma_\mu,\gamma_\nu]. $$

With the Euclidean definition of the gamma matrices $\sigma_{\mu\nu}$ satisfies

$$ \sigma_{\mu\nu}^\dagger = \sigma_{\nu\mu} = -\sigma_{\mu\nu} = \sigma_{\mu\nu}^{-1}. $$

For the Hermitian Dirac operator $\gamma^5D$ we can make the following replacement without affecting any of the calculations presented here.

$$ \sigma_{\mu\nu} \to \bar{\sigma}_{\mu\nu} = \gamma_5\sigma_{\mu\nu}. $$

The field strength tensor is defined as

$$ F_{\mu\nu}(x) = \frac{1}{8}\left\{Q_{\mu\nu}(x) - Q_{\mu\nu}^\dagger(x)\right\} $$

with

$$ \begin{aligned}
\begin{split}
 Q_{\mu\nu}(x)
 &= U_\mu(x)U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x) \\
 &+ U_\nu(x)U_\mu^\dagger(x-\hat{\mu}+\hat{\nu})U_\nu^\dagger(x-\hat{\mu})U_\mu(x-\hat{\mu}) \\
 &+ U_\mu^\dagger(x-\hat{\mu})U_\nu^\dagger(x-\hat{\mu}-\hat{\nu})U_\mu(x-\hat{\mu}-\hat{\nu})U_\nu(x-\hat{\nu}) \\
 &+ U_\nu^\dagger(x-\hat{\nu})U_\mu(x-\hat{\nu})U_\nu(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x)
\end{split}
\end{aligned} $$

Because $Q_{\mu\nu}^\dagger=Q_{\nu\mu}$ we have
$F_{\mu\nu}=-F_{\nu\mu}$. For this reason we can change the sum over
$\mu,\nu$ in equation [\[eq:clover\]](#eq:clover){reference-type="eqref"
reference="eq:clover"} to a sum over $\mu<\nu$ and a factor of two.

$$ D_{sw} = -\frac{c_{sw}}{2}\sum_x\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$(eq:clover2)

The quantity $\sigma_{\mu\nu}F_{\mu\nu}$ is Hermitian and block diagonal.
It can be written as

$$ \begin{aligned}
 \sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu} =
 \begin{pmatrix}
 A & B & 0 & 0 \\
 B^\dagger & -A & 0 & 0 \\
 0  & 0 & C & D \\
 0 & 0 & D^\dagger & -C
 \end{pmatrix}
\end{aligned} $$

with the definitions

$$ \begin{aligned}
 A &= -iF_{03}+iF_{12} \\
 B &= -iF_{01}-F_{02}-F_{13}+iF_{23} \\
 C &= iF_{03}+iF_{12} \\
 D &= iF_{01}+F_{02}-F_{13}+iF_{23}
\end{aligned} $$

## Pseudofermion Forces

For the forces we use the following short-hand notation for the
derivative wrt. the link variables.

$$ \dot{S} = \partial_{x,\mu}^a S $$

### Preliminaries

To calculate the pseudofermion forces let us write down the action as

$$ S = \phi^\dagger(H^{-2})\phi, $$

where $H=\gamma^5D$ is the Hermitian Dirac operator.
When differentiating the action we obtain

$$ \dot{S} = -2\mathrm{Re}~\xi^\dagger\dot{H}\eta, $$(eq:dotS)

with the definitions
 $$\begin{aligned}
 \eta &= H^{-2}\phi, \\
 \xi &= H\eta.
\end{aligned} $$

### Forces

Here we will only consider the forces from the clover term and not the hopping term.
The clover part of the Dirac operator is given by
$$ H_{sw} = -\frac{c_{sw}}{2}\sum_{\mu<\nu}\bar{\sigma}_{\mu\nu}F_{\mu\nu}(x) $$(eq:Hsw)

When inserting equation [\[eq:Hsw\]](#eq:Hsw){reference-type="eqref" reference="eq:Hsw"} in equation [\[eq:dotS\]](#eq:dotS){reference-type="eqref" reference="eq:dotS"} we obtain

$$ \dot{S} = c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{F}_{\mu\nu}\eta). $$

From the definition of $F_{\mu\nu}$ it follows that

$$ \begin{aligned}
 \dot{S} &= \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \xi^\dagger\bar{\sigma}_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger\eta), \\
 &=
 \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \eta^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\xi).
\end{aligned} $$

This can in be written as

$$ \dot{S} = \frac{1}{8}c_{sw}\sum_{\mu<\nu} \mathrm{Re}~\mathrm{tr}\left[\dot{Q}_{\mu\nu}\left\{\bar{\sigma}_{\mu\nu}\eta(x)\otimes\xi^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi(x)\otimes\eta^\dagger(x)\right\}\right] $$

In a short hand notation we need to calculate

$$ \dot{S} = \frac{1}{8}c_{sw}\mathrm{Re}~\mathrm{tr}[\dot{Q}_{\mu\nu}(x)X_{\mu\nu}(x)] $$(eq:force)

with

$$ X_{\mu\nu}(x) = \bar{\sigma}_{\mu\nu}\eta(x)\otimes\xi^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi(x)\otimes\eta^\dagger(x) $$

This matrix has the properties $X_{\mu\nu}^\dagger=X_{\nu\mu}=-X_{\mu\nu}$.
The expression for $\dot{Q}_{\mu\nu}(x)$ contains eight different terms (two from each of the four leafs).
The eight contributions to the force can be written as

$$ \begin{aligned}
 F_1(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)X_{\mu\nu}(x)] \\
 F_2(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})X_{\mu\nu}^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})] \\
 F_3(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})X_{\mu\nu}^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})] \\
 F_4(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)X_{\mu\nu}(x+\hat{\mu})U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)] \\
 F_5(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)X_{\mu\nu}^\dagger(x+\hat{\mu})U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})] \\
 F_6(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu(x+\hat{\mu})X_{\mu\nu}(x+\hat{\mu}+\hat{\nu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)] \\
 F_7(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})X_{\mu\nu}(x+\hat{\nu})U_\nu^\dagger(x)] \\
 F_8(x) &=
 \mathrm{Re}~\mathrm{tr}[\dot{U}_\mu(x)U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})X_{\mu\nu}^\dagger(x)]
\end{aligned} $$

where each term should be multiplied by $c_{sw}/8$.
The calculation can be done efficiently by noticing that several products and terms appear in multiple places.
Introduce the intermediate variables

$$ \begin{aligned}
 Z_0 &= X_{\mu\nu}(x) \\
 Z_1 &= X_{\mu\nu}(x+\hat{\mu}) \\
 Z_2 &= X_{\mu\nu}(x-\hat{\nu}) \\
 Z_3 &= X_{\mu\nu}(x+\hat{\mu}-\hat{\nu}) \\
 Z_4 &= X_{\mu\nu}(x+\hat{\mu}+\hat{\nu}) \\
 Z_5 &= X_{\mu\nu}(x+\hat{\nu}) \\
 W_0 &= U_\mu^\dagger(x-\hat{\nu}) \\
 W_1 &= U_\nu(x-\hat{\nu}) \\
 W_2 &= U_\nu(x+\hat{\mu}) \\
 W_3 &= U_\mu^\dagger(x+\hat{\nu}) \\
 W_4 &= U_\nu^\dagger(x) \\
 W_5 &= U_\nu^\dagger(x+\hat{\mu}-\hat{\nu}) \\
 W_6 &= W_0W_1 \\
 W_7 &= W_2W_3 \\
 W_8 &= W_7W_4-W_5W_6
\end{aligned} $$

The total force can now be written as

$$ F(x) = \frac{c_{sw}}{8}\dot{U}_\mu(x)\left\{W_8Z_0 + Z_1W_8 - W_5(W_0Z_2W_1+Z_3W_6) + (W_2Z_4W_3+W_7Z_5)W_4\right\} $$

This brings us down to a total of 15 matrix multiplications and 6 additions.

## Logarithmic Forces

In the case of even-odd preconditioning (see the next section) the action of the small determinant $D_{oo}$ can be written as

$$ S_{sw} = -N_f\log\det D_{oo} = -N_f~\mathrm{tr}\log D_{oo} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\log D_{oo}(x) $$

The derivative is given by

$$ \dot{S} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\left[D_{oo}^{-1}(x)\dot{D}_{oo}(x)\right] $$

with $D_{oo}(x)$ given by

$$ D_{oo}(x) = 4+m_0-\frac{c_{sw}}{2}\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$

Both the determinant and the inverse of $D_{oo}(x)$ can be calculated from an LDL decomposition.
If we insert the above definition we obtain

$$ \dot{S} = \frac{N_fc_{sw}}{2}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{F}_{\mu\nu}(x)) $$

$$ \dot{S} = \frac{N_fc_{sw}}{2\cdot8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}
 \mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{Q}_{\mu\nu}(x) + D_{oo}^{-1}(x)\sigma_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger(x)) $$

Since $D_{oo}^{-1}$ is Hermitian we can write the result as two times the real part.
To simplify the result we define $X_{\mu\nu}(x)=D_{oo}^{-1}(x)\sigma_{\mu\nu}$ such that

$$ \dot{S} = \frac{N_fc_{sw}}{8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{Re}~\mathrm{tr}[X_{\mu\nu}(x)\dot{Q}_{\mu\nu}(x)] $$

This is equivalent to Eq. [\[eq:force\]](#eq:force){reference-type="eqref" reference="eq:force"} except from the factor $N_f$ and the definition of $X_{\mu\nu}(x)$.
Notice that we still have the identity $X_{\mu\nu}^\dagger=-X_{\mu\nu}$.
The sum over $x$ can be extended to all sites by setting $X_{\mu\nu}$ to zero on the even sites.
To calculate the inverse $D_{oo}^{-1}$ we introduce the definitions:

$$ D_{oo} = D_{oo}^\dagger = \begin{pmatrix}
 D_+ & 0 \\
 0 & D_-
 \end{pmatrix} $$

$$ D_{oo}^{-1} = \begin{pmatrix}
 D_{+}^{-1} & 0 \\
 0 & D_{-}^{-1}
 \end{pmatrix} $$

$$ D_{+}^{-1} =
 \begin{pmatrix}
  D_{11} &  D_{12} \\
  D_{21} &  D_{22} \\
 \end{pmatrix} $$

$$ D_{-}^{-1} =
 \begin{pmatrix}
  D_{33} &  D_{34} \\
  D_{43} &  D_{44} \\
 \end{pmatrix} $$

Because of hermiticity we know that $D_{12} = D_{21}^\dagger$ and $D_{34} = D_{43}^\dagger$.
The six independent elements of $X_{\mu\nu}$ can now be written as

$$ \begin{aligned}
 X_{01} &= i(D_{34}+D_{43}) - i(D_{12}+D_{21}) \\
 X_{02} &= ~(D_{12}-D_{21}) + ~(D_{43}-D_{34}) \\
 X_{03} &= i(D_{22}-D_{11}) + i(D_{33}-D_{44}) \\
 X_{12} &= i(D_{11}-D_{22}) + i(D_{33}-D_{44}) \\
 X_{13} &= ~(D_{12}-D_{21}) + ~(D_{34}-D_{43}) \\
 X_{23} &= i(D_{12}+D_{21}) + i(D_{34}+D_{43})
\end{aligned} $$

## Even-odd Preconditioning

### Method 1

We can write the determinant as

$$ \det D = \det(D_{oo})\det(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe}) $$

Use the notation

$$ \begin{aligned}
 Q_{oo} &= \gamma_5D_{oo} \\
 Q &= \gamma_5(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe})
\end{aligned} $$

The action is

$$ S = S_1 + S_2 = \phi_1^\dagger Q_{oo}^{-2}\phi_1 + \phi_2^\dagger Q^{-2}\phi_2 $$

### Forces for $\phi_1$-term

The derivative is

$$ \dot{S}_1 = -2\mathrm{Re}\left[\phi_1^\dagger(Q_{oo}^{-2}Q_{oo}\dot{Q}_{oo}Q_{oo}^{-2})\phi_1\right] $$

and we can write it as

$$ \dot{S}_1 = -2\mathrm{Re}\left[\xi^\dagger\dot{Q}_{oo}\eta\right] $$

with

$$ \begin{aligned}
 \eta &= Q_{oo}^{-2}\phi_1, \\
 \xi &= Q_{oo}\eta.
\end{aligned} $$

### Forces for $\phi_2$-term

The derivative is

$$ \dot{S}_2 = -2\mathrm{Re}\left[\phi_2^\dagger(Q^{-2}Q\dot{Q}Q^{-2})\phi_2\right] $$

and we can write it as

$$ \dot{S}_2 = -2\mathrm{Re}\left[\xi^\dagger\dot{Q}\eta\right] $$

with

$$ \begin{aligned}
 \eta &= Q^{-2}\phi_2, \\
 \xi &= Q\eta.
\end{aligned} $$

The explicit expression for $\xi^\dagger\dot{Q}\eta$ is given by

$$ \xi^\dagger \dot{Q}\eta
 = \xi^\dagger\gamma_5\dot{D}_{ee}\eta
 - \xi^\dagger\gamma_5\dot{D}_{eo}D_{oo}^{-1}D_{oe}\eta
 + \xi^\dagger\gamma_5D_{eo}D_{oo}^{-1}\dot{D}_{oo}D_{oo}^{-1}D_{oe}\eta
 - \xi^\dagger\gamma_5D_{eo}D_{oo}^{-1}\dot{D}_{oe}\eta $$

and it can be written as

$$ \xi^\dagger \dot{Q}\eta
 = \xi^\dagger\gamma_5\dot{D}_{ee}\eta
 - \xi^\dagger\gamma_5\dot{D}_{eo}\eta_1
 + \xi_1^\dagger\gamma_5\dot{D}_{oo}\eta_1
 - \xi_1^\dagger\gamma_5\dot{D}_{oe}\eta $$

with

$$ \begin{aligned}
 \eta_1 &= D_{oo}^{-1}D_{oe}\eta \\
  \xi_1 &= D_{oo}^{-1}D_{oe}\xi
\end{aligned} $$

### Method 2

The action of $D_{oo}$ can also be expressed directly as the logarithm of the determinant.

$$ S = -2\log\det D_{oo} + \phi^\dagger Q^{-2}\phi $$

This is the approach implemented in the code.

## LDL factorization

With even-odd preconditioning we need to calculate the inverse $D_{oo}^{-1}$ when applying the dirac operator and when calculating the forces.
Because this matrix is Hermitian and block diagonal it can be inverted locally with an exact solver.
The most practical solver is via an LDL decomposition.

$$ A = LDL^\dagger $$

Sum over $j$

$$ D_j = A_{jj} - \sum_{k=1}^{j-1}L_{jk}L_{jk}^*D_k $$

Sum over $i>j$.

$$ L_{ij} = \frac{1}{D_j}\left(A_{ij}-\sum_{k=1}^{j-1}L_{ik}L_{jk}^*D_k\right) $$

The determinant is given by

$$ \det(A) = \prod_k D_k $$

### LDL Decomposition

Calculates the LDL decomposition $A=LDL^\dagger$ in-place.
After the decomposition, the lower triangular part of $A$ is $L$ and the diagonal is $D$.

```fortran
do i=0, N-1
    do k=0, i-1
        A_ii = A_ii - L_ik * conj(L_ik) * A_kk
    enddo
    do j=i+1, N-1
        do k=0, i-1
	    A_ji = A_ji - A_jk * conj(L_ik) * A_kk
	enddo
	A_ji = A_ji/A_ii
    enddo
enddo
```

### Forward substitution

Calculates $x=L^{-1}b$.

```fortran
do i=0, N-1
    x_i = b_i
    do k=0, i-1
        x_i = x_i - A_ik * x_k
    enddo
enddo
```

### Backward substitution with diagonal

Calculates $x=(L^\dagger)^{-1}D^{-1}x$.

```fortran
do i=N-1, 0
    x_i = x_i/A_ii
    do k=i+1, N-1
        x_i = x_i - conj(A_ki*) * x_k
    enddo
enddo
```

### Full inversion

This algorithm calculates the inverse $B=A^{-1}$ from the LDL decomposition.
Because the inverse is Hermitian we only calculate the lower triangular part.

```fortran
do i=0, N-1
    B_ii = 1
    do j=i, N-1
        do k=i, j-1
	    B_ji = L_jk * B_ki
	enddo
    enddo
    do j=N-1, i
        B_ji = B_ji/L_ii
	do k=j+1, N-1
	    B_ji = conj(L_kj*) * B_ki
	enddo
    enddo
enddo
```
