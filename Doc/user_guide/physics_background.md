# Supported Features

## Introduction

### Conventions

This section summarizes the main formulae that are used for implementing
the HMC for dynamical Wilson fermions in higher representations. The
Dirac operator is constructed following Ref. [@Luscher:1996sc], but
using Hermitian generators 

$$T^{a\dagger}=T^a.$$ 

For the fundamental
representation, the normalization of the generators is such that:

$$\mathrm{tr}\, \left(T^a T^b \right) = \frac12 \delta^{ab}.$$

For a generic representation $R$, we define: 

$$\mathrm{tr }_R \left(T^a T^b \right) = T_R \delta^{ab},$$

$$\sum_a \left(T^a T^a \right)_{AB} = C_2(R) \delta_{AB},$$

which implies: 

$$T_R = \frac{1}{N^2-1} C_2(R) d_R$$ 

where $d_R$ is the dimension of the representation $R$. The relevant group factors may be computed from the Young tableaux of the representation of $SU(N)$ by using the formula:

$$C_2(R) =\frac{1}{2}\left(nN+ \sum_{i=1}^{m} n_i \left( n_i+1-2i
\right) - \frac{n^2}{N}\right)$$ 

where $n$ is the number of boxes in the
diagram, $i$ ranges over the rows of the Young tableau, $m$ is the
number of rows, and $n_i$ is the number of boxes in the $i$-th row.

|    R |        $d_R$        |       $T_R$     |           $C_2(R)$          |
|------|---------------------|-----------------|-----------------------------|
| fund |         $N$         |     $\frac12$   |     $\frac{N^2-1}{2 N}$     |
| Adj  |      $N^2-1$        |       $N$       |            $N$              |
|  2S  | $\frac{1}{2}N(N+1)$ | $\frac{N+2}{2}$ |  $C_2(f)\frac{2(N+2)}{N+1}$ |
|  2AS | $\frac{1}{2}N(N-1)$ |  $\frac{N-2}{2}$|  $C_2(f)\frac{2(N-2)}{N-1}$ |


A generic element of the algebra is written as: $X=i X^a T^a$, and the
scalar product of two elements of the algebra is defined as:

$$(X,Y)= \mathrm{tr\ } \left(X^\dagger Y\right) = T_f X^a Y^a,$$

$$\Vert X \Vert^2 = \mathrm{tr } \left(X^\dagger X\right)
 = \sum_{ij} \left| X_{ij} \right|^2$$
 
#### $\gamma$ matrices

We use the chiral representation for the Dirac $\gamma$ matrices:

$$\gamma_\mu=\begin{pmatrix}0&e_\mu\\e_\mu^\dagger&0\\\end{pmatrix}\, ,$$ 

where $e_\mu$ are $2\times 2$ matrices given by $e_0=-1$, $e_k=-i\sigma_k$, 

$$\sigma_1=
\begin{pmatrix}
0&1\\
1&0\\
\end{pmatrix},\,\,
\sigma_2=
\begin{pmatrix}
0&-i\\
i&0\\
\end{pmatrix},\,\,
\sigma_3=
\begin{pmatrix}
1&0\\
0&-1\\
\end{pmatrix}\, .$$ 

We have:

$$\gamma_5=\gamma_0\gamma_1\gamma_2\gamma_3=
\begin{pmatrix}
1&0\\
0&-1\\
\end{pmatrix}\, .$$

### The Dirac operator

The massless Dirac operator is written as in Ref. [@Luscher:1996sc]:

$$D = \frac12 \left\{\gamma_\mu \left(\nabla_\mu + \nabla^*_\mu \right) 
- \nabla^*_\mu \nabla_\mu \right\}$$ 

with 

$$\nabla_\mu\phi(x) = U^R (x,\mu)\phi(x+\mu) - \phi(x)$$
$$\nabla_\mu^*\phi(x) = \phi(x) - U^R (x-\mu,\mu)^\dagger\phi(x-\mu)$$

and therefore the action of the massive Dirac operator yields:

$$\begin{aligned}
 D_m \phi(x) =& (D+m) \phi(x)\\ 
=& - \frac12 \left\{ \left(1-\gamma_\mu\right) U^R(x,\mu) \phi(x+\mu) \right.\\
&+
\left.\left(1+\gamma_\mu\right) U^R(x-\mu,\mu)^\dagger \phi(x-\mu)-(8+2m) \phi(x) \right\}, \end{aligned}$$(eq:DM) 

where $U^R$ are the link variables in the representation $R$.

Rescaling the fermion fields by $\sqrt{\kappa}=\left(\frac{2}{8+2m}\right)^{1/2}$, we can write the fermionic action as:

$$S_f = \sum_{x,y} \phi^\dagger(x) D_m(x,y) \phi(y),$$ 

where

$$D_m(x,y) = \delta_{x,y} - \frac{\kappa}{2}
\left[(1-\gamma_\mu) U^R(x,\mu) \delta_{y,x+\mu} + 
(1+\gamma_\mu) U^R(x-\mu,\mu)^\dagger \delta_{y,x-\mu} \right],$$ 

and the Hermitian Dirac operator is obtained as:

$$Q_m = \gamma_5 D_m. $$(eq:QM) 

The fermionic determinant in the path integral can be represented by introducing complex pseudofermionic fields: 

$$\left(\det D_m\right)^{N_f} = \int \mathcal D \phi \mathcal D \phi^\dagger e^{-\phi^\dagger Q_m^{-N_f} \phi} \equiv \int \mathcal D \phi \mathcal D \phi^\dagger e^{-S_\mathrm{pf}}.$$

### Force for the HMC molecular dynamics 

The HMC Hamiltonian is given by:

$$\mathcal{H}=\mathcal{H}_\pi+\mathcal{H}_G+\mathcal{H}_F \, ,$$ 

where

$$\mathcal{H}_\pi = \frac{1}{2} \sum_{x,\mu} ( \pi(x,\mu) , \pi(x,\mu) ) = \frac{1}{2} T_f \sum_{a,x,\mu} \pi^a(x,\mu)^2 \, ,$$

$$\mathcal{H}_G = \beta \sum_{\mu<\nu} \left( 1- \frac{1}{N} \mathrm{Re\ tr\ } \mathcal{P}_{\mu\nu}\right) \, ,$$

$$\mathcal{H}_F = \phi^\dagger ( Q_m^2 - \beta )^{-l} \phi \, , \,\,\,\, l=\frac{N_f}{2}>0 \, ,$$(eq:HF)

and we have introduced for each link variable a conjugate momentum in
the algebra of the gauge group: $\pi(x,\mu)=i \pi^a(x,\mu) T_f^a$. In
the expression of $\mathcal{H}_F$ we omitted the sum over position, spin
and color indices and we have also introduced an arbitrary shift $\beta$
for the matrix $Q_m^2$, as this will be useful in the discussion for the
RHMC algorithm.

The equation of motion for the link variables are given by (the
$\dot{\square}$ indicates the derivative with respect to the molecular
dynamics time): 

$$\dot U(x\mu) = \pi(x,\mu) U(x,\mu)\, ,$$ 

while the equation of motion for the momenta can be obtain as follows from the
requirement that the hamiltonian $\mathcal{H}$ is a conserved quantity:

$$0 = \dot{\mathcal{H}} = \dot{\mathcal{H}}_\pi + \dot{\mathcal{H}}_G + \dot{\mathcal{H}_F} \, .$$(eq:HCONS)

For the first two derivatives we have: 

$$\dot{\mathcal{H}}_\pi = \sum_{x,\mu} ( \pi(x,\mu) , \dot\pi(x,\mu) ) = T_f \sum_{x,\mu} \sum_a \pi^a(x,\mu) \dot\pi^a(x,\mu) \, $$(eq:HDOTPI)

$$\begin{aligned}
\dot{\mathcal{H}}_{G} 
&= \sum_{x,\mu} -\frac{\beta}{N} \mathrm{Re\, tr} \left(\dot U(x,\mu) V^\dagger(x,\mu) \right) \\
&= \sum_{x,\mu} -\frac{\beta}{N} \mathrm{Re\, tr} \left(\pi(x,\mu) U(x,\mu) V^\dagger(x,\mu) \right) \\ 
&= \sum_{x,\mu} \sum_a -\frac{\beta}{N} \pi^a(x,\mu) \mathrm{Re\, tr} \left(i T^a_f U(x,\mu) V^\dagger(x,\mu) \right) \, , 
\end{aligned}$$(eq:HDOTG)

where $V(x,\mu)$ is the sum of the staples around the link $U(x,\mu)$.

The computation of the fermionic force goes as follows. We only consider
the case $l=1$ since this is the only case relevant both for the HMC
algorithm and the RHMC algorithm (see below). We have: 

$$\begin{aligned} \dot{\mathcal{H}}_F = -\ \phi^\dagger (Q_m^2 - \beta)^{-1} \dot{(Q_m^2)} (Q_m^2 - \beta)^{-1} \phi \, . \end{aligned}$$(eq:FF1)

Defining: 

$$\eta = (Q_m^2 - \beta)^{-1} \phi \, , $$(eq:HMCETA)
$$\xi = Q_m \eta \, ,$$ 

and using the fact that the matrix $(Q_m^2-\beta)$ is hermitian, we can rewrite {eq}`eq:FF1` as

$$\begin{aligned} \dot{\mathcal{H}}_F = - 2 \ \xi^\dagger \dot{(Q_m)} \eta \, .\end{aligned}$$(eq:FF2)

Inserting the explicit form of $Q_m$, eqs. {eq}`eq:QM` and {eq}`eq:DM` into eq. {eq}`eq:FF2` we obtain 

$$\begin{aligned}\dot{\mathcal{H}}_F &= \mathrm{Re\ }\sum_{x,\mu} \xi(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \eta(x+\mu) + \xi(x+\mu)^\dagger \dot U^R(x,\mu)^\dagger \gamma_5 (1+\gamma_\mu) \eta(x) \\&= \mathrm{Re\ }\sum_{x,\mu} \xi(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \eta(x+\mu) + \eta(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \xi(x+\mu)\end{aligned}$$

where the sum over spin and color indices is intended and we made
explicit the fact the the whole expression is real. We now use the fact that

$$\dot U^R (x,\mu) = \pi^R(x,\mu) U^R(x,\mu) = i \pi^a(x,\mu) T^a_R U^R(x,\mu) $$(eq:URDOT)

Notice that, since we define $T^a_R(x,\mu) = R_* T^a(x,\mu)$, the
$\pi^a(x,\mu)$ in the above equation are the same as those appearing in
the expressions for $\dot{\mathcal{H}}_{\pi,G}$. Using eq. {eq}`eq:URDOT` in the
expression for $\dot{\mathcal{H}}_{F}$ we find: 

$$\begin{aligned}
\dot{\mathcal{H}}_F = \sum_{x,\mu} \sum_a  \pi^a(x,\mu) & \mathrm{Re\ Tr\ } \left[ iT^a_R U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \right. \\
&\left. \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right] \, .\end{aligned}$$(eq:HDOTF)

Note that capitalized $\mathrm{Tr}$ indicates the trace over both color and spin indices as opposed to the lower case $\mathrm{tr}$, which is the trace over color only.

Inserting eq.s {eq}`eq:HDOTPI`, {eq}`eq:HDOTG` into eq. {eq}`eq:HCONS` we obtain the equations of motion for the momenta $\pi^a(x,\mu)$ 

$$\begin{aligned}
\dot\pi^a(x,\mu) &= \dot\pi^a_G(x,\mu) + \dot\pi^a_F(x,\mu) \, , \\
\dot\pi^a_G(x,\mu) &= \frac{\beta}{N} \frac{1}{T_f} \mathrm{Re\ tr\ } \left[ i T^a_f U(x,\mu) V^\dagger(x,\mu) \right] \, ,\\
\dot\pi^a_F(x,\mu) &=-\frac{1}{T_f} \mathrm{Re\ Tr\ } \left[ iT^a_R U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \right. \nonumber\\
                                        &\quad\quad\quad    \left. \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right]\, . \end{aligned}$$(eq:PIDOT3)

For sake of convenience we introduce the following projectors $P^a_R$ over the algebra in the representation $R$

$$P^a_R ( F ) = - \frac{1}{T_R} \mathrm{Re\ tr\ } \left[ i T^a_R F \right] \, ,$$

which can be used to rewrite eq.s {eq}`eq:PIDOT2` and {eq}`eq:PIDOT3` in a more compact form: 

$$\begin{aligned}
\dot\pi^a_G(x,\mu) &= - \frac{\beta}{N} P^a_f \left( U(x,\mu) V^\dagger(x,\mu) \right) \, ,\\
\dot\pi^a_F(x,\mu) &= \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \right. \right. \nonumber\\
&\quad\quad\quad    \left. \left. \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right] \right)\, . \end{aligned}$$(eq:HFFORCE)

### Checks of the MD force

The formulae derived in the previous Section can be checked against two
known examples. The first, and almost trivial, check is obtained by
assuming that the representation $R$ is again the fundamental
representation. The well-known expression for the MD force for the usual
HMC is then recovered.

The second case that has already been studied in the literature is the
case of fermions in the adjoint representationof the gauge group
SU($2$) [@Donini:1996nr]. We agree with eq. (16) in
Ref.  [@Donini:1996nr], provided that we exchange the indices $a$ and
$b$ in that formula.

### HMC Algorithm

We briefly review the construction of the HMC algorithm [@{??}].

Given the action $S(\phi)$ of a system of bosonic fields $\phi$, our
goal is to generate a Markov process with fixed probability distribution
$P_S(\phi) = 1/Z \exp[-S(\phi) ]$. A sufficient condition to have such a
Markov process is that it is ergodic and it satifies detailed balance:

$$P_S(\phi)P_M(\phi\rightarrow \phi') = P_S(\phi')P_M(\phi' \rightarrow \phi) \, .$$

We define $P_M(\phi \rightarrow \phi')$ with the following three-step
process:

1.  We expand the configuration space with additional fields, the
    "momenta\" $\pi$ randomly chosen with probability $P_k(\pi)$ such
    that $P_k(\pi)=P_k(-\pi)$ -- usually one takes
    $P_k(\pi)\propto \exp[-\pi^2/2]$;

2.  In the extended configuration space $(\phi, \pi)$, we generate a new
    configuration $(\phi',\pi')$ with probability
    $P_h((\phi,\pi)\rightarrow(\phi',\pi'))$ such that

    $$P_h((\phi,\pi)\rightarrow(\phi',\pi')) = P_h((\phi',-\pi')\rightarrow(\phi,-\pi))$$
    (reversibility condition)

3.  We accept the new configuration $\phi'$ with probability

    $$P_A((\phi,\pi)\rightarrow(\phi',\pi')) = min \left\{ 1, \frac{P_S(\phi')P_k(\pi')}{P_S(\phi)P_k(\pi)} \right\} \, .$$

    It is easy to see that the resulting probability

    $$P_M(\phi\rightarrow\phi') = \int d\pi d\pi' P_k(\pi) P_h((\phi,\pi)\rightarrow(\phi',\pi')) P_A((\phi,\pi)\rightarrow(\phi',\pi')) \, ,$$

    satisfies detailed balance. Care must be taken to ensure ergodicity.

As already stated, the distribution $P_k(\pi)$ is generally taken to be
gaussian (this should also guarantee ergodicity). The process $P_h$ is
instead identified with the Hamiltonian flow of a yet unspecified
Hamiltonian $H$ in the phase space $(\phi,\pi)$ (giving to $\pi$ the
meaning of "momenta"). The time reversal symmetry of classical dynamics
equation of motion guarantees the reversibility condition. The resulting
probability $P_h$ is then a delta function (the process is completely
deterministic). Numerical integration to a given accuracy will result in
a broader distribution and care must be taken to guarantee the
reversibility condition in this case. Since we want a high acceptance
rate (low correlation among the configurations), we must carefully
choose the Hamiltonian $H$. One simple way is to take $P_k$ to be
gaussian and define $H(\pi,\phi)=-\ln [P_k(\pi) P_S(\phi)] = \pi^2/2 + S(\phi)$ (omitting irrelevant constants). If $H$ is exactly conserved by the process $P_h$
then the acceptance probability is 1.

When fermionic degrees of freedom are present in the action $S$, we can
first integrate them out, resulting in a non local bosonic action and
then apply the above scheme. In practice, to deal with a non-local
action is not convienent from a numerical point a view and stochastic
estimates are used.

Consider a quadratic fermionic term in the action:
$S(\bar\psi,\psi) = \bar\psi M \psi$ with a generic interaction matrix
$M(\phi)$ function of the bosonic fields $\phi$. The contribution of
this term to the partition function is
$\int d\bar\psi d\psi \exp [ -S(\bar\psi,\psi)] = \mathrm{det}[M(\phi)]$.

Assuming that the matrix $M(\phi)$ is positive definite, we can rewrite
$\mathrm{det}[M]=\int d\bar\eta d\eta \exp[ \bar\eta (M)^{-1} \eta ]$, where
$\bar\eta$,$\eta$ are two new complex bosonic fields, called
pseudofermions. This term can be taken into account generating random
pseudofermions $\bar\eta$, $\eta$ with the desidered probability
distribution and keeping then fixed during the above HMC configuration
generation for the remaining bosonic fields $\phi$.

### RHMC formulation

The fermionic part of the HMC hamiltonian, for $N_f$ degenerate quarks
and $N_{pf}$ pseudofermions, can be written as:

$$\mathcal{H}_F = \sum_{k=1}^{N_{pf}} \phi_k^\dagger ( Q_m^2 )^{-l_k} \phi_k \,\, ;\,\, \sum_k l_k = \frac{N_f}{2}\, , $$(eq:HFN)

and $l_k>0$. For the sake of simplicity we will set all the $l_k$ to be
equal: 

$$\forall k,\,\, l_k = \frac{N_f}{2N_{pf}}\, .$$

In the RHMC algorithm [@Clark:2005sq] rational approximations are used
whenever we need to take some fractional power of the positive definite
fermion matrix $Q_m^2$.

In this implementation we use three different rational approximations.

The first one is used to approximate eq. {eq}`eq:HFN` (we need
only one approximation because all $l_k$ are equal): 

$$\mathcal{H}_F = \sum_{k=1}^{N_{pf}} \phi_k^\dagger r_{a}( Q_m^2 )\phi_k \, , $$(eq:HFRHMC)
 
$$( Q_m^2 )^{-\frac{N_f}{2N_{pf}}} \simeq r_{a}(Q_m^2) = \alpha_0^a + \sum_{n=1}^{d_{1}} \alpha_n^a ( Q^2_m - \beta_n^a )^{-1} \, .$$

Using the formulas derived in the previous sections, it is easy to write the force corresponding to eq. {eq}`eq:HFRHMC`. In fact, eq. {eq}`eq:HFRHMC` is nothing but a sum of terms of the form eq. {eq}`eq:HFRHMC` once we put $l=1$, $\beta=\beta_n^a$. The RHMC force will be then given by a sum over $n=1,\dots,d_1$ of terms given by
eq. {eq}`eq:HFFORCE` multiplied by a factor $\alpha_n^a$. Notice that since $l=1$, to compute $\eta$ as in eq. {eq}`eq:HMCETA` a simple shifted inversion is required.

The second rational approximation is required in the heat bath update of
pseudofermions. In order to generate pseudofermions distributed as in
eq. {eq}`eq:HFN`, a simple two-step process is used. For each pseudofermion we first generate a gaussian distributed field $\tilde\phi_k$:

$$P(\tilde\phi_k)\propto \exp [ -\tilde\phi_k^\dagger \tilde\phi_k ] \, ,$$

and then we set: 

$$\phi_k = (Q_m^2)^{\frac{l_k}{2}} \tilde\phi_k \, ,$$

making use of the fact that $(Q_m^2)$ is hermitean (notice the plus sign
in the exponent.) The RHMC algorithm uses a rational approximation to
compute the above quantities (again we need only one approximation since
all $l_k$ are equal): 

$$\begin{aligned}
( Q_m^2 )^{\frac{l_k}{2}} \simeq r_{b}(Q_m^2) &=& \alpha_0^b + \sum_{n=1}^{d_{2}} \alpha_n^b ( Q^2_m - \beta_n^b )^{-1} \, .\end{aligned}$$

The third rational approximation is used in the code for the Metropolis
test. Starting from eq. {eq}`eq:HFN` for each pseudofermion we can rewrite:

$$\phi_k^\dagger ( Q_m^2 )^{-l_k}\phi_k = \left\| (Q_m^2)^{-\frac{l_k}{2}} \phi_k \right\|^2\, ,$$

where we used the property that $Q_m^2$ is hermitean. The rational
approximation needed in this case is: 

$$\begin{aligned}
( Q_m^2 )^{-\frac{l_k}{2}} \simeq r_{c}(Q_m^2) &=& \alpha_0^c + \sum_{n=1}^{d_{3}} \alpha_n^c ( Q^2_m - \beta_n^c )^{-1} \, .\end{aligned}$$

Notice that if $d_2=d_3$ the coefficients for the two approximations
$r_b$ and $r_c$ can each be obtained from the other.

In order to compute the coefficients $\alpha_n$, $\beta_n$ appearing in
the rational approximations the Remez algorithm is needed. In this
implementation we do not compute those coefficients "on the fly", but
rather we use a precomputation step to generate a table of coefficients
form which we pick up the right values when needed. The generation of
this table goes as follows.

First note that we need to compute rational approximations for a
function $f(x)$ of the form $f(x)=x^l$ and the approximation must be
accurate over the spectral range of the operator $Q_m^2$. To simplify
the computation of the table we note that the following proposition
holds: if $f(x)$ is a homogeneous function of degree $l$ and $r(x)$ is
an optimal (in the sense of relative error) rational approximation to
$f(x)$ over the interval $[\epsilon,\mathrm{h}]$ to a given accuracy
then $r(kx)/k^l$ is an optimal rational approximation for the same
function and the same accuracy over the interval
$[\epsilon/k,\mathrm{h}/k]$. Notice that the coefficients of the
"rescaled" rational approximation are easily obtained from that of the
original approximation. A simple corollary is that, given a homogeneuos
function $f(x)$, we can divide the rational approximations with the same
accuracy in classes distinguished by the ratio $\epsilon/\mathrm{h}$;
within each class the coefficients of the rational approximations are
easily related to each other, so that we only need to compute one
rational approximation in each class. This is what is done in our
implementation.

In detail: we generate a table containing the coefficients for the
rational approximations belonging in different classes distinguished by
the function $f(x)$ which we want to approximate and the accuracy which
is required. We arbitrary set $\mathrm{h}$ to a fixed value equal to the
absolute upper bound on the spectrum of the matrix $Q_m^2$. This choice
fixes the representative of each class, because the lower bound of the
approximation is now a function of $\mathrm{h}$.

At run-time this table is used to generate optimal rational
approximations rescaling the precomputed coefficients to the desired
interval containing the spectrum of the matrix $Q_m^2$. This interval is
obtained by computing the maximum and minimum eigenvalue of $Q_m^2$ on
each configuration when needed. In our code we update this interval only
before the metropolis test, while we keep it fixed during the molecular dynamics.

### Even-Odd preconditioning

It is a very well know fact that the time spend for a simulation with
dynamical fermions is dominated by the time required for the inversions
of the Dirac operator. The convergence of such inversions can be
improved using an appropriate precondining. The idea is to rewrite the
fermionic determinant as a determinant (or product of determinants) of
better conditioned matrix (matrices) than the original Dirac operator.
For the non-improved Wilson action this can be easily done using the
*even-odd* preconditioning. We start rewriting the Dirac operator $D_m$
as a block matrix: 

$$D_m = \begin{pmatrix}
4+m& D_{eo}\\
D_{oe} &4+m\\
\end{pmatrix}\,\,\, ,$$ 

where each block has a dimension half that of the original Dirac matrix. The diagonal blocks connecting sites with the same parity are proportional to the identity matrix, while off-diagonal blocks connect sites with opposite parity. We have (since $D_m$ is $\gamma_5$-hermitean):

$$\gamma_5 D_{eo} \gamma_5 = D_{oe}^\dagger\,\, .$$ 

The determinant of the Dirac matrix $D_m$ can be rewritten as:

$${\rm det\ } D_m = {\rm det\ } ( (4+m)^2 - D_{oe} D_{eo} ) = {\rm det\ } ( (4+m)^2 - D_{eo} D_{oe} ) \equiv {\rm det\ } D^{eo}_m\,\, ,$$

using the well known formula for the determinant of a block matrix.
Since the determinant of $D_m$ and of $D_m^{eo}$ are the same the latter
can be used in numercal simulations. Note that the even-odd
preconditioned matrix only connects sites with the same parity thus it
have only half of the size of the original Dirac matrix and as $D_m$ it
is $\gamma_5$-hermitean. We define as before the hermitean matrix
$Q_m^{eo}\equiv \gamma_5 D_m^{eo}$, which will be used in practice.

The formulation of the HMC algorithm does not change and the only
difference is that pseudofermions fields are now only defined on half of
the lattice sites, conventionally the even sites in what follows. We now
give the explicit expression for the fermionic force for the
preconditioned system described by the hamiltonian: 

$$
\mathcal{H}_F = \phi_e^\dagger ( (Q^{eo}_m)^2 - \beta )^{-1} \phi_e \,\, ,$$

where as before we are assuming $N_f=2$ or a rational approximation of
the actual fractional power function, and where we made explicit that
$\phi_e$ is only defined on even sites. Eq. {eq}`eq:FF2` is unchanged: 

$$\dot{\mathcal{H}}_F = - 2 \ \xi_e^\dagger \dot{(Q^{eo}_m)} \eta_e \, ,$$(eq:FFPRE)

where as before we have defined: 

$$\eta_e = ((Q^{eo}_m)^2 - \beta)^{-1} \phi_e \, ,$$

$$\xi_e = Q^{eo}_m \eta_e \, .$$ 

The explicit form of $Q_m^{eo}$ must be used at this point. We have: 

$$\dot{(Q^{eo}_m)} = -\gamma_5 (\dot{D_{eo}} D_{oe} + D_{eo}\dot{D_{oe}} )\,\, .$$(eq:QPREDOT)

Defining 

$$\sigma_o = D_{oe} \eta_e \, , $$

$$\rho_o = D_{oe} \xi_e \, ,$$ 

and inserting eq. {eq}`eq:QPREDOT` into eq. {eq}`eq:FFPRE` we find:

$$\begin{aligned}
\dot{\mathcal{H}}_F = - \sum_{\mu,x\in even} {\rm Tr}_{x,\mu}  \left[ \sigma_o(x+\mu)\otimes\xi_e(x)^\dagger + \rho_o(x+\mu)\otimes\eta_e(x)^\dagger \right] \\
- \sum_{\mu,x\in odd} {\rm Tr}_{x,\mu}  \left[ \xi_e(x+\mu)\otimes\sigma_o(x)^\dagger + \eta_e(x+\mu)\otimes\rho_o(x)^\dagger \right] \end{aligned}$$(eq:FORPRE)

and for convenience we use the shorthand notation:

$${\rm Tr}_{x,\mu} \left[ \Phi \right] \equiv \mathrm{Re\ Tr\ } \left[ \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \Phi \right]\,\, .$$

From eq. {eq}`eq:FORPRE` it is clear that the fermionic force has a different expression on sites of different parities. Proceeding as before we arrive at the final
expressions. For $x\in even$: 

$$\begin{aligned} \dot\pi^a_F(x,\mu) &= - \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \right. \right. \\
&\quad\quad\quad    \left. \left. \left\{  \sigma_o(x+\mu)\otimes\xi_e(x)^\dagger + \rho_o(x+\mu)\otimes\eta_e(x)^\dagger \right\} \right] \right)\, ,\end{aligned}$$

while for $x\in odd$: 

$$\begin{aligned}
\dot\pi^a_F(x,\mu) &= - \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \right. \right. \nonumber\\
&\quad\quad\quad    \left. \left. \left\{\xi_e(x+\mu)\otimes\sigma_o(x)^\dagger + \eta_e(x+\mu)\otimes\rho_o(x)^\dagger   \right\} \right] \right)\, .\end{aligned}$$

## Two-point functions

This is a summary of the formulae used for the mesonic two-pt
functions.

Let $\Gamma$ and $\Gamma^\prime$ be two generic matrices in the Clifford
algebra, we define the two-pt function:

$$f_{\Gamma\Gamma^\prime}(t) = \sum_{\bf x}\langle \bar\psi({\bf x},t) \Gamma
\psi({\bf x},t) \bar\psi(0) \Gamma^\prime \psi(0) \rangle$$ 

Performing the Wick contractions yields: 

$$\begin{aligned}\langle \bar\psi({\bf x},t) \Gamma\psi({\bf x},t) \bar\psi(0) \Gamma^\prime \psi(0) \rangle =&- \mathrm{tr} \left[ \Gamma S(x-y) \Gamma^\prime S(y-x) \right]  \\
=& - \mathrm{tr} \left[ \Gamma S(x-y) \Gamma^\prime \gamma_5 S^\dagger(x-y) \gamma_5 \right] \end{aligned}$$ 

In practice we invert the Hermitean Dirac operator $\gamma_5 D$ by solving the equation:

$$Q_{AB}(x-y) \eta^{\bar A,x_0}_B(y) = \delta_{A,\bar A} \delta_{x,x_0}$$

where $A=\{a,\alpha\}$ is a collective index for colour and spin, and
$\bar A$, $x_0$ are the position of the source for the inverter.

Using the field $\eta$ that we obtain from the inverter, the correlator
above becomes:

$$\langle \ldots \rangle = - \tilde \Gamma_{AB} \eta^{C,y}_B(x)
\tilde \Gamma^\prime_{CD} \eta^{D,y}_A(x)^*$$ 

where

$\tilde \Gamma= \gamma_5 \Gamma$, and $\tilde \Gamma^\prime =
\gamma_5 \Gamma^\prime$.

## Hasenbusch acceleration

Let us summarize the Hasenbusch trick (for two flavours)

$$\mathcal{H}_F = \phi^\dagger ( Q_m^2 )^{-1} \phi \,$$ 

where $Q_m =\gamma_5 D_m$ is the hermitian Dirac operator. After integration over the pseudofermions it gives the determinant:

$$\det{ Q_m ^2 } = \det{D_m^{\dagger} D_m}$$

The Hasenbusch trick can be rewritten in the following form :

$$\det{ Q_m ^2 }  = \det{W_- W_+} \det{\frac{ Q_m^2}{W_- W_+}}$$

Where $W_{\pm}$ can be chosen arbitrarily as long as the determinant is
well defined. We discuss in the next subsections various choices of
$W_{\pm}$.

In any case the two term can be evaluated independently, and we have:

$$\mathcal{H}_{F_1} =   \phi_1^\dagger ( W_- W_+ )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger Q_m^{-1} W_- W_+ Q_m^{-1} \phi_2$$

This can be combined with even-odd preconditioning.

### Wilson Mass Shift

Assume 

$$W_{+} = \left( D_m + \delta_m\right) ,\quad W_{-} =
W_{+}^\dagger=  \left( D^{\dagger}_m + \delta_m\right)$$

Note that, as written in a comment in the code, $W_+ Q_m^{-1} = (a D +
b ) D^{-1} \gamma_5$.

Then 

$$Q_m^{-1} \left( D^\dagger_m + \delta_m\right)  \left( D_m +
   \delta_m  \right) Q_m^{-1}   =  \left( \gamma_5 + \delta_m Q^{-1}
   \right)  \left( \gamma_5 +    \delta_m Q^{-1}\right)$$

The force can then be computed : 

$$\begin{aligned}
 \dot{\mathcal{H}_{F_2}} =&  - \delta_m \phi_2^\dagger \left[  \left( \gamma_5 + \delta_m Q^{-1}
   \right) \dot{Q^{-1}} + \dot{Q^{-1}} \left( \gamma_5 + \delta_m Q^{-1}
   \right)   \right] \phi_2 \\
=&  - \delta_m \phi_2^\dagger \left[  \left( \gamma_5 + \delta_m
    Q^{-1}\right) Q_m^{-1} \dot{Q} Q_m^{-1}  \right]\phi_2 + \rm{h.c}
  \end{aligned}$$

Note that the equation as now the standard form of the forces for the
HMC algorithm provided that:

$$X\equiv Q^{-1}\phi_2,\quad\textrm{and}\quad Y^{\dagger}=\phi_2^\dagger
(\gamma_5 + \delta_m Q^{-1}) Q_m^{-1}$$

From which we deduce

$$Y =   Q_m^{-1}(\gamma_5 + \delta_m Q^{-1})  \phi_2 =  D^{-1} ( \phi_2  +
\delta_m \gamma_5 X)$$

Which matches one comment in the the force_hmc.c file.

#### Even-Odd Preconditioning

Writing 

$$D_m = \begin{pmatrix}
4 + m & D_{eo} \\
 D_{oe} & 4+m \\
\end{pmatrix}$$ 

The determinant in the 2 flavour case can be written as follows: 

$$\det D_m^2 = \det Q^2 \propto  \det D^\dagger_{eo} D_{oe}$$

$$Q = \gamma_5 \begin{pmatrix} 
1 + 4m & M_{\rm{eo}} \\
M_{\rm{oe}} & 1 +4m\\
\end{pmatrix} \equiv \gamma_5 \begin{pmatrix} 
M_{\rm{ee}}  & M_{\rm{eo}} \\
M_{\rm{oe}} & M_{\rm{oo}}\\
\end{pmatrix}$$ 

Note that $M_{\rm{ee}}^{-1}$ can be computed:

$$M_{\rm{ee}}^{-1} = \frac{1}{1+4m}$$

Now we can conveniently rewrite 

$$Q_{\pm} = \begin{pmatrix} 
\gamma_5 M_{\rm{ee}}  & 0 \\ \gamma_5 M_{\rm{oe}} & 1\\
\end{pmatrix} \begin{pmatrix} 
1  &  \left(M_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \\0 & \gamma_5
  \left(M^{\rm{oo}} - \frac{1}{4+m}M_{\rm{oe}} M_{\rm{eo}} \right)\\
\end{pmatrix}$$

From the last equation we deduce that:

$$\det{Q} = \det{\gamma_5 M^{\rm{ee}}} \det{\gamma_5  \left(M_{\rm{oo}}
    -\frac{1}{4+m} M_{\rm{oe}}     M_{\rm{eo}} \right)} \propto
\det{\gamma_5  \left((4+m) M_{\rm{oo}}
    - M_{\rm{oe}}     M_{\rm{eo}} \right)}$$

Note that the first determinant is a constant that could be computed.

In the following we will denote 

$$\hat{Q}_{m,eo} \equiv \gamma_5
  \left((4+m)^2 - M_{\rm{oe}}  M_{\rm{eo}} \right)$$ 

where $\hat{Q}_{m,eo}$ is defined on the odd sites of the lattice.

Now defining 

$$W_{+} = D_{m+\delta m} ,\quad W_{-} =
W_{+}^\dagger=   D^{\dagger}_{m+\delta m}$$

$$\begin{aligned}
\det{ Q_m (W_-  W_+)^{-1}  Q_m} \propto \det{ Q_{m,eo} (\hat{D}_{m+\delta_m}
  \hat{D}_{m+\delta_m,eo}  )^{-1}  Q_{m,eo}}\end{aligned}$$

We thus have 

$$\begin{aligned}
\mathcal{H}_{F_1} = \phi_1^\dagger \left( \hat{D}_{m+\delta m,eo} \hat{D}_{m+\delta_m,eo} \right)^{-1}\phi_1\end{aligned}$$

and 

$$\begin{aligned}
\mathcal{H}_{F_2} = \phi_2^\dagger Q_{m,eo}^{-1} \hat{D}_{m+\delta_m,eo} \hat{D}_{m+\delta_m,eo} Q_{m,eo}^{-1}\phi_2 \end{aligned}$$

Note that as in the non-even-odd case this can be rewritten as:

$$\begin{aligned}
\mathcal{H}_{F_2} =\phi_2^{\dagger} (\gamma_5 + \delta_m (1 +
\delta_m (4+m) )Q_{m,eo}^{-1}   (\gamma_5 + \delta_m (1 +
\delta_m (4+m) )Q_{m,eo}^{-1}\phi_2 \end{aligned}$$

### Twisted Mass Shift

Assume

$$W_{+} = \left( Q_m + i \mu_2  \right) ,\quad W_{-} = W_+^{\dagger}=\left( Q_m - i \mu_2 \right)$$

Note that $W_- W_+ = Q_m ^2 + \mu_2^2$ and that $W_\pm^{\dagger}=
W_{\mp}$.

Instead of dealing with $\det{ Q_m (W_- 
 W_+)^{-1} Q_m }$, we consider the slightly more general case where the
determinant to evaluate is 

$$\begin{aligned}
\det{ (Q_m + i \mu_1) (W_- 
 W_+)^{-1}  (Q_m - i \mu_1)}\propto&\int D\phi_2 D\phi_2^\dagger e^{-\phi^\dagger_2 \Big(Q_+ (W_- 
 W_+)^{-1}  Q_- \Big)^{-1}\phi_2 }\big] \\
=& \int D\phi_2 D\phi_2^\dagger e^{-\phi_2Q_-^{-1} W_- W_+
  Q_+^{-1}  \phi_2 }\end{aligned}$$ 

The following formulae can then be used for the case of several hasenbusch masses. The case of the determinant $\det{ Q_m (W_- W_+)^{-1} Q_m }$ can be recovered by setting $\mu_1=0$ in the following equations.

We have: 

$$\begin{aligned}
(Q_m&-i\mu_1)^{-1} W_- W_+ (Q_m+i\mu_1)^{-1} \\
 =& ( 1 - i(\mu_2 - \mu_1)
(Q_m - i \mu_1)^{-1}) (1+ i(\mu_2 - \mu_1)(Q_m - i \mu_1)^{-1}) \\
=& 1+ i(\mu_2 - \mu_1) (Q_m+i\mu_1)^{-1}  - i(\mu_2 - \mu_1)
(Q_m - i \mu_1)^{-1} + (\mu_2- \mu_1)^2\big((Q_m + i \mu_1)(Q_m - i
\mu_1)\big)^{-1} \\
=& 1+ (\mu_2- \mu_1)^2\big(Q_m^2 + \mu_1^2\big)^{-1}+ i(\mu_2 -
\mu_1) (Q_m^2 +\mu_1^2)^{-1}  (Q_m-i\mu_1) \\
\end{aligned}$$

$$\begin{aligned}
\qquad -i(\mu_2 &- \mu_1) (Q_m^2 +  \mu_1^2)^{-1} (Q_m+i\mu_1) \\
=& 1 +(\mu_2- \mu_1) \big(Q_m^2 + \mu_1^2\big)^{-1} \big( (\mu_2- \mu_1)
+ 2 \mu_1 \big) \\
=& 1 +(\mu_2^2- \mu_1^2) \big(Q_m^2 + \mu_1^2\big)^{-1} \end{aligned}$$

The force can then be computed: (global sign and factor $i$ have to be
checked) 

$$\begin{aligned}
\dot{\mathcal{H}_{F_2}} =& i(\mu_2-\mu_1) \phi_2^{\dagger} \Big[ ( 1 - i(\mu_2 - \mu_1) 
  (Q_m - i \mu_1)^{-1}) \dot{(Q_m+i\mu_1)^{-1}} 
 - \dot{(Q_m - i    \mu_1)^{-1}} (1+ i (\mu_2 - \mu_1) (Q_m+i \mu_1)^{-1})\Big] \phi_2 \\
=&  i(\mu_2-\mu_1) \phi_2^{\dagger} \left[  ( 1 - i(\mu_2 - \mu_1) 
  (Q_m - i\mu_1)^{-1}) ( Q_m+i\mu_1)^{-1} \dot{Q_m} (Q_m+i\mu_1)^{-1}
\right] \phi_2 +\rm{h.c}\end{aligned}$$

$$X\equiv (Q_m+i\mu_1)^{-1}\phi_2,\quad\textrm{and}\quad Y^{\dagger}=i\phi_2^\dagger
(1 -i (\mu_2-\mu_1) (Q-i\mu_1)^{-1}) (Q_m+i\mu_1)^{-1}$$

From which we deduce 

$$\begin{aligned}
Y =&   -i (Q_m - i \mu_1)^{-1}(1 +  i(\mu_2-\mu_1)(Q+i\mu_1)^{-1})
\phi_2 \\ =&-i (Q_m - i \mu_1)^{-1}  ( \phi_2  + i(\mu_2-\mu_1)  X)  \end{aligned}$$

Note that in the particular case where $\mu_1=0$,

$$Q_m^{-1} W_- W_+ Q_m^{-1} = ( 1 - i\mu_2Q_m^{-1}) (1+ i\mu_2Q_m^{-1}))
= 1 + \mu_2^2 Q_m^{-2}$$

Which leads to 

$$\begin{aligned}
\dot{\mathcal{H}_{F_2}}  &=& \mu_2^2 \phi_2^{\dagger} \dot{Q_m^{-2}} \phi_2\end{aligned}$$

Note also that the forces are explicitly proportional to $\mu_2^2$.

#### Even-Odd Preconditioning

Note that we have : $\widetilde{\mu} \equiv 2 \kappa \mu$.

$$Q_{\pm}= \gamma_5 \begin{pmatrix} 
1 \pm i \widetilde{\mu} \gamma_5 & M_{\rm{eo}} \\
M_{\rm{oe}} & 1 \pm i \widetilde{\mu} \gamma_5\\
\end{pmatrix} \equiv \gamma_5 \begin{pmatrix} 
M^{\pm}_{\rm{ee}}  & M_{\rm{eo}} \\
M_{\rm{oe}} & M^{\pm}_{\rm{oo}}\\
\end{pmatrix}$$ 

Note that $M_{\rm{ee}}^{-1}$ can be computed:

$$M_{\rm{ee}}^{-1} = ( 1 \pm i\widetilde{\mu} \gamma_5)^{-1} = \frac{1\mp
i \widetilde{\mu}\gamma_5 }{ 1 + \widetilde{\mu}^2}$$

Now we can conveniently rewrite 

$$Q_{\pm} =  \begin{pmatrix} \gamma_5 M^{\pm}_{\rm{ee}}  & 0 \\ \gamma_5 M_{\rm{oe}} & 1\\
\end{pmatrix} \begin{pmatrix} 
1  &  \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \\0 & \gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)\\
\end{pmatrix}$$

From the last equation we deduce that:

$$\det{Q_{\pm}} = \det{\gamma_5 M^{\pm}_{\rm{ee}} } \det{\gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)}$$

Note that the first determinant is a constant that could be computed.

In the following we will denote 

$$\hat{Q}_{\pm} \equiv \gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)$$ 

where $\hat{Q}_{\pm}$ is defined on the odd sites of the lattice.

We thus have

$$\det{Q_+ Q_-} = \det{Q_+}\det{Q_{-}}\propto\det{\hat{Q}_+ \hat{Q}_-}$$

and we thus get the following Hamiltonian:

$$\mathcal{H}_{F_1} = \phi_1^{\dagger} \left(\hat{Q}_+
    \hat{Q}_-\right)^{-1} \phi_1$$

The corresponding force then reads :

$$\dot{\mathcal{H}_{F_1}} = - \phi_{0}^\dagger\left(    \hat{Q}_-^{-1}
  \hat{Q}_+^{-1}  \dot{\hat{Q}}_+  \hat{Q}_+^{-1}  + \hat{Q}_{-}^{-1}
  \dot{\hat{Q}}_{-}  \hat{Q}_-^{-1}   \hat{Q}_+^{-1}    \right)    \phi_0$$

Now using that $Q_{\pm} ^{\dagger} = Q_{\mp}$, the previous equation can
be written:

$$\dot{\mathcal{H}_{F_1}} = - \left(Y_{\rm{o}}^{\dagger} \dot{\hat{Q}}_{+} X_{\rm{o}}  + \rm{h.c}\right)$$

with

$$X_{\rm{o}}=\hat{Q}_+^{-1}\phi_0,\quad. Y_{\rm{o}}= \left(\hat{Q}_+\hat{Q}_-\right)^{-1} \phi_0,$$ 

where we have used that

$$\hat{Q}_\pm^{\dagger} = \hat{Q}_\mp.$$ 

Furthermore we have

$$\dot{\hat{Q}}_{\pm} =  \gamma_5 \left( -  \dot{M}_{\rm{oe}}
  \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} -  M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} \dot{M}_{\rm{eo}}\right)$$

Now noting that 

$$\dot{Q}_{\pm} =  \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}$$

We have 

$$\begin{aligned}
  Y ^{\dagger} \dot{Q} X  =& \begin{pmatrix} A^\dagger &
    B^\dagger \end{pmatrix} \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}  \begin{pmatrix} C \\  D \\\end{pmatrix} \\
=&  A^\dagger \gamma_5 \dot{M}_{\rm{oe}}  C + B^\dagger \gamma_5 \dot{M}_{\rm{eo}} D\end{aligned}$$

Now chosing $A^\dagger= Y^\dagger_0$,
$C=\left(M^{+}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} X_0$, $B^\dagger=Y_0^\dagger \gamma_5 M_{\rm{oe}}\left(M^{+}_{\rm{ee}}\right)^{-1} \gamma_5$, and $D= X_0$ allows to
write: 

$$\dot{\mathcal{H}_{F_1}} =   Y^\dagger \dot{Q} X + \rm{h.c}$$ 

with

$$ X=\begin{pmatrix}\left(M^{+}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} X_0
  \\ X_0 \end{pmatrix},\quad \rm{and} \quad Y=\begin{pmatrix}
   \left(M^{-}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} Y_0
  \\ Y_0 \end{pmatrix}$$ 

We have used that $\dot{Q}_+= \dot{Q}_{-}$ and

$$ M_{\rm{eo}}^{\dagger}=\gamma_5 M_{\rm{oe}} \gamma_5 $$

##### Determinant Ratio

We use that $\det{ Q_+ (W_- 
 W_+)^{-1} Q_-} = \det{W_+^{-1} Q_+ Q_- W_-^{-1}} \propto\det{\hat{W}_+^{-1} \hat{Q}_+ \hat{Q}_- \hat{W}_-^{-1}}$

We thus have to compute 

$$\begin{aligned}
\dot{\mathcal{H}_{F_2}}  =& \phi_2^\dagger\big[ \delta \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+  +  \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \delta \hat{W}_+ \\
+&  \hat{W}_- \delta \hat{Q}_-^{-1} \hat{Q}_+^{-1} \hat{W}_+ +
\hat{W}_- \hat{Q}_-^{-1} \delta \hat{Q}_+^{-1} \hat{W}_+   \big] \phi_2 \\
=& \phi_2^\dagger\big[ \dot{\hat{W}}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+  +  \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \delta \hat{W}_+ \\
-&  \hat{W}_- \hat{Q}_-^{-1} \dot{\hat{Q}}_- \hat{Q}_-^{-1} \hat{Q}_+^{-1} \hat{W}_+ -
\hat{W}_- \hat{Q}_-^{-1} \hat{Q}_+^{-1} \dot{\hat{Q}}_+ \hat{Q}_+^{-1} \hat{W}_+   \big] \phi_2 \\\end{aligned}$$

Now we introduce

$$X_W =  (\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+ \phi_2,  Y_W =
\hat{Q}_+^{-1} \hat{W}_+\phi_2 = \hat{Q}_- X_W$$

such that 

$$\begin{aligned}
\dot{\mathcal{H}_{F_2}} &=& \phi_2^\dagger  \dot{\hat{W}}_- X_W +
X_W^\dagger  \delta \hat{W}_+ \phi_2 \\
&-& Y_W^\dagger \dot{\hat{Q}}_- X_W - X_W^\dagger   \dot{\hat{Q}}_+  Y_W\end{aligned}$$

Now recalling that 

$$\begin{aligned}
\dot{\hat{Q}}_{\pm} =&  -\gamma_5 \left(   \dot{M}_{\rm{oe}}
  \left(1\pm i\mu_1 \gamma_5\right)^{-1} M_{\rm{eo}}   + M_{\rm{oe}}
    \left( 1\pm i\mu_1 \gamma_5\right)^{-1} \dot{M}_{\rm{eo}}\right)  \\
\dot{\hat{W}}_{\pm} =&  -\gamma_5 \left(   \dot{M}_{\rm{oe}}
  \left(1\pm i\mu_2 \gamma_5\right)^{-1} M_{\rm{eo}} +  M_{\rm{oe}}
    \left( 1\pm i\mu_2 \gamma_5\right)^{-1} \dot{M}_{\rm{eo}}\right)\end{aligned}$$

Now can write the last expression in terms of $\dot{Q} \equiv
\dot{Q}_{\pm}$. 

$$\begin{aligned}
\dot{\mathcal{H}_{F_2}} =& Y_1^\dagger \dot{Q} X_1 + X_1^\dagger
\dot{Q} Y_1 - X_2^\dagger \dot{Q} Y_2  - Y_2^\dagger \dot{Q} X_2 \\
=& 2~ \rm{Re}\Big[ Y_1^\dagger \dot{Q} X_1 - Y_2^\dagger \dot{Q}
X_2 \Big],\end{aligned}$$ 

with 

$$\begin{aligned}
Y_1 =& \begin{pmatrix}  (1+i\mu_1\gamma_5)^{-1} M_{\rm{eo}} Y_W \\
  Y_W \end{pmatrix},\quad Y_2  =  \begin{pmatrix}  (1+i\mu_2\gamma_5)^{-1} M_{\rm{eo}} \phi_2 \\
  \phi_2 \end{pmatrix},\\
X_{1,2} =&  \begin{pmatrix}  (1-i\mu_{1,2}\gamma_5)^{-1} M_{\rm{eo}} X_W \\
  X_W \end{pmatrix}, \dot{Q} \equiv  \dot{Q}_{\pm} =  \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}\end{aligned}$$

##### Twisted Wilson-Dirac Operator

Instead of applying the eo preconditionning to the twisted mass operator
we can use the wilson dirac eo operator and do a different splitting.

We define:

$$Q_{\rm{eo}} \equiv \gamma_5 \big( (4+m)^2 - D_{eo} D_{oe}\big)$$

Now we split the determinant as follows:

$$\det{ Q_{\rm{eo}} ^2 }  = \det{W_- W_+} \det{\frac{  Q_{\rm{eo}} ^2
  }{W_- W_+}}$$

And we choose 

$$W_{\pm} = Q_{\rm{eo}} \pm i \mu$$

The corresponding Hamiltonian read:

$$\mathcal{H}_{F_1} = \phi_1^\dagger ( W_- W_+ )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger Q_{\rm{eo}}^{-1} W_- W_+ Q_{\rm{eo}}^{-1} \phi_2$$

Since the operator are now very similar to the non even-odd case, we can reuse
some formulae. In particular, we can rewrite the Hamiltonian as follows:

$$\mathcal{H}_{F_1} =  \phi_1^\dagger (Q_{\rm{eo}}  + \mu^2  )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger \big( 1  + \mu^2 Q_{\rm{eo}}^{-1}\big)  \phi_2$$

From which we have the following forces: 

$$\begin{aligned}
\dot{\mathcal{H}}_{F_1} =&   \phi_1^\dagger  W_+^{-1} \delta
W_-^{-1}\phi_1 + \rm{h.c} \\
=&  \phi_1^\dagger  W_+^{-1}  W_-^{-1} \dot{Q}_{\rm{eo}} W_-^{-1} \phi_1+ \rm{h.c}\end{aligned}$$

Now we want to rewrite the last equation as a function of

$$\dot{Q} = \gamma_5  \begin{pmatrix}  0 & \dot{D}_{\rm{eo}} \\
  \dot{D}_{\rm{oe}} & 0  \end{pmatrix}$$

$$X_{\rm{o}}=W_+^{-1}\phi_1,\quad. Y_{\rm{o}}= \left(W_+
    W_-\right)^{-1} \phi_1,$$ 

where we have used that

$$W_\pm^{\dagger} = W_\mp.$$ 

Furthermore we have

$$\dot{Q}_{\rm{eo}} =  -\gamma_5 \left(  \dot{M}_{\rm{oe}}
   M_{\rm{eo}} +  M_{\rm{oe}} \dot{M}_{\rm{eo}}\right)$$

noting that

$$\begin{aligned}
  Y ^{\dagger} \dot{Q} X  &=& \begin{pmatrix} A^\dagger &
    B^\dagger \end{pmatrix} \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}  \begin{pmatrix} C \\  D \end{pmatrix} \\
&=&  A^\dagger \gamma_5 \dot{M}_{\rm{oe}}  C + B^\dagger \gamma_5 \dot{M}_{\rm{eo}} D\end{aligned}$$

and chosing $A^\dagger= Y^\dagger_0$, $C =  M_{\rm{eo}} X_0$,
$B^\dagger=
 Y_0^\dagger \gamma_5 M_{\rm{oe}}$, and $D= X_0$ allows us to write: 
 
$$\dot{\mathcal{H}_{F_1}} = -Y^\dagger \dot{Q} X + \rm{h.c}$$ 
 
with 
 
$$X=\begin{pmatrix} M_{\rm{eo}} X_0 \\ X_0 \end{pmatrix},\quad \rm{and} \quad Y=\begin{pmatrix} M_{\rm{eo}} Y_0 \\ Y_0 \end{pmatrix}$$ 

We have used that $M_{\rm{eo}}^{\dagger}=\gamma_5 M_{\rm{oe}} \gamma_5$.

Similarly for the second Hamiltonian we get: 

$$\begin{aligned} \dot{\mathcal{H}}_{F_2} = \mu_2\phi_2^\dagger  \dot{Q_{\rm{eo}}^{-1}}  \phi_2\end{aligned}$$

which is exactly the force that appears in case of a pure Wilson-Dirac even-ddd preconditioned operator up to a multiplicative factor.

## Clover Term

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
 Q_{\mu\nu}(x)
 &= U_\mu(x)U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x) \\
 &+ U_\nu(x)U_\mu^\dagger(x-\hat{\mu}+\hat{\nu})U_\nu^\dagger(x-\hat{\mu})U_\mu(x-\hat{\mu}) \\
 &+ U_\mu^\dagger(x-\hat{\mu})U_\nu^\dagger(x-\hat{\mu}-\hat{\nu})U_\mu(x-\hat{\mu}-\hat{\nu})U_\nu(x-\hat{\nu}) \\
 &+ U_\nu^\dagger(x-\hat{\nu})U_\mu(x-\hat{\nu})U_\nu(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x)
\end{aligned} $$

Because $Q_{\mu\nu}^\dagger=Q_{\nu\mu}$ we have $F_{\mu\nu}=-F_{\nu\mu}$.
For this reason we can change the sum over $\mu,\nu$ in Eq. {eq}`eq:clover` to a sum over $\mu<\nu$ and a factor of two.

$$ D_{sw} = -\frac{c_{sw}}{2}\sum_x\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$

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

### Pseudofermion Forces

For the forces we use the following short-hand notation for the derivative with respect to the link variables.

$$ \dot{S} = \partial_{x,\mu}^a S $$

To calculate the pseudofermion forces let us write down the action as

$$ S = \phi^\dagger(H^{-2})\phi, $$

where $H=\gamma^5D$ is the Hermitian Dirac operator.
When differentiating the action we obtain

$$ \dot{S} = -2\mathrm{Re}~\xi^\dagger\dot{H}\eta, $$(eq:dotS)

with the definitions

$$ \begin{aligned}
 \eta &= H^{-2}\phi, \\
 \xi &= H\eta.
\end{aligned} $$

#### Forces

Here we will only consider the forces from the clover term and not the hopping term.
The clover part of the Dirac operator is given by

$$ H_{sw} = -\frac{c_{sw}}{2}\sum_{\mu<\nu}\bar{\sigma}_{\mu\nu}F_{\mu\nu}(x) $$(eq:Hsw)

When inserting Eq. {eq}`eq:dotS` we obtain

$$ \dot{S} = c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{F}_{\mu\nu}\eta). $$

From the definition of $F_{\mu\nu}$ it follows that

$$ \begin{aligned}
 \dot{S} &= \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \xi^\dagger\bar{\sigma}_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger\eta), \\
         &= \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \eta^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\xi).
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

#### Logarithmic Forces

In the case of even-odd preconditioning (see the next section) the action of the small determinant $D_{oo}$ can be written as

$$ S_{sw} = -N_f\log\det D_{oo} = -N_f~\mathrm{tr}\log D_{oo} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\log D_{oo}(x) $$

The derivative is given by

$$ \dot{S} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\left[D_{oo}^{-1}(x)\dot{D}_{oo}(x)\right] $$

with $D_{oo}(x)$ given by

$$ D_{oo}(x) = 4+m_0-\frac{c_{sw}}{2}\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$

Both the determinant and the inverse of $D_{oo}(x)$ can be calculated from an LDL decomposition.
If we insert the above definition we obtain

$$ \dot{S} = \frac{N_fc_{sw}}{2}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{F}_{\mu\nu}(x)) $$

$$ \dot{S} = \frac{N_fc_{sw}}{2\cdot8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu} \mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{Q}_{\mu\nu}(x) + D_{oo}^{-1}(x)\sigma_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger(x)) $$

Since $D_{oo}^{-1}$ is Hermitian we can write the result as two times the real part.
To simplify the result we define $X_{\mu\nu}(x)=D_{oo}^{-1}(x)\sigma_{\mu\nu}$ such that

$$ \dot{S} = \frac{N_fc_{sw}}{8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{Re}~\mathrm{tr}[X_{\mu\nu}(x)\dot{Q}_{\mu\nu}(x)] $$

This is equivalent to Eq. {eq}`eq:force` except from the factor $N_f$ and the definition of $X_{\mu\nu}(x)$.
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
  D_{11} & D_{12} \\
  D_{21} & D_{22} \\
 \end{pmatrix} $$

$$ D_{-}^{-1} =
 \begin{pmatrix}
  D_{33} & D_{34} \\
  D_{43} & D_{44} \\
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


### Even-odd Preconditioning

#### Method 1

We can write the determinant as

$$ \det D = \det(D_{oo})\det(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe}) $$

Use the notation

$$ \begin{aligned}
 Q_{oo} &= \gamma_5D_{oo} \\
 Q &= \gamma_5(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe})
\end{aligned} $$

The action is

$$ S = S_1 + S_2 = \phi_1^\dagger Q_{oo}^{-2}\phi_1 + \phi_2^\dagger Q^{-2}\phi_2 $$

##### Forces for $\phi_1$-term

The derivative is

$$ \dot{S}_1 = -2\mathrm{Re}\left[\phi_1^\dagger(Q_{oo}^{-2}Q_{oo}\dot{Q}_{oo}Q_{oo}^{-2})\phi_1\right] $$

and we can write it as

$$ \dot{S}_1 = -2\mathrm{Re}\left[\xi^\dagger\dot{Q}_{oo}\eta\right] $$

with

$$ \begin{aligned}
 \eta &= Q_{oo}^{-2}\phi_1, \\
 \xi &= Q_{oo}\eta.
\end{aligned} $$

##### Forces for $\phi_2$-term

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

#### Method 2

The action of $D_{oo}$ can also be expressed directly as the logarithm of the determinant.

$$ S = -2\log\det D_{oo} + \phi^\dagger Q^{-2}\phi $$

This is the approach implemented in the code.

### LDL factorization

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

#### LDL Decomposition

Calculates the LDL decomposition $A=LDL^\dagger$ in-place.
After the decomposition, the lower triangular part of $A$ is $L$ and the diagonal is $D$.

```
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

#### Forward substitution

Calculates $x=L^{-1}b$.

```
do i=0, N-1
    x_i = b_i
    do k=0, i-1
        x_i = x_i - A_ik * x_k
    enddo
enddo
```

#### Backward substitution with diagonal

Calculates $x=(L^\dagger)^{-1}D^{-1}x$.

```
do i=N-1, 0
    x_i = x_i/A_ii
    do k=i+1, N-1
        x_i = x_i - conj(A_ki) * x_k
    enddo
enddo
```

#### Full inversion

This algorithm calculates the inverse $B=A^{-1}$ from the LDL decomposition.
Because the inverse is Hermitian we only calculate the lower triangular part.

```
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
	    B_ji = conj(L_kj) * B_ki
	enddo
    enddo
enddo
```


## Exponential Clover Term

The exponential version of the clover term (including mass term) can be
written as

$$D_{sw} = \sum_x  (4+m_0) \exp ( A(x) ), \text{ with } A(x) =  -\frac{c_{sw}}{4(4+m_0)}\sum_{\mu,\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$ (eq:clover_exp)

where $\sigma_{\mu\nu}$ is again defined by

$$\sigma_{\mu\nu} = \frac{1}{2}[\gamma_\mu,\gamma_\nu].$$

As for the clover term above, we can simplify the sum over $\mu\nu$ to a sum over $\mu < \nu$ and introduce a factor of two. We define

$$A(x) = -\frac{c_{sw}}{2(4+m_0)}\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) $$ (eq:clover2)

The quantity $\sigma_{\mu\nu}F_{\mu\nu}$ is Hermitian and block
diagonal. It can be written as 

$$\begin{aligned} A(x)=-\frac{c_{sw}}{2(4+m_0)} \sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu} = 
 \begin{pmatrix}
 a & b & 0 & 0 \\
 b^\dagger & -a & 0 & 0 \\
 0  & 0 & c & d \\
 0 & 0 & d^\dagger & -c
 \end{pmatrix} \equiv \begin{pmatrix}
 A^+ & 0 \\ 
 0 & A^- 
 \end{pmatrix}, \end{aligned}$$ (eq:blocks)

where $A^\pm$ are $2 \times 2$ matrices in spin space and $a,b,c,d$ are $N_F \times N_F$.

This formulation of $A(x)$ as a block matrix will be useful for the
exponentiation.

### Evaluation of the operator

The evaluation of the exponential of $A(x)$ can be split as:

$$\exp A(x)  = \begin{pmatrix}
\exp A^+ & 0 \\
0& \exp A^-
\end{pmatrix}$$ 

and so, the problem is reduced to the exponential of two
$(2 N_F) \times  (2 N_F)$ matrices. The evaluation can be performed in
two ways.

1.  Using the Taylor expansion:

    $$\exp(A^\pm) = \sum_{k=0}^N  \frac{1}{k!} (A^\pm)^k.$$

2.  Using the Horner scheme:

$$\exp(A^\pm) = \sum_{k=0}^{\dim A^\pm -1} b_k(A^\pm)  (A^\pm)^k, $$(eq:exphorner)

where $b_k$ are computed recursively as follows. We start with

$$\begin{aligned} q_{N,0} = 1/N!, q_{N,1 \cdots (\dim A^\pm)-1} = 0.\end{aligned}$$

Then, the recursion proceeds: 

$$q_{n,0} = - p_0 q_{n+1, \dim A^\pm-1}  + 1/n!,$$
$$q_{n,i} = - p_i q_{n+1, \dim A^\pm-1}  + q_{n+1,i-1},$$
$$\text{ with } i=1 \cdots (\dim A^\pm) -1,$$(eq:horner) 

where $p_i$ represent the coefficients of the characteristic polynomial of the matrix $A^\pm$:

$$P(A^\pm) = \sum_{n=0}^{\dim A\pm} p_n (A^\pm)^n.$$ 

For instance, the characteristic polynomial of a $4 \times 4$ traceless matrix has the following coefficients:

$$p_0=\frac{1}{8 } \left(\mathrm{tr}A^2\right)^2 - \frac{1}{4} \mathrm{tr}A^4 ,\ p_1 = -\frac{1}{3}\mathrm{tr}A^3, \ p_2= -\frac{1}{2}\mathrm{tr}A^2, \ p_3=0, \ p_3=1.$$

Finally, the coefficients of eq. {eq}`eq:exphorner` are $b_k  =q_{0,k}$.

The Horner scheme method is currently implemented only for $SU(2)$ and
$SU(3)$ with fundamental fermions.

### Pseudofermion Forces

For the forces we use the following short-hand notation for the
derivative with respect to the link variables.

$$\dot{S} = \partial_{x,\mu}^a S.$$ 

To calculate the pseudofermion
forces let us write down the action as 

$$S = \phi^\dagger(H^{-2})\phi,$$

where $H=\gamma^5D$ is the Hermitian Dirac operator. When differentiating the action we obtain

$$\dot{S} = -2\mathrm{Re}~\xi^\dagger\dot{H}\eta, $$ (eq:dotS)

with the definitions 

$$\begin{aligned}
 \eta &= H^{-2}\phi, \\
 \xi &= H\eta.\end{aligned}$$
 
 #### Forces

Here we will only consider the forces from the clover term. For the
exponential version of the clover term, the implementation is very
similar to the traditional clover term.

The clover part of the Dirac operator is given by

$$H_{sw} = (4+m_0) \gamma_5 \exp \left(- \frac{c_{sw}}{2(4+m_0)}  \sum_{\mu<\nu}{\sigma}_{\mu\nu}F_{\mu\nu}(x)  \right) = (4+m) \gamma_5 \exp A(x).$$ (eq:Hsw)

An optimized way of calculating the derivative is provided by the double
Horner scheme. The basic idea is that the derivative of a matrix can be
expressed as:

$$d e^A = \sum_{k=0}^{\dim A -1} \sum_{l=0}^{\dim A-1} C_{kl}(A) A^l dA A^k, $$ (eq:dexphorner)

where the $C_{kl}$ coefficients depend on the matrix $A$, similarly to
the ones eq. {eq}`eq:exphorner`. They are calculated performing first the
iteration in eq. {eq}`horner`, and then repeating the iteration process on the result of the first iteration. For compactness, we shall omit the limits of the sum henceforth. When inserting eq. {eq}`eq:Hsw` in eq. {eq}`eq:dotS`, and using eq. {eq}`eq:dexphorner` we obtain 

$$\dot{S} = c_{sw} \sum_{k}  \sum_{\mu<\nu}\mathrm{Re}(  \xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{F}_{\mu\nu}\eta_k),$$
with $$\xi_k = \sum_l \begin{pmatrix}
C^+_{kl} (A^+)^l \xi^+ \\
C^-_{kl} (A^-)^l \xi^-
\end{pmatrix}, \text{ and } 
\eta_k =  \begin{pmatrix}
 (A^+)^k \eta^+ \\
(A^-)^k \eta^-
\end{pmatrix}.$$

From the definition of $F_{\mu\nu}$ it follows that 

$$\begin{aligned}
 \dot{S} &=  \frac{1}{8}c_{sw}\sum_k \sum_{\mu<\nu}\mathrm{Re}(\xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta_k + \xi_k^\dagger\bar{\sigma}_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger\eta_k), \\
 &=
 \frac{1}{8}c_{sw}\sum_k \sum_{\mu<\nu}\mathrm{Re}(\xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta_k + \eta_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\xi_k).\end{aligned}$$

This can in be written as

$$\dot{S} = \frac{1}{8}c_{sw}\sum_{\mu<\nu} \mathrm{Re}~\mathrm{tr}\left[\dot{Q}_{\mu\nu} \sum_k\left\{\bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x)\right\}\right]$$

As for the clover term above we need to calculate now

$$\dot{S} = \frac{1}{8}c_{sw}\mathrm{Re}~\mathrm{tr}[\dot{Q}_{\mu\nu}(x)X_{\mu\nu}(x)] $$(eq:force)

now with

$$X_{\mu\nu}(x) = \sum_k \bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x).$$

The total force can now be expressed as in the clover term above. 

### Even-odd preconditioning

Even-odd preconditioning is particularly simple for the exponential
case, since the force coming from the little determinant vanished. This
can be seen because of the fact that:

$$\det D_{oo}  = \exp(\log \det D_{oo} ) =\exp( \mathrm{tr}\log D_{oo}) =  1,$$

and so it is a constant term in the action that does not contribute to
the force.

### Implementation of $X_{\mu\nu}$ using Taylor series 

In the current version of the code, the horner scheme is only implemeted
for $SU(2)$ and $SU(3)$ with fundamental fermions. For other theories, a
less efficient, but more flexible, alternative is used. For this,
we use the Taylor series:

$$dA = \sum_{k=0}^N \sum_{l=0}^{N-k} \frac{1}{(k+l+1)!} A^{k} dA A^{l},$$

with $N$ sufficiently large. The implementation changes only in the definition of $X_{\mu\nu}$:

$$X_{\mu\nu}(x) = \sum_{k=0}^N \bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x),$$

where now: 

$$\xi_k = \sum_l \frac{1}{(k+l+1)!} \begin{pmatrix}
(A^+)^l \xi^+ \\
 (A^-)^l \xi^-
\end{pmatrix}, \text{ and } 
\eta_k =  \begin{pmatrix}
 (A^+)^k \eta^+ \\
(A^-)^k \eta^-
\end{pmatrix}.$$

## Stout smearing 

The implementation follows \[hep-lat/0311018\] closely. We define the
smeared links as 

$$U'_\mu(x) = e^{Q_\mu(x)}U_\mu(x)$$

where $\Sigma_\mu(x)$ is an element of the Lie algebra, defined via the
projection 

$$Q_\mu(x) = \mathcal{P}\{\Omega_\mu(x)\}.$$

The projection operator is not unique, but the most common choice is

$$\mathcal{P}(\Omega) = \frac{1}{2}(\Omega-\Omega^\dagger) - \frac{1}{2N}\mathrm{tr}(\Omega-\Omega^\dagger).$$

However, in a setup with mixed representations, it is convenient to use
the following two-step procedure for the projection. This allows us to
project matrices from different representations onto the generators of
the fundamental representation. 

$$\begin{aligned} A_\mu^a(x) &= -\frac{1}{T_f}\mathrm{tr}[iT^a_R\Omega_\mu(x)] \\ Q_\mu(x) &= iT^a_F A_\mu^a(x)\end{aligned}$$

The matrix $\Omega_\mu(x)$ is defined as

$$\Omega_\mu(x) = U_\mu(x)V_\mu(x)$$ 

$$V_\mu(x) = \sum_{\nu\neq\mu}
 \rho_{\mu\nu}\left(
 U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x) +
 U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})
 \right)$$

For the force calculation we use the chain rule.

$$\frac{dS}{dU} = \frac{dS}{dU'}\frac{dU'}{dU}$$

The first derivative on the right-hand side is the usual force
$\Sigma_\mu'(x)$ evaluated using the smeared links. The second term is
the derivative of the smeared links with respect to the fundamental
links. This can be written in the following way, because the derivative
of the action is surrounded by a trace.

$$\frac{dS}{dU} = e^Q\Sigma' + \frac{d\Omega}{dU}\mathcal{P}(X)$$

When using a Taylor expansion to define the exponential function, we can
use the following definition of $X$.

$$X = \sum_{n=0}^\infty\sum_{k=0}^n Q^kU\Sigma' Q^{n-k}$$

The derivative of the $\Omega$ matrix is the last missing piece. Define
$\Lambda=\mathcal{P}(X)$ and consider

$$\frac{d}{d U_\mu(x)} U_\mu(x)V_\mu(x)\Lambda_\mu(x)$$

Here we have a sum over $\nu\neq\mu$. There are eight contributions to
the above derivative. 

$$\rho_{\mu\nu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\Lambda_\mu(x)$$ 

$$\rho_{\nu\mu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\Lambda_\nu(x-\hat{\nu})U_\nu(x-\hat{\nu})$$

$$\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\Lambda_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})$$

$$\rho_{\nu\mu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\Lambda_\nu^\dagger(x)$$


$$\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})\Lambda_\mu(x)$$

$$\rho_{\nu\mu}U_\nu^\dagger(x-\hat{\nu}+\hat{\mu})\Lambda_\nu^\dagger(x-\hat{\nu}+\hat{\mu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})$$

$$\rho_{\mu\nu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})\Lambda_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)$$

$$\rho_{\nu\mu}\Lambda_\nu(x+\hat{\mu})U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)$$

This can be simplified because several products appear more than once
and we can use $\Lambda^\dagger=-\Lambda$ to remove some of the
Hermitian conjugates. In the following we also assume that
$\rho_{\mu\nu}=\rho_{\nu\mu}$. 

$$\rho_{\mu\nu}W_2U_\nu^\dagger(x)\{\Lambda_\mu(x)-\Lambda_\nu(x)\}$$

$$\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\{\Lambda_\nu(x-\hat{\nu})-\Lambda_\mu(x-\hat{\nu})\}U_\nu(x-\hat{\nu})$$

$$\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})\{W_1\Lambda_\mu(x)-\Lambda_\nu(x+\hat{\mu}-\hat{\nu})W_1\}$$

$$\rho_{\mu\nu}\{\Lambda_\nu(x+\hat{\mu})W_2-W_2\Lambda_\mu(x+\hat{\nu})\}U_\nu^\dagger(x)$$

Here 

$$\begin{aligned}
 W_1 &= U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu}) \\
 W_2 &= U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})\end{aligned}$$

This brings us down to 13 multiplications.