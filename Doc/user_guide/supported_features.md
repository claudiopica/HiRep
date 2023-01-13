@page supported_features Supported Features
[TOC]
# Conventions {#conventions}

This section summarizes the main formulae that are used for implementing
the HMC for dynamical Wilson fermions in higher representations. The
Dirac operator is constructed following @cite Luscher:1996sc, but
using Hermitian generators 

\f{equation}{T^{a\dagger}=T^a.\f} 

For the fundamental
representation, the normalization of the generators is such that:

\f{equation}{\mathrm{tr}\, \left(T^a T^b \right) = \frac12 \delta^{ab}.\f}

For a generic representation \f$ R\f$, we define: 

\f{equation}{\mathrm{tr }_R \left(T^a T^b \right) = T_R \delta^{ab},\f}

\f{equation}{\sum_a \left(T^a T^a \right)_{AB} = C_2(R) \delta_{AB},\f}

which implies

\f{equation}{T_R = \frac{1}{N^2-1} C_2(R) d_R\,,\f} 

where \f$ d_R\f$ is the dimension of the representation \f$ R\f$. The relevant group factors may be computed from the Young tableaux of the representation of \f$ SU(N)\f$ using

\f{equation}{C_2(R) =\frac{1}{2}\left(nN+ \sum_{i=1}^{m} n_i \left( n_i+1-2i
\right) - \frac{n^2}{N}\right)\, .\f} 

Here \f$ n\f$ is the number of boxes in the diagram, \f$ i\f$ ranges over the rows of the Young tableau, \f$ m\f$ is the number of rows, and \f$ n_i\f$ is the number of boxes in the \f$ i\f$-th row.

|    R        |        \f$ d_R\f$        |       \f$ T_R\f$     |           \f$ C_2(R)\f$          |
|-------------|---------------------|-----------------|-----------------------------|
| fundamental (fund)|         \f$ N\f$         |     \f$ \frac12\f$   |     \f$ \frac{N^2-1}{2 N}\f$     |
| adjoint (adj) |      \f$ N^2-1\f$        |       \f$ N\f$       |            \f$ N\f$              |
| two-index symmetric (2S) | \f$ \frac{1}{2}N(N+1)\f$ | \f$ \frac{N+2}{2}\f$ |  \f$ C_2(f)\frac{2(N+2)}{N+1}\f$ |
| two-index antisymmetric (2AS) | \f$ \frac{1}{2}N(N-1)\f$ |  \f$ \frac{N-2}{2}\f$|  \f$ C_2(f)\frac{2(N-2)}{N-1}\f$ |


A generic element of the algebra is written as \f$ X=i X^a T^a\f$ and the
scalar product of two elements of the algebra is defined as

\f{equation}{(X,Y)= \mathrm{tr\ } \left(X^\dagger Y\right) = T_f X^a Y^a,\f}

\f{equation}{\Vert X \Vert^2 = \mathrm{tr } \left(X^\dagger X\right)
 = \sum_{ij} \left| X_{ij} \right|^2\, .\f}
 
### \f$ \gamma\f$ matrices

We use the chiral representation for the Dirac \f$ \gamma\f$ matrices where

\f{equation}{\gamma_\mu=\begin{pmatrix}0&e_\mu \\\ e_\mu^\dagger&0 \end{pmatrix}\, .\f} 

Then \f$ e_\mu\f$ are \f$ 2\times 2\f$ matrices given by \f$ e_0=-\mathbb{1}\f$, \f$ e_k=-i\sigma_k\f$ corresponding to

\f{equation}{\sigma_1=
\begin{pmatrix}
0&1\\\
1&0
\end{pmatrix},\,\,
\sigma_2=
\begin{pmatrix}
0&-i\\\
i&0
\end{pmatrix},\,\,
\sigma_3=
\begin{pmatrix}
1&0\\\
0&-1
\end{pmatrix}\, .\f} 

Finally

\f{equation}{\gamma_5=\gamma_0\gamma_1\gamma_2\gamma_3=
\begin{pmatrix}
1&0\\\
0&-1
\end{pmatrix}\, .\f}

# Representations

The hermitean generators \f$T^a_f\f$ for the fundamental representation used are of the form:
\f{equation}{
\begin{pmatrix} 
0&1&0&\dots\\
1&0&0&\dots\\
0&0&0&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, ,
\begin{pmatrix} 
0&i&0&\dots\\
-i&0&0&\dots\\
0&0&0&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, ,
\begin{pmatrix} 
1&0&0&\dots\\
0&1&0&\dots\\
0&0&-2&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, ,
\f}
normalized so that \f$T_f=1/2\f$. The generators for the other representations 
will be obtained in the following.

We first give the explicit form for the representation functions \f$R\f$ which map 
$U\rightarrow U^R$. We define for each representation an orthonormal base \f$e_R\f$ for 
the appropriate vector space of matrices. 

For the Adjoint representation we define the base \f$e_{Adj}\f$ for the \f$N\times N\f$ 
traceless hermitean matrices to be \f$e_{Adj}^a=T^a_f/\sqrt{T_f}\f$, \f$a=1,\dots,N^2-1\f$ 
(i.e. proportional to the generators of the fundamental representation and 
normalized to 1.)

For the two-index Symmetric representation the base \f$e^{(ij)}_{S}\f$, with \f$i\le j\f$, for 
the \f$N\times N\f$ symmetric matrices is given by:
\f{eqnarray}{
i\neq j \, ,\,\,\,&e^{(ij)}_S&=\frac{1}{\sqrt{2}}\begin{pmatrix} 
0&1&0&\dots\\
1&0&0&\dots\\
0&0&0&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, , \\
i=j \, ,\,\,\,&e^{(ii)}_S&=\begin{pmatrix} 
0&0&0&\dots\\
0&1&0&\dots\\
0&0&0&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, , 
\f}
where the non zero entries are at position \f$(i,j)\f$, etc.

For the two-index Antisymmetric representation the base \f$e^{(ij)}_{AS}\f$, with \f$i<j\f$, for 
the \f$N\times N\f$ symmetric matrices is given by:
\f{equation}{
e^{(ij)}_{AS}=\frac{1}{\sqrt{2}}\begin{pmatrix} 
0&1&0&\dots\\
-1&0&0&\dots\\
0&0&0&\dots\\
\dots&\dots&\dots&\dots
\end{pmatrix}\, , 
\f}
where, as above, the non zero entries are at position \f$(i,j)\f$.

The maps \f$R\f$ are explicitly given by:
\f{eqnarray}{
(R^{Adj} U)_{ab} &=& U^{Adj}_{ab} = \mathrm{tr\ }\left[ e^a_{Adj} U e^b_{Adj} U^\dagger\right]\,\, , a,b=1,\dots,N^2-1\, ,\\ 
(R^{S} U)_{(ij)(lk)} &=& U^{S}_{(ij)(lk)} = \mathrm{tr\ }\left[ (e^{(ij)}_{S})^\dagger U e^{(lk)}_S U^T\right]\,\, , i\le j, l\le k\, ,\\ 
(R^{A} U)_{(ij)(lk)} &=& U^{A}_{(ij)(lk)} = \mathrm{tr\ }\left[ (e^{(ij)}_{A})^\dagger U e^{(lk)}_A U^T\right]\,\, , i< j, l< k\, .
\f}

The generators \f$T_R^a\f$ used are defined as the image of the generators in the fundamental
under the differential of the maps \f$R\f$ defined above: \f$T^a_R = R_* T^a_f\f$.
Explicit expression can easily be worked out form the definition above.
The invariants \f$T_R\f$ and \f$C_2(R)\f$ for the generators defined in this way are given in 
Sect. [Conventions](#conventions)


# The Dirac operator

The massless Dirac operator is written as in @cite Luscher:1996sc

\f{equation}{
  D = \frac12 \left\{\gamma_\mu \left(\nabla_\mu + \nabla^*_\mu \right) - \nabla^*_\mu \nabla_\mu \right\}
\f} 

with 

\f{equation}{\nabla_\mu\phi(x) = U^R (x,\mu)\phi(x+\mu) - \phi(x)\f}
\f{equation}{\nabla_\mu^*\phi(x) = \phi(x) - U^R (x-\mu,\mu)^\dagger\phi(x-\mu)\f}

and therefore the action of the massive Dirac operator yields

\f{equation}{
  \begin{aligned}
 D_m \phi(x) =& (D+m) \phi(x)\\ 
=& - \frac12 \left\{ \left(1-\gamma_\mu\right) U^R(x,\mu) \phi(x+\mu) \right.
+
\left.\left(1+\gamma_\mu\right) U^R(x-\mu,\mu)^\dagger \phi(x-\mu)-(8+2m) \phi(x) \right\}, \end{aligned}
\label{eq:DM}
\f} 

where \f$ U^R\f$ are the link variables in the representation \f$ R\f$.

Rescaling the fermion fields by \f$ \sqrt{\kappa}=\left(\frac{2}{8+2m}\right)^{1/2}\f$, we can write the fermionic action as:

\f{equation}{S_f = \sum_{x,y} \phi^\dagger(x) D_m(x,y) \phi(y),\f} 

where

\f{equation}{D_m(x,y) = \delta_{x,y} - \frac{\kappa}{2}
\left[(1-\gamma_\mu) U^R(x,\mu) \delta_{y,x+\mu} + 
(1+\gamma_\mu) U^R(x-\mu,\mu)^\dagger \delta_{y,x-\mu} \right],\f} 

and the Hermitian Dirac operator is obtained as

\f{equation}{Q_m = \gamma_5 D_m.  \label{eq:QM} \f} 

The fermionic determinant in the path integral can be represented by introducing complex pseudofermionic fields: 

\f{equation}{\left(\det D_m\right)^{N_f} = \int \mathcal D \phi \mathcal D \phi^\dagger e^{-\phi^\dagger Q_m^{-N_f} \phi} \equiv \int \mathcal D \phi \mathcal D \phi^\dagger e^{-S_\mathrm{pf}}.\f}

# Forces for the HMC molecular dynamics 

The HMC Hamiltonian is given by

\f{equation}{\mathcal{H}=\mathcal{H}_\pi+\mathcal{H}_G+\mathcal{H}_F \, ,\f} 

where

\f{equation}{\mathcal{H}_\pi = \frac{1}{2} \sum_{x,\mu} ( \pi(x,\mu) , \pi(x,\mu) ) = \frac{1}{2} T_f \sum_{a,x,\mu} \pi^a(x,\mu)^2 \, ,\f}

\f{equation}{\mathcal{H}_G = \beta \sum_{\mu<\nu} \left( 1- \frac{1}{N} \mathrm{Re\ tr\ } \mathcal{P}_{\mu\nu}\right) \, ,\f}

\f{equation}{\mathcal{H}_F = \phi^\dagger ( Q_m^2 - \beta )^{-l} \phi \, , \,\,\,\, l=\frac{N_f}{2}>0 \, , \label{eq:HF} \f}

and we have introduced for each link variable a conjugate momentum in
the algebra of the gauge group, defined as 

\f{equation}{\pi(x,\mu)=i \pi^a(x,\mu) T_f^a\, .\f}

In the expression of \f$ \mathcal{H}_F\f$ we omitted the sum over position, spin
and color indices and we have also introduced an arbitrary shift \f$ \beta\f$
for the matrix \f$ Q_m^2\f$, as this will be useful in the discussion for the
RHMC algorithm.

The equations of motion for the link variables are given by 

\f{equation}{\dot U(x,\mu) = \pi(x,\mu) U(x,\mu)\, .\f}

The notation \f$\dot{\square}\f$ indicates the derivative with respect to the molecular dynamics time. 

We obtain the equations of motion for the momenta from the requirement that the Hamiltonian \f$ \mathcal{H}\f$ is a conserved quantity

\f{equation}{
  0 = \dot{\mathcal{H}} = \dot{\mathcal{H}}_\pi + \dot{\mathcal{H}}_G + \dot{\mathcal{H}_F} \, .
  \label{eq:HCONS}
\f}

For the first two derivatives we have

\f{equation}{\dot{\mathcal{H}}_\pi = \sum_{x,\mu} ( \pi(x,\mu) , \dot\pi(x,\mu) ) = T_f \sum_{x,\mu} \sum_a \pi^a(x,\mu) \dot\pi^a(x,\mu) \,  \label{eq:HDOTPI} \f}

\f{equation}{
\begin{aligned}
\dot{\mathcal{H}}_{G} 
=& \sum_{x,\mu} -\frac{\beta}{N} \mathrm{Re\, tr} \left(\dot U(x,\mu) V^\dagger(x,\mu) \right) = \sum_{x,\mu} -\frac{\beta}{N} \mathrm{Re\, tr} \left(\pi(x,\mu) U(x,\mu) V^\dagger(x,\mu) \right) =\\ &\sum_{x,\mu} \sum_a -\frac{\beta}{N} \pi^a(x,\mu) \mathrm{Re\, tr} \left(i T^a_f U(x,\mu) V^\dagger(x,\mu) \right) \, , 
\end{aligned}
\label{eq:HDOTG}
\f}

where \f$ V(x,\mu)\f$ is the sum of the staples around the link \f$ U(x,\mu)\f$.

The computation of the fermionic force goes as follows. We only consider
the case \f$ l=1\f$ since this is the only case relevant both for the HMC
algorithm and the RHMC algorithm (see below). We have

\f{equation}{\begin{aligned} \dot{\mathcal{H}}_F = -\ \phi^\dagger (Q_m^2 - \beta)^{-1} \dot{(Q_m^2)} (Q_m^2 - \beta)^{-1} \phi \, . \end{aligned} \label{eq:FF1} \f}

Defining

\f{equation}{\eta = (Q_m^2 - \beta)^{-1} \phi \, ,  \label{eq:HMCETA} \f}
\f{equation}{\xi = Q_m \eta \, ,\f} 

and using the fact that the matrix \f$ (Q_m^2-\beta)\f$ is Hermitian, we can rewrite eq.\f$(\ref{eq:FF1})\f$ as

\f{equation}{\begin{aligned} \dot{\mathcal{H}}_F = - 2 \ \xi^\dagger \dot{(Q_m)} \eta \, .\end{aligned} \label{eq:FF2} \f}

Inserting the explicit form of \f$Q_m\f$, eq.\f$(\ref{eq:QM})\f$ and eq.\f$(\ref{eq:DM})\f$ into eq.\f$(\ref{eq:FF2})\f$ we obtain 

\f{equation}{\begin{aligned}\dot{\mathcal{H}}_F &= \mathrm{Re\ }\sum_{x,\mu} \xi(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \eta(x+\mu) + \xi(x+\mu)^\dagger \dot U^R(x,\mu)^\dagger \gamma_5 (1+\gamma_\mu) \eta(x) \\&= \mathrm{Re\ }\sum_{x,\mu} \xi(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \eta(x+\mu) + \eta(x)^\dagger \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \xi(x+\mu)\end{aligned}\,,\f}

where the sum over spin and color indices is intended and we made
explicit the fact that the whole expression is real. Further

\f{equation}{\dot U^R (x,\mu) = \pi^R(x,\mu) U^R(x,\mu) = i \pi^a(x,\mu) T^a_R U^R(x,\mu) \,. \label{eq:URDOT} \f}

Notice, that since we define \f$  T^a_R(x,\mu) = R_* T^a(x,\mu)\f$, the
\f$ \pi^a(x,\mu)\f$ in the above equation are the same as those appearing in
the expressions for \f$ \dot{\mathcal{H}}_{\pi,G}\f$. Using eq.\f$(\ref{eq:URDOT})\f$ in the
expression for \f$ \dot{\mathcal{H}}_{F}\f$ we find

\f{equation}{\begin{aligned}
\dot{\mathcal{H}}_F = \sum_{x,\mu} \sum_a  \pi^a(x,\mu) & \mathrm{Re\ Tr\ } \left[ iT^a_R U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \right. 
\left. \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right] \, .\end{aligned}\label{eq:HDOTF}\f}

Note that capitalized \f$ \mathrm{Tr}\f$ indicates the trace over both color and spin indices as opposed to the lower case \f$ \mathrm{tr}\f$, which is the trace over color only.

Inserting eq.\f$(\ref{eq:HDOTPI})\f$, eq.\f$(\ref{eq:HDOTG})\f$ into eq.\f$(\ref{eq:HCONS})\f$ we obtain the equations of motion for the momenta \f$ \pi^a(x,\mu)\f$ 

\f{equation}{\begin{aligned}
\dot\pi^a(x,\mu) &= \dot\pi^a_G(x,\mu) + \dot\pi^a_F(x,\mu) \, , 
\end{aligned}
\label{eq:PIDOT1} \f}

\f{equation}{\begin{aligned}
\dot\pi^a_G(x,\mu) &= \frac{\beta}{N} \frac{1}{T_f} \mathrm{Re\ tr\ } \left[ i T^a_f U(x,\mu) V^\dagger(x,\mu) \right] \, , 
\end{aligned}
\label{eq:PIDOT2} \f}

\f{equation}{\begin{aligned}
\dot\pi^a_F(x,\mu) &=-\frac{1}{T_f} \mathrm{Re\ Tr\ } \left[ iT^a_R U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right]\, . 
\end{aligned}
\label{eq:PIDOT3} \f}

For sake of convenience we introduce the following projectors \f$ P^a_R\f$ over the algebra in the representation \f$ R\f$

\f{equation}{P^a_R ( F ) = - \frac{1}{T_R} \mathrm{Re\ tr\ } \left[ i T^a_R F \right] \, ,\f}

which can be used to rewrite eqs eq.\f$(\ref{eq:PIDOT2})\f$ and eq.\f$(\ref{eq:PIDOT3})\f$ in a more compact form: 

\f{equation}{\begin{aligned}
\dot\pi^a_G(x,\mu) &= - \frac{\beta}{N} P^a_f \left( U(x,\mu) V^\dagger(x,\mu) \right) \, ,\\
\dot\pi^a_F(x,\mu) &= \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \left\{ \eta(x+\mu)\otimes\xi(x)^\dagger + \xi(x+\mu)\otimes\eta(x)^\dagger \right\} \right] \right)\, . \end{aligned}\label{eq:HFFORCE}\f}

### Checks of the MD force

The formulae derived in the previous section can be checked against two
known examples. The first, and almost trivial, check is obtained by
assuming that the representation \f$ R\f$ is again the fundamental
representation. The well-known expression for the MD force for the usual
HMC is then recovered.

The second case that has already been studied in the literature is the
case of fermions in the adjoint representation of the gauge group
SU(\f$ 2\f$) @cite Donini:1996nr. We agree with eq. (16) in
@cite Donini:1996nr, provided that we exchange the indices \f$ a\f$ and
\f$ b\f$ in that formula.

# HMC Algorithm

Given the action \f$ S(\phi)\f$ of a system of bosonic fields \f$ \phi\f$, our
goal is to generate a Markov process with fixed probability distribution
\f$ P_S(\phi) = Z^{-1} \exp[-S(\phi)]\f$. A sufficient condition to have such a
Markov process is that it is ergodic and it satifies detailed balance:

\f{equation}{P_S(\phi)P_M(\phi\rightarrow \phi') = P_S(\phi')P_M(\phi' \rightarrow \phi) \, .\f}

We define \f$ P_M(\phi \rightarrow \phi')\f$ with the following three-step
process:

1.  We expand the configuration space with additional fields, the
    momenta \f$ \pi\f$ randomly chosen with probability \f$ P_k(\pi)\f$ such
    that \f$ P_k(\pi)=P_k(-\pi)\f$ -- usually one takes
    \f$ P_k(\pi)\propto \exp[-\pi^2/2]\f$;

2.  In the extended configuration space \f$ (\phi, \pi)\f$, we generate a new
    configuration \f$ (\phi',\pi')\f$ with probability
    \f$ P_h((\phi,\pi)\rightarrow(\phi',\pi'))\f$ such that

    \f{equation}{P_h((\phi,\pi)\rightarrow(\phi',\pi')) = P_h((\phi',-\pi')\rightarrow(\phi,-\pi))\f}
    (reversibility condition)

3.  We accept the new configuration \f$ \phi'\f$ with probability

    \f{equation}{P_A((\phi,\pi)\rightarrow(\phi',\pi')) = \mathrm{min} \left\{ 1, \frac{P_S(\phi')P_k(\pi')}{P_S(\phi)P_k(\pi)} \right\} \, .\f}

    It is easy to see that the resulting probability

    \f{equation}{P_M(\phi\rightarrow\phi') = \int d\pi d\pi' P_k(\pi) P_h((\phi,\pi)\rightarrow(\phi',\pi')) P_A((\phi,\pi)\rightarrow(\phi',\pi')) \, ,\f}

    satisfies detailed balance. Care must be taken to ensure ergodicity.

As already stated, the distribution \f$ P_k(\pi)\f$ is generally taken to be
Gaussian (this should also guarantee ergodicity). The process \f$ P_h\f$ is
instead identified with the Hamiltonian flow of a yet unspecified
Hamiltonian \f$ H\f$ in the phase space \f$ (\phi,\pi)\f$ (giving to \f$ \pi\f$ the
meaning of "momenta"). The time reversal symmetry of classical dynamics
equation of motion guarantees the reversibility condition. The resulting
probability \f$ P_h\f$ is then a delta function (the process is completely
deterministic). Numerical integration to a given accuracy will result in
a broader distribution and care must be taken to guarantee the
reversibility condition in this case. Since we want a high acceptance
rate (low correlation among the configurations), we must carefully
choose the Hamiltonian \f$ H\f$. One simple way is to take \f$ P_k\f$ to be
Gaussian and define 

\f{equation}{H(\pi,\phi)=-\ln [P_k(\pi) P_S(\phi)] = \pi^2/2 + S(\phi)\f}

(omitting irrelevant constants). If \f$ H\f$ is exactly conserved by the process \f$ P_h\f$
then the acceptance probability is 1.

When fermionic degrees of freedom are present in the action \f$ S\f$, we can
first integrate them out, resulting in a non-local bosonic action and
then apply the above scheme. In practice, to deal with a non-local
action is not convenient from a numerical point a view and stochastic
estimates are used.

Consider a quadratic fermionic term in the action

\f{equation}{S(\bar\psi,\psi) = \bar\psi M \psi\f}

with a generic interaction matrix \f$ M(\phi)\f$ depending on the bosonic fields \f$ \phi\f$. The contribution of this term to the partition function is

\f{equation}{\int d\bar\psi d\psi \exp [ -S(\bar\psi,\psi)] = \mathrm{det}[M(\phi)]\, .\f}

Assuming that the matrix \f$ M(\phi)\f$ is positive definite, we can rewrite

\f{equation}{\mathrm{det}[M]=\int d\bar\eta d\eta \exp[ \bar\eta (M)^{-1} \eta ]\, ,\f}

where \f$ \bar\eta\f$,\f$ \eta\f$ are two new complex bosonic fields, called
pseudofermions. This term can be taken into account generating random
pseudofermions \f$ \bar\eta\f$, \f$ \eta\f$ with the desired probability
distribution and keeping then fixed during the above HMC configuration
generation for the remaining bosonic fields \f$ \phi\f$.

# RHMC formulation

The fermionic part of the HMC Hamiltonian, for \f$ N_f\f$ degenerate quarks
and \f$ N_{pf}\f$ pseudofermions, can be written as:

\f{equation}{\mathcal{H}_F = \sum_{k=1}^{N_{pf}} \phi_k^\dagger ( Q_m^2 )^{-l_k} \phi_k \,\, ;\,\, \sum_k l_k = \frac{N_f}{2}\, ,  \label{eq:HFN} \f}

and \f$ l_k>0\f$. For the sake of simplicity we will set all the \f$ l_k\f$ to be
equal: 

\f{equation}{\forall k,\,\, l_k = \frac{N_f}{2N_{pf}}\, .\f}

In the RHMC algorithm @cite Clark:2005sq rational approximations are used
whenever we need to take some fractional power of the positive definite
fermion matrix \f$ Q_m^2\f$.

In this implementation we use three different rational approximations. The first one is used to approximate eq.\f$(\ref{eq:HFN})\f$. Here, we need only one approximation because all \f$ l_k\f$ are equal yielding

\f{equation}{\mathcal{H}_F = \sum_{k=1}^{N_{pf}} \phi_k^\dagger r_{a}( Q_m^2 )\phi_k \, ,  \label{eq:HFRHMC} \f}
 
\f{equation}{( Q_m^2 )^{-\frac{N_f}{2N_{pf}}} \simeq r_{a}(Q_m^2) = \alpha_0^a + \sum_{n=1}^{d_{1}} \alpha_n^a ( Q^2_m - \beta_n^a )^{-1} \, .\f}

Using the formulae derived in the previous sections, it is easy to write the force corresponding to eq.\f$(\ref{eq:HFRHMC})\f$. In fact, eq.\f$(\ref{eq:HFRHMC})\f$ is nothing but a sum of terms of the form eq.\f$(\ref{eq:HFRHMC})\f$ once we put \f$ l=1\f$, \f$ \beta=\beta_n^a\f$. The RHMC force will be then given by a sum over \f$ n=1,\dots,d_1\f$ of terms given by
eq.\f$(\ref{eq:HFFORCE})\f$ multiplied by a factor \f$ \alpha_n^a\f$. Notice that since \f$ l=1\f$, to compute \f$ \eta\f$ as in eq.\f$(\ref{eq:HMCETA})\f$ a simple shifted inversion is required.

The second rational approximation is required in the heat bath update of
pseudofermions. In order to generate pseudofermions distributed as in
eq.\f$(\ref{eq:HFN})\f$, a simple two-step process is used. For each pseudofermion we first generate a gaussian distributed field \f$ \tilde\phi_k\f$

\f{equation}{P(\tilde\phi_k)\propto \exp [ -\tilde\phi_k^\dagger \tilde\phi_k ] \, ,\f}

and then set

\f{equation}{\phi_k = (Q_m^2)^{\frac{l_k}{2}} \tilde\phi_k \, ,\f}

making use of the fact that \f$ (Q_m^2)\f$ is Hermitian (notice the plus sign
in the exponent). The RHMC algorithm uses a rational approximation to
compute the above quantities. Again we need only one approximation since
all \f$ l_k\f$ are equal. 

\f{equation}{\begin{aligned}
( Q_m^2 )^{\frac{l_k}{2}} \simeq r_{b}(Q_m^2) &=& \alpha_0^b + \sum_{n=1}^{d_{2}} \alpha_n^b ( Q^2_m - \beta_n^b )^{-1} \, .\end{aligned}\f}

The third rational approximation is used in the code for the Metropolis
test. Starting from eq.\f$(\ref{eq:HFN})\f$ for each pseudofermion we can rewrite:

\f{equation}{\phi_k^\dagger ( Q_m^2 )^{-l_k}\phi_k = \left\| (Q_m^2)^{-\frac{l_k}{2}} \phi_k \right\|^2\, ,\f}

where we used the property that \f$ Q_m^2\f$ is Hermitian. The rational
approximation needed in this case is: 

\f{equation}{\begin{aligned}
( Q_m^2 )^{-\frac{l_k}{2}} \simeq r_{c}(Q_m^2) &=& \alpha_0^c + \sum_{n=1}^{d_{3}} \alpha_n^c ( Q^2_m - \beta_n^c )^{-1} \, .\end{aligned}\f}

Notice that if \f$ d_2=d_3\f$ the coefficients for the two approximations
\f$ r_b\f$ and \f$ r_c\f$ can each be obtained from the other.

In order to compute the coefficients \f$ \alpha_n\f$, \f$ \beta_n\f$ appearing in
the rational approximations the Remez algorithm is needed. In this
implementation we do not compute those coefficients "on the fly", but
rather we use a precomputation step to generate a table of coefficients
from which we pick up the right values when needed. The generation of
this table goes as follows:

First note that we need to compute rational approximations for a
function \f$ f(x)\f$ of the form \f$ f(x)=x^l\f$ and the approximation must be
accurate over the spectral range of the operator \f$ Q_m^2\f$. To simplify
the computation of the table we note that the following proposition
holds: if \f$ f(x)\f$ is a homogeneous function of degree \f$ l\f$ and \f$ r(x)\f$ is
an optimal (in the sense of relative error) rational approximation to
\f$ f(x)\f$ over the interval \f$ [\epsilon,\mathrm{h}]\f$ to a given accuracy
then \f$ r(kx)/k^l\f$ is an optimal rational approximation for the same
function and the same accuracy over the interval
\f$ [\epsilon/k,\mathrm{h}/k]\f$. Notice that the coefficients of the
"rescaled" rational approximation are easily obtained from that of the
original approximation. A simple corollary is that, given a homogeneuos
function \f$ f(x)\f$, we can divide the rational approximations with the same
accuracy in classes distinguished by the ratio \f$ \epsilon/\mathrm{h}\f$;
within each class the coefficients of the rational approximations are
easily related to each other, so that we only need to compute one
rational approximation in each class. This is what is done in our
implementation.

In detail: we generate a table containing the coefficients for the
rational approximations belonging in different classes distinguished by
the function \f$ f(x)\f$ which we want to approximate and the accuracy which
is required. We arbitrarily set \f$ \mathrm{h}\f$ to a fixed value equal to the
absolute upper bound on the spectrum of the matrix \f$ Q_m^2\f$. This choice
fixes the representative of each class, because the lower bound of the
approximation is now a function of \f$ \mathrm{h}\f$.

At run-time this table is used to generate optimal rational
approximations rescaling the precomputed coefficients to the desired
interval containing the spectrum of the matrix \f$ Q_m^2\f$. This interval is
obtained by computing the maximum and minimum eigenvalue of \f$ Q_m^2\f$ on
each configuration when needed. In our code we update this interval only
before the metropolis test, while we keep it fixed during the molecular dynamics.

# Even-Odd preconditioning

It is a very well know fact that the time spend for a simulation with
dynamical fermions is dominated by the time required for the inversions
of the Dirac operator. The convergence of such inversions can be
improved using appropriate preconditioning. The idea is to rewrite the
fermionic determinant as a determinant (or product of determinants) of
better conditioned matrix (matrices) than the original Dirac operator.
For the non-improved Wilson action this can be easily done using the
*even-odd* preconditioning. We start rewriting the Dirac operator \f$ D_m\f$
as a block matrix: 

\f{equation}{D_m = \begin{pmatrix}
4+m& D_{eo}\\
D_{oe} &4+m\\
\end{pmatrix}\,\,\, ,\f} 

where each block has a dimension half that of the original Dirac matrix. The diagonal blocks connecting sites with the same parity are proportional to the identity matrix, while off-diagonal blocks connect sites with opposite parity. We have (since \f$ D_m\f$ is \f$ \gamma_5\f$-hermitian):

\f{equation}{\gamma_5 D_{eo} \gamma_5 = D_{oe}^\dagger\,\, .\f} 

The determinant of the Dirac matrix \f$ D_m\f$ can be rewritten as:

\f{equation}{{\rm det\ } D_m = {\rm det\ } ( (4+m)^2 - D_{oe} D_{eo} ) = {\rm det\ } ( (4+m)^2 - D_{eo} D_{oe} ) \equiv {\rm det\ } D^{eo}_m\,\, ,\f}

using the well known formula for the determinant of a block matrix.
Since the determinant of \f$ D_m\f$ and of \f$ D_m^{eo}\f$ are the same the latter
can be used in numerical simulations. Note that the even-odd
preconditioned matrix only connects sites with the same parity thus it
have only half of the size of the original Dirac matrix and as \f$ D_m\f$ it
is \f$ \gamma_5\f$-Hermitian. We define as before the Hermitian matrix
\f$ Q_m^{eo}\equiv \gamma_5 D_m^{eo}\f$, which will be used in practice.

The formulation of the HMC algorithm does not change and the only
difference is that pseudofermionic fields are now only defined on half of
the lattice sites, conventionally the even sites in what follows. We now
give the explicit expression for the fermionic force for the
preconditioned system described by the Hamiltonian: 

\f{equation}{
\mathcal{H}_F = \phi_e^\dagger ( (Q^{eo}_m)^2 - \beta )^{-1} \phi_e \,\, ,\f}

where as before we are assuming \f$ N_f=2\f$ or a rational approximation of
the actual fractional power function, and where we made explicit that
\f$ \phi_e\f$ is only defined on even sites. Eq. eq.\f$(\ref{eq:FF2})\f$ is unchanged: 

\f{equation}{\dot{\mathcal{H}}_F = - 2 \ \xi_e^\dagger \dot{(Q^{eo}_m)} \eta_e \, , \label{eq:FFPRE} \f}

where as before we have defined: 

\f{equation}{\eta_e = ((Q^{eo}_m)^2 - \beta)^{-1} \phi_e \, ,\f}

\f{equation}{\xi_e = Q^{eo}_m \eta_e \, .\f} 

The explicit form of \f$ Q_m^{eo}\f$ must be used at this point. We have: 

\f{equation}{(\dot{Q}^{eo}_m) = -\gamma_5 (\dot{D}_{eo} D_{oe} + D_{eo}\dot{D}_{oe} )\,\, . \label{eq:QPREDOT} \f}

Defining 

\f{equation}{\sigma_o = D_{oe} \eta_e \, , \f}

\f{equation}{\rho_o = D_{oe} \xi_e \, ,\f} 

and inserting eq.\f$(\ref{eq:QPREDOT})\f$ into eq.\f$(\ref{eq:FFPRE})\f$ we find:

\f{equation}{\begin{aligned}
\dot{\mathcal{H}}_F = - \sum_{\mu,x\in \mathrm{even}} {\rm Tr}_{x,\mu}  \left[ \sigma_o(x+\mu)\otimes\xi_e(x)^\dagger + \rho_o(x+\mu)\otimes\eta_e(x)^\dagger \right] \\
- \sum_{\mu,x\in \mathrm{odd}} {\rm Tr}_{x,\mu}  \left[ \xi_e(x+\mu)\otimes\sigma_o(x)^\dagger + \eta_e(x+\mu)\otimes\rho_o(x)^\dagger \right] \end{aligned}\label{eq:FORPRE}\f}

employing the shorthand notation:

\f{equation}{{\rm Tr}_{x,\mu} \left[ \Phi \right] \equiv \mathrm{Re\ Tr\ } \left[ \dot U^R(x,\mu) \gamma_5 (1-\gamma_\mu) \Phi \right]\,\, .\f}

From eq.\f$(\ref{eq:FORPRE})\f$ it is clear that the fermionic force has a different expression on sites of different parities. Proceeding as before we arrive at the final expressions. For \f$ x\in \mathrm{even}\f$: 

\f{equation}{ \dot\pi^a_F(x,\mu) = - \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \left\{  \sigma_o(x+\mu)\otimes\xi_e(x)^\dagger + \rho_o(x+\mu)\otimes\eta_e(x)^\dagger \right\} \right] \right)\, ,\f}

while for \f$ x\in \mathrm{odd}\f$: 

\f{equation}{\begin{aligned}
\dot\pi^a_F(x,\mu) &= - \frac{T_R}{T_f} P^a_R \left( U^R(x,\mu) \mathrm{tr_{spin}} \left[ \gamma_5 (1-\gamma_\mu) \left\{\xi_e(x+\mu)\otimes\sigma_o(x)^\dagger + \eta_e(x+\mu)\otimes\rho_o(x)^\dagger   \right\} \right] \right)\, .\end{aligned}\f}

# Two-point functions

This is a summary of the formulae used for the mesonic two-point
functions. Let \f$ \Gamma\f$ and \f$ \Gamma^\prime\f$ be two generic matrices in the Clifford
algebra. Then we define the two-point function

\f{equation}{f_{\Gamma\Gamma^\prime}(t) = \sum_{\bf x}\langle \bar\psi({\bf x},t) \Gamma
\psi({\bf x},t) \bar\psi(0) \Gamma^\prime \psi(0) \rangle\, .\f} 

Performing the Wick contractions yields

\f{equation}{\langle \bar\psi({\bf x},t) \Gamma\psi({\bf x},t) \bar\psi(0) \Gamma^\prime \psi(0) \rangle = - \mathrm{tr} \left[ \Gamma S(x-y) \Gamma^\prime S(y-x) \right]
= - \mathrm{tr} \left[ \Gamma S(x-y) \Gamma^\prime \gamma_5 S^\dagger(x-y) \gamma_5 \right]\,.\f} 

In practice we invert the Hemitian Dirac operator \f$ \gamma_5 D\f$ by solving the equation:

\f{equation}{Q_{AB}(x-y) \eta^{\bar A,x_0}_B(y) = \delta_{A,\bar A} \delta_{x,x_0}\f}

where \f$ A=\{a,\alpha\}\f$ is a collective index for colour and spin, and
\f$ \bar A\f$, \f$ x_0\f$ are the position of the source for the inverter. Using the field \f$ \eta\f$ that we obtain from the inverter, the correlator
above becomes:

\f{equation}{\langle \ldots \rangle = - \tilde \Gamma_{AB} \eta^{C,y}_B(x)
\tilde \Gamma^\prime_{CD} \eta^{D,y}_A(x)^*\f} 

where \f$ \tilde \Gamma= \gamma_5 \Gamma\f$, and $\tilde \Gamma^\prime =
\gamma_5 \Gamma^\prime$.

# Hasenbusch acceleration

Let us summarize the Hasenbusch trick (for two flavours)

\f{equation}{\mathcal{H}_F = \phi^\dagger ( Q_m^2 )^{-1} \phi \,\f} 

where \f$ Q_m =\gamma_5 D_m\f$ is the hermitian Dirac operator. After integration over the pseudofermions it gives the determinant:

\f{equation}{\det{ Q_m ^2 } = \det{D_m^{\dagger} D_m}\f}

The Hasenbusch trick can be rewritten in the following form :

\f{equation}{\det{ Q_m ^2 }  = \det{W_- W_+} \det{\frac{ Q_m^2}{W_- W_+}}\f}

Where \f$ W_{\pm}\f$ can be chosen arbitrarily as long as the determinant is
well defined. We discuss in the next subsections various choices of
\f$ W_{\pm}\f$.

In any case the two term can be evaluated independently, and we have:

\f{equation}{\mathcal{H}_{F_1} =   \phi_1^\dagger ( W_- W_+ )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger Q_m^{-1} W_- W_+ Q_m^{-1} \phi_2\f}

This can be combined with even-odd preconditioning.

## Wilson Mass Shift

Assume 

\f{equation}{W_{+} = \left( D_m + \delta_m\right) ,\quad W_{-} =
W_{+}^\dagger=  \left( D^{\dagger}_m + \delta_m\right)\f}

Note that, as written in a comment in the code, $W_+ Q_m^{-1} = (a D +
b ) D^{-1} \gamma_5$.

Then 

\f{equation}{Q_m^{-1} \left( D^\dagger_m + \delta_m\right)  \left( D_m +
   \delta_m  \right) Q_m^{-1}   =  \left( \gamma_5 + \delta_m Q^{-1}
   \right)  \left( \gamma_5 +    \delta_m Q^{-1}\right)\f}

The force can then be computed : 

\f{equation}{\begin{aligned}
 \dot{\mathcal{H}_{F_2}} =&  - \delta_m \phi_2^\dagger \left[  \left( \gamma_5 + \delta_m Q^{-1}
   \right) \dot{Q^{-1}} + \dot{Q^{-1}} \left( \gamma_5 + \delta_m Q^{-1}
   \right)   \right] \phi_2 \\
=&  - \delta_m \phi_2^\dagger \left[  \left( \gamma_5 + \delta_m
    Q^{-1}\right) Q_m^{-1} \dot{Q} Q_m^{-1}  \right]\phi_2 + \rm{h.c}
  \end{aligned}\f}

Note that the equation as now the standard form of the forces for the
HMC algorithm provided that:

\f{equation}{X\equiv Q^{-1}\phi_2,\quad\textrm{and}\quad Y^{\dagger}=\phi_2^\dagger
(\gamma_5 + \delta_m Q^{-1}) Q_m^{-1}\f}

From which we deduce

\f{equation}{Y =   Q_m^{-1}(\gamma_5 + \delta_m Q^{-1})  \phi_2 =  D^{-1} ( \phi_2  +
\delta_m \gamma_5 X)\f}

Which matches one comment in the the force_hmc.c file.

### Even-Odd Preconditioning

Writing 

\f{equation}{D_m = \begin{pmatrix}
4 + m & D_{eo} \\
 D_{oe} & 4+m \\
\end{pmatrix}\f} 

The determinant in the 2 flavour case can be written as follows: 

\f{equation}{\det D_m^2 = \det Q^2 \propto  \det D^\dagger_{eo} D_{oe}\f}

\f{equation}{Q = \gamma_5 \begin{pmatrix} 
1 + 4m & M_{\rm{eo}} \\
M_{\rm{oe}} & 1 +4m\\
\end{pmatrix} \equiv \gamma_5 \begin{pmatrix} 
M_{\rm{ee}}  & M_{\rm{eo}} \\
M_{\rm{oe}} & M_{\rm{oo}}\\
\end{pmatrix}\f} 

Note that \f$ M_{\rm{ee}}^{-1}\f$ can be computed:

\f{equation}{M_{\rm{ee}}^{-1} = \frac{1}{1+4m}\f}

Now we can conveniently rewrite 

\f{equation}{Q_{\pm} = \begin{pmatrix} 
\gamma_5 M_{\rm{ee}}  & 0 \\ \gamma_5 M_{\rm{oe}} & 1\\
\end{pmatrix} \begin{pmatrix} 
1  &  \left(M_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \\0 & \gamma_5
  \left(M_{\rm{oo}} - \frac{1}{4+m}M_{\rm{oe}} M_{\rm{eo}} \right)\\
\end{pmatrix}\f}

From the last equation we deduce that:

\f{equation}{\det{Q} = \det{\gamma_5 M_{\rm{ee}}} \det{\gamma_5  \left(M_{\rm{oo}}
    -\frac{1}{4+m} M_{\rm{oe}}     M_{\rm{eo}} \right)} \propto
\det{\gamma_5  \left((4+m) M_{\rm{oo}}
    - M_{\rm{oe}}     M_{\rm{eo}} \right)}\f}

Note that the first determinant is a constant that could be computed.

In the following we will denote 

\f{equation}{\hat{Q}_{m,eo} \equiv \gamma_5
  \left((4+m)^2 - M_{\rm{oe}}  M_{\rm{eo}} \right)\f} 

where \f$ \hat{Q}_{m,eo}\f$ is defined on the odd sites of the lattice.

Now defining 

\f{equation}{W_{+} = D_{m+\delta m} ,\quad W_{-} =
W_{+}^\dagger=   D^{\dagger}_{m+\delta m}\f}

\f{equation}{\begin{aligned}
\det{ Q_m (W_-  W_+)^{-1}  Q_m} \propto \det{ Q_{m,eo} (\hat{D}_{m+\delta_m}
  \hat{D}_{m+\delta_m,eo}  )^{-1}  Q_{m,eo}}\end{aligned}\f}

We thus have 

\f{equation}{\begin{aligned}
\mathcal{H}_{F_1} = \phi_1^\dagger \left( \hat{D}_{m+\delta m,eo} \hat{D}_{m+\delta_m,eo} \right)^{-1}\phi_1\end{aligned}\f}

and 

\f{equation}{\begin{aligned}
\mathcal{H}_{F_2} = \phi_2^\dagger Q_{m,eo}^{-1} \hat{D}_{m+\delta_m,eo} \hat{D}_{m+\delta_m,eo} Q_{m,eo}^{-1}\phi_2 \end{aligned}\f}

Note that as in the non-even-odd case this can be rewritten as:

\f{equation}{\begin{aligned}
\mathcal{H}_{F_2} =\phi_2^{\dagger} (\gamma_5 + \delta_m (1 +
\delta_m (4+m) )Q_{m,eo}^{-1}   (\gamma_5 + \delta_m (1 +
\delta_m (4+m) )Q_{m,eo}^{-1}\phi_2 \end{aligned}\f}

## Twisted-Mass Shift

Assume

\f{equation}{W_{+} = \left( Q_m + i \mu_2  \right) ,\quad W_{-} = W_+^{\dagger}=\left( Q_m - i \mu_2 \right)\f}

Note that \f$ W_- W_+ = Q_m ^2 + \mu_2^2\f$ and that $W_\pm^{\dagger}=
W_{\mp}$.

Instead of dealing with $\det{ Q_m (W_- 
 W_+)^{-1} Q_m }$, we consider the slightly more general case where the
determinant to evaluate is 

\f{equation}{\begin{aligned}
\det{ (Q_m + i \mu_1) (W_- 
 W_+)^{-1}  (Q_m - i \mu_1)}\propto&\int D\phi_2 D\phi_2^\dagger e^{-\phi^\dagger_2 \Big(Q_+ (W_- 
 W_+)^{-1}  Q_- \Big)^{-1}\phi_2 }\big] \\
=& \int D\phi_2 D\phi_2^\dagger e^{-\phi_2Q_-^{-1} W_- W_+
  Q_+^{-1}  \phi_2 }\end{aligned}\f} 

The following formulae can then be used for the case of several Hasenbusch masses. The case of the determinant \f$ \det{ Q_m (W_- W_+)^{-1} Q_m }\f$ can be recovered by setting \f$ \mu_1=0\f$ in the following equations.

We have: 

\f{equation}{\begin{aligned}
(Q_m-&i\mu_1)^{-1} W_- W_+ (Q_m+i\mu_1)^{-1} = ( 1 - i(\mu_2 - \mu_1)(Q_m - i \mu_1)^{-1}) (1+ i(\mu_2 - \mu_1)(Q_m - i \mu_1)^{-1}) \\
=& 1+ i(\mu_2 - \mu_1) (Q_m+i\mu_1)^{-1} - i(\mu_2 - \mu_1)(Q_m - i \mu_1)^{-1} + (\mu_2- \mu_1)^2\big((Q_m + i \mu_1)(Q_m - i\mu_1)\big)^{-1} \\
=& 1+ (\mu_2- \mu_1)^2\big(Q_m^2 + \mu_1^2\big)^{-1}+ i(\mu_2 -
\mu_1) (Q_m^2 +\mu_1^2)^{-1}  (Q_m-i\mu_1) -i(\mu_2 - \mu_1) (Q_m^2 +  \mu_1^2)^{-1} (Q_m+i\mu_1) \\
=& 1 +(\mu_2- \mu_1) \big(Q_m^2 + \mu_1^2\big)^{-1} \big( (\mu_2- \mu_1)
+ 2 \mu_1 \big) \\
=& 1 +(\mu_2^2- \mu_1^2) \big(Q_m^2 + \mu_1^2\big)^{-1} \end{aligned}\f}

The force can then be computed: (global sign and factor \f$ i\f$ have to be
checked) 

\f{equation}{\begin{aligned}
\dot{\mathcal{H}_{F_2}} =& i(\mu_2-\mu_1) \phi_2^{\dagger} \Big[ ( 1 - i(\mu_2 - \mu_1) 
  (Q_m - i \mu_1)^{-1}) \dot{(Q_m+i\mu_1)^{-1}} 
 - \dot{(Q_m - i    \mu_1)^{-1}} (1+ i (\mu_2 - \mu_1) (Q_m+i \mu_1)^{-1})\Big] \phi_2 \\
=&  i(\mu_2-\mu_1) \phi_2^{\dagger} \left[  ( 1 - i(\mu_2 - \mu_1) 
  (Q_m - i\mu_1)^{-1}) ( Q_m+i\mu_1)^{-1} \dot{Q_m} (Q_m+i\mu_1)^{-1}
\right] \phi_2 +\rm{h.c}\end{aligned}\f}

with

\f{equation}{X\equiv (Q_m+i\mu_1)^{-1}\phi_2,\quad\textrm{and}\quad Y^{\dagger}=i\phi_2^\dagger
(1 -i (\mu_2-\mu_1) (Q-i\mu_1)^{-1}) (Q_m+i\mu_1)^{-1}\,.\f}

From this we deduce 

\f{equation}{\begin{aligned}
Y =&   -i (Q_m - i \mu_1)^{-1}(1 +  i(\mu_2-\mu_1)(Q+i\mu_1)^{-1})
\phi_2 \\ =&-i (Q_m - i \mu_1)^{-1}  ( \phi_2  + i(\mu_2-\mu_1)  X)  \,.\end{aligned}\f}

Note that in the particular case where \f$ \mu_1=0\f$,

\f{equation}{Q_m^{-1} W_- W_+ Q_m^{-1} = ( 1 - i\mu_2Q_m^{-1}) (1+ i\mu_2Q_m^{-1}))
= 1 + \mu_2^2 Q_m^{-2}\f}

Which leads to 

\f{equation}{\begin{aligned}
\dot{\mathcal{H}_{F_2}}  &=& \mu_2^2 \phi_2^{\dagger} \dot{Q_m^{-2}} \phi_2\end{aligned}\f}

Note also that the forces are explicitly proportional to \f$ \mu_2^2\f$.

### Even-Odd Preconditioning

Note that we have \f$ \widetilde{\mu} \equiv 2 \kappa \mu\f$.

\f{equation}{Q_{\pm}= \gamma_5 \begin{pmatrix} 
1 \pm i \widetilde{\mu} \gamma_5 & M_{\rm{eo}} \\
M_{\rm{oe}} & 1 \pm i \widetilde{\mu} \gamma_5\\
\end{pmatrix} \equiv \gamma_5 \begin{pmatrix} 
M^{\pm}_{\rm{ee}}  & M_{\rm{eo}} \\
M_{\rm{oe}} & M^{\pm}_{\rm{oo}}\\
\end{pmatrix}\f} 

\f$ M_{\rm{ee}}^{-1}\f$ can be computed as

\f{equation}{M_{\rm{ee}}^{-1} = ( 1 \pm i\widetilde{\mu} \gamma_5)^{-1} = \frac{1\mp
i \widetilde{\mu}\gamma_5 }{ 1 + \widetilde{\mu}^2}\,.\f}

Now we can conveniently rewrite 

\f{equation}{Q_{\pm} =  \begin{pmatrix} \gamma_5 M^{\pm}_{\rm{ee}}  & 0 \\ \gamma_5 M_{\rm{oe}} & 1\\
\end{pmatrix} \begin{pmatrix} 
1  &  \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \\0 & \gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)\\
\end{pmatrix}\f}

From the last equation we deduce that:

\f{equation}{\det{Q_{\pm}} = \det{\gamma_5 M^{\pm}_{\rm{ee}} } \det{\gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)}\f}

Note that the first determinant is a constant that could be computed. In the following we will denote 

\f{equation}{\hat{Q}_{\pm} \equiv \gamma_5
  \left(M^{\pm}_{\rm{oo}} - M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} \right)\f} 

where \f$ \hat{Q}_{\pm}\f$ is defined on the odd sites of the lattice.

We thus have

\f{equation}{\det{Q_+ Q_-} = \det{Q_+}\det{Q_{-}}\propto\det{\hat{Q}_+ \hat{Q}_-}\f}

and we thus get the following Hamiltonian:

\f{equation}{\mathcal{H}_{F_1} = \phi_1^{\dagger} \left(\hat{Q}_+
    \hat{Q}_-\right)^{-1} \phi_1\f}

The corresponding force then reads :

\f{equation}{\dot{\mathcal{H}_{F_1}} = - \phi_{0}^\dagger\left(    \hat{Q}_-^{-1}
  \hat{Q}_+^{-1}  \dot{\hat{Q}}_+  \hat{Q}_+^{-1}  + \hat{Q}_{-}^{-1}
  \dot{\hat{Q}}_{-}  \hat{Q}_-^{-1}   \hat{Q}_+^{-1}    \right)    \phi_0\f}

Now using that \f$ Q_{\pm} ^{\dagger} = Q_{\mp}\f$, the previous equation can
be written:

\f{equation}{\dot{\mathcal{H}_{F_1}} = - \left(Y_{\rm{o}}^{\dagger} \dot{\hat{Q}}_{+} X_{\rm{o}}  + \rm{h.c}\right)\f}

with

\f{equation}{X_{\rm{o}}=\hat{Q}_+^{-1}\phi_0,\quad Y_{\rm{o}}= \left(\hat{Q}_+\hat{Q}_-\right)^{-1} \phi_0,\f} 

where we have used that

\f{equation}{\hat{Q}_\pm^{\dagger} = \hat{Q}_\mp.\f} 

Furthermore we have

\f{equation}{\dot{\hat{Q}}_{\pm} =  \gamma_5 \left( -  \dot{M}_{\rm{oe}}
  \left(M^{\pm}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} -  M_{\rm{oe}}
    \left(M^{\pm}_{\rm{ee}}\right)^{-1} \dot{M}_{\rm{eo}}\right)\,.\f}

Now noting that 

\f{equation}{\dot{Q}_{\pm} =  \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}\f}

we have 

\f{equation}{\begin{aligned}
  Y ^{\dagger} \dot{Q} X  =& \begin{pmatrix} A^\dagger &
    B^\dagger \end{pmatrix} \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}  \begin{pmatrix} C \\  D \\\end{pmatrix} \\
=&  A^\dagger \gamma_5 \dot{M}_{\rm{oe}}  C + B^\dagger \gamma_5 \dot{M}_{\rm{eo}} D\end{aligned}\f}

Choosing \f$ A^\dagger= Y^\dagger_0\f$,
\f$ C=\left(M^{+}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} X_0\f$, \f$ B^\dagger=Y_0^\dagger \gamma_5 M_{\rm{oe}}\left(M^{+}_{\rm{ee}}\right)^{-1} \gamma_5\f$, and \f$ D= X_0\f$ allows to
write

\f{equation}{\dot{\mathcal{H}_{F_1}} =   Y^\dagger \dot{Q} X + \rm{h.c}\f} 

with

\f{equation}{ X=\begin{pmatrix}\left(M^{+}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} X_0
  \\ X_0 \end{pmatrix},\quad \rm{and} \quad Y=\begin{pmatrix}
   \left(M^{-}_{\rm{ee}}\right)^{-1} M_{\rm{eo}} Y_0
  \\ Y_0 \end{pmatrix}\,.\f} 

Here, we have used that \f$ \dot{Q}_+= \dot{Q}_{-}\f$ and

\f{equation}{ M_{\rm{eo}}^{\dagger}=\gamma_5 M_{\rm{oe}} \gamma_5 \,.\f}

#### Determinant Ratio

We use that 

\f{equation}{\det{ Q_+ (W_- 
 W_+)^{-1} Q_-} = \det{W_+^{-1} Q_+ Q_- W_-^{-1}} \propto\det{\hat{W}_+^{-1} \hat{Q}_+ \hat{Q}_- \hat{W}_-^{-1}}\,.\f}
  
We thus have to compute 

\f{equation}{\begin{aligned}
\dot{\mathcal{H}_{F_2}}  &=\, \phi_2^\dagger\big[ \delta \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+  +  \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \delta \hat{W}_+
+  \hat{W}_- \delta \hat{Q}_-^{-1} \hat{Q}_+^{-1} \hat{W}_+ +
\hat{W}_- \hat{Q}_-^{-1} \delta \hat{Q}_+^{-1} \hat{W}_+   \big] \phi_2 \\
&=\, \phi_2^\dagger\big[ \dot{\hat{W}}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+  +  \hat{W}_-
(\hat{Q}_+ \hat{Q}_-)^{-1} \delta \hat{W}_+
-  \hat{W}_- \hat{Q}_-^{-1} \dot{\hat{Q}}_- \hat{Q}_-^{-1} \hat{Q}_+^{-1} \hat{W}_+ -
\hat{W}_- \hat{Q}_-^{-1} \hat{Q}_+^{-1} \dot{\hat{Q}}_+ \hat{Q}_+^{-1} \hat{W}_+   \big] \phi_2 \\\end{aligned}\f}

Now we introduce

\f{equation}{X_W =  (\hat{Q}_+ \hat{Q}_-)^{-1} \hat{W}_+ \phi_2,  Y_W =
\hat{Q}_+^{-1} \hat{W}_+\phi_2 = \hat{Q}_- X_W\f}

such that 

\f{equation}{\dot{\mathcal{H}_{F_2}} = \phi_2^\dagger  \dot{\hat{W}}_- X_W +
X_W^\dagger  \delta \hat{W}_+ \phi_2 - Y_W^\dagger \dot{\hat{Q}}_- X_W - X_W^\dagger   \dot{\hat{Q}}_+  Y_W\f}

Now recalling that 

\f{equation}{\begin{aligned}
\dot{\hat{Q}}_{\pm} =&  -\gamma_5 \left(   \dot{M}_{\rm{oe}}
  \left(1\pm i\mu_1 \gamma_5\right)^{-1} M_{\rm{eo}}   + M_{\rm{oe}}
    \left( 1\pm i\mu_1 \gamma_5\right)^{-1} \dot{M}_{\rm{eo}}\right)  \\
\dot{\hat{W}}_{\pm} =&  -\gamma_5 \left(   \dot{M}_{\rm{oe}}
  \left(1\pm i\mu_2 \gamma_5\right)^{-1} M_{\rm{eo}} +  M_{\rm{oe}}
    \left( 1\pm i\mu_2 \gamma_5\right)^{-1} \dot{M}_{\rm{eo}}\right)\end{aligned}\f}

Now we can write the last expression in terms of \f$\dot{Q} \equiv \dot{Q}_{\pm}\f$. 

\f{equation}{\dot{\mathcal{H}_{F_2}} = Y_1^\dagger \dot{Q} X_1 + X_1^\dagger
\dot{Q} Y_1 - X_2^\dagger \dot{Q} Y_2  - Y_2^\dagger \dot{Q} X_2
= 2~ \rm{Re}\Big[ Y_1^\dagger \dot{Q} X_1 - Y_2^\dagger \dot{Q}
X_2 \Big]\,,\f} 

with 

\f{equation}{\begin{aligned}
Y_1 =& \begin{pmatrix}  (1+i\mu_1\gamma_5)^{-1} M_{\rm{eo}} Y_W \\
  Y_W \end{pmatrix},\quad Y_2  =  \begin{pmatrix}  (1+i\mu_2\gamma_5)^{-1} M_{\rm{eo}} \phi_2 \\
  \phi_2 \end{pmatrix},\\
X_{1,2} =&  \begin{pmatrix}  (1-i\mu_{1,2}\gamma_5)^{-1} M_{\rm{eo}} X_W \\
  X_W \end{pmatrix}, \qquad\dot{Q} \equiv  \dot{Q}_{\pm} =  \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}\end{aligned}\,.\f}

#### Twisted Wilson-Dirac Operator

Instead of applying the even-odd preconditioning to the twisted-mass operator
we can use the Wilson-Dirac even-odd operator and do a different splitting.

We define

\f{equation}{Q_{\rm{eo}} \equiv \gamma_5 \big( (4+m)^2 - D_{eo} D_{oe}\big)\,.\f}

Now split the determinant 

\f{equation}{\det{ Q_{\rm{eo}} ^2 }  = \det{W_- W_+} \det{\frac{  Q_{\rm{eo}} ^2
  }{W_- W_+}}\f}

and choose 

\f{equation}{W_{\pm} = Q_{\rm{eo}} \pm i \mu\,.\f}

The corresponding Hamiltonian reads

\f{equation}{\mathcal{H}_{F_1} = \phi_1^\dagger ( W_- W_+ )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger Q_{\rm{eo}}^{-1} W_- W_+ Q_{\rm{eo}}^{-1} \phi_2\f}

Since the operators are now very similar to the non even-odd case, we can reuse
some formulae. In particular, we can rewrite the Hamiltonian

\f{equation}{\mathcal{H}_{F_1} =  \phi_1^\dagger (Q_{\rm{eo}}  + \mu^2  )^{-1} \phi_1,\quad,
\mathcal{H}_{F_2} = \phi_2^\dagger \big( 1  + \mu^2 Q_{\rm{eo}}^{-1}\big)  \phi_2\,.\f}

From this we have the following forces: 

\f{equation}{\dot{\mathcal{H}}_{F_1} = \phi_1^\dagger  W_+^{-1} \delta W_-^{-1}\phi_1 + \rm{h.c}
= \phi_1^\dagger  W_+^{-1}  W_-^{-1} \dot{Q}_{\rm{eo}} W_-^{-1} \phi_1+ \rm{h.c}\,.\f}

Now we want to rewrite the last equation as a function of

\f{equation}{\dot{Q} = \gamma_5  \begin{pmatrix}  0 & \dot{D}_{\rm{eo}} \\
  \dot{D}_{\rm{oe}} & 0  \end{pmatrix}\f}

\f{equation}{X_{\rm{o}}=W_+^{-1}\phi_1,\quad Y_{\rm{o}}= \left(W_+
    W_-\right)^{-1} \phi_1,\f} 

where we have used that

\f{equation}{W_\pm^{\dagger} = W_\mp\,.\f} 

Furthermore we have

\f{equation}{\dot{Q}_{\rm{eo}} =  -\gamma_5 \left(  \dot{M}_{\rm{oe}}
   M_{\rm{eo}} +  M_{\rm{oe}} \dot{M}_{\rm{eo}}\right)\f}

noting that

\f{equation}{\begin{aligned}
  Y ^{\dagger} \dot{Q} X  =& \begin{pmatrix} A^\dagger &
    B^\dagger \end{pmatrix} \gamma_5 \begin{pmatrix} 
 0  & \dot{M}_{\rm{eo}} \\
\dot{M}_{\rm{oe}} & 0\\
\end{pmatrix}  \begin{pmatrix} C \\  D \end{pmatrix} \\
=&\,  A^\dagger \gamma_5 \dot{M}_{\rm{oe}}  C + B^\dagger \gamma_5 \dot{M}_{\rm{eo}} D\end{aligned}\f}

and chosing \f$ A^\dagger= Y^\dagger_0\f$, \f$ C =  M_{\rm{eo}} X_0\f$,
$B^\dagger=
 Y_0^\dagger \gamma_5 M_{\rm{oe}}\f$ , and \f$D= X_0$ allows us to write: 
 
\f{equation}{\dot{\mathcal{H}_{F_1}} = -Y^\dagger \dot{Q} X + \rm{h.c}\f} 
 
with 
 
\f{equation}{X=\begin{pmatrix} M_{\rm{eo}} X_0 \\ X_0 \end{pmatrix},\quad \rm{and} \quad Y=\begin{pmatrix} M_{\rm{eo}} Y_0 \\ Y_0 \end{pmatrix}\f} 

We have used that \f$ M_{\rm{eo}}^{\dagger}=\gamma_5 M_{\rm{oe}} \gamma_5\f$.

Similarly for the second Hamiltonian we get

\f{equation}{\begin{aligned} \dot{\mathcal{H}}_{F_2} = \mu_2\phi_2^\dagger  \dot{Q}_{\rm{eo}}^{-1}  \phi_2\,,\end{aligned}\f}

which is exactly the force that appears in case of a pure Wilson-Dirac even-odd preconditioned operator up to a multiplicative factor.

# Clover Term

The clover term can be written as

\f{equation}{ D_{sw} = -\frac{c_{sw}}{4}\sum_x\sum_{\mu,\nu}\sigma_{\mu\nu}F_{\mu\nu}(x),  \label{eq:clover} \f}

with the (unconventional) definition of \f$ \sigma_{\mu\nu}\f$ given by

\f{equation}{ \sigma_{\mu\nu} = \frac{1}{2}[\gamma_\mu,\gamma_\nu]. \f}

With the Euclidean definition of the gamma matrices, \f$ \sigma_{\mu\nu}\f$ satisfies

\f{equation}{ \sigma_{\mu\nu}^\dagger = \sigma_{\nu\mu} = -\sigma_{\mu\nu} = \sigma_{\mu\nu}^{-1}. \f}

For the Hermitian Dirac operator \f$ \gamma^5D\f$ we can make the following replacement without affecting any of the calculations presented here

\f{equation}{ \sigma_{\mu\nu} \to \bar{\sigma}_{\mu\nu} = \gamma_5\sigma_{\mu\nu}. \f}

The field strength tensor is defined as

\f{equation}{ F_{\mu\nu}(x) = \frac{1}{8}\left\{Q_{\mu\nu}(x) - Q_{\mu\nu}^\dagger(x)\right\} \f}

with

\f{equation}{ \begin{aligned}
 Q_{\mu\nu}(x)
 &= U_\mu(x)U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x) \\
 &+ U_\nu(x)U_\mu^\dagger(x-\hat{\mu}+\hat{\nu})U_\nu^\dagger(x-\hat{\mu})U_\mu(x-\hat{\mu}) \\
 &+ U_\mu^\dagger(x-\hat{\mu})U_\nu^\dagger(x-\hat{\mu}-\hat{\nu})U_\mu(x-\hat{\mu}-\hat{\nu})U_\nu(x-\hat{\nu}) \\
 &+ U_\nu^\dagger(x-\hat{\nu})U_\mu(x-\hat{\nu})U_\nu(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x)
\end{aligned} \f}

Because \f$ Q_{\mu\nu}^\dagger=Q_{\nu\mu}\f$ we have \f$ F_{\mu\nu}=-F_{\nu\mu}\f$.
For this reason we can change the sum over \f$ \mu,\nu\f$ in Eq. eq.\f$(\ref{eq:clover})\f$ to a sum over \f$ \mu<\nu\f$ and a factor of two.

\f{equation}{ D_{sw} = -\frac{c_{sw}}{2}\sum_x\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) \f}

The quantity \f$ \sigma_{\mu\nu}F_{\mu\nu}\f$ is Hermitian and block diagonal.
It can be written as

\f{equation}{ \begin{aligned}
 \sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu} =
 \begin{pmatrix}
 A & B & 0 & 0 \\
 B^\dagger & -A & 0 & 0 \\
 0  & 0 & C & D \\
 0 & 0 & D^\dagger & -C
 \end{pmatrix}
\end{aligned} \f}

with the definitions

\f{equation}{ \begin{aligned}
 A &= -iF_{03}+iF_{12} \\
 B &= -iF_{01}-F_{02}-F_{13}+iF_{23} \\
 C &= iF_{03}+iF_{12} \\
 D &= iF_{01}+F_{02}-F_{13}+iF_{23}
\end{aligned} \f}

## Pseudofermion Forces

For the forces we use the following short-hand notation for the derivative with respect to the link variables.

\f{equation}{ \dot{S} = \partial_{x,\mu}^a S \f}

To calculate the pseudofermion forces let us write down the action as

\f{equation}{ S = \phi^\dagger(H^{-2})\phi, \f}

where \f$ H=\gamma^5D\f$ is the Hermitian Dirac operator.
When differentiating the action we obtain

\f{equation}{ \dot{S} = -2\mathrm{Re}~\xi^\dagger\dot{H}\eta,  \label{eq:dotS} \f}

with the definitions

\f{equation}{ \begin{aligned}
 \eta &= H^{-2}\phi, \\
 \xi &= H\eta.
\end{aligned} \f}

### Forces

Here we will only consider the forces from the clover term and not the hopping term.
The clover part of the Dirac operator is given by

\f{equation}{ H_{sw} = -\frac{c_{sw}}{2}\sum_{\mu<\nu}\bar{\sigma}_{\mu\nu}F_{\mu\nu}(x)\,.  \label{eq:Hsw} \f}

When inserting Eq. eq.\f$(\ref{eq:dotS})\f$ we obtain

\f{equation}{ \dot{S} = c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{F}_{\mu\nu}\eta). \f}

From the definition of \f$ F_{\mu\nu}\f$ it follows that

\f{equation}{ \begin{aligned}
 \dot{S} &= \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \xi^\dagger\bar{\sigma}_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger\eta), \\
         &= \frac{1}{8}c_{sw}\sum_{\mu<\nu}\mathrm{Re}(\xi^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta + \eta^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\xi).
\end{aligned} \f}

This can in be written as

\f{equation}{ \dot{S} = \frac{1}{8}c_{sw}\sum_{\mu<\nu} \mathrm{Re}~\mathrm{tr}\left[\dot{Q}_{\mu\nu}\left\{\bar{\sigma}_{\mu\nu}\eta(x)\otimes\xi^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi(x)\otimes\eta^\dagger(x)\right\}\right] \f}

In a short hand notation we need to calculate

\f{equation}{ \dot{S} = \frac{1}{8}c_{sw}\mathrm{Re}~\mathrm{tr}[\dot{Q}_{\mu\nu}(x)X_{\mu\nu}(x)]  \label{eq:force} \f}

with

\f{equation}{ X_{\mu\nu}(x) = \bar{\sigma}_{\mu\nu}\eta(x)\otimes\xi^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi(x)\otimes\eta^\dagger(x) \f}

This matrix has the properties \f$ X_{\mu\nu}^\dagger=X_{\nu\mu}=-X_{\mu\nu}\f$.
The expression for \f$ \dot{Q}_{\mu\nu}(x)\f$ contains eight different terms (two from each of the four leafs).
The eight contributions to the force can be written as

\f{equation}{ \begin{aligned}
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
\end{aligned} \f}

where each term should be multiplied by \f$ c_{sw}/8\f$.
The calculation can be done efficiently by noticing that several products and terms appear in multiple places.
Introduce the intermediate variables

\f{equation}{ \begin{aligned}
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
\end{aligned} \f}

The total force can now be written as

\f{equation}{ F(x) = \frac{c_{sw}}{8}\dot{U}_\mu(x)\left\{W_8Z_0 + Z_1W_8 - W_5(W_0Z_2W_1+Z_3W_6) + (W_2Z_4W_3+W_7Z_5)W_4\right\} \f}

This brings us down to a total of 15 matrix multiplications and 6 additions.

### Logarithmic Forces

In the case of even-odd preconditioning the action of the small determinant \f$ D_{oo}\f$ can be written as

\f{equation}{ S_{sw} = -N_f\log\det D_{oo} = -N_f~\mathrm{tr}\log D_{oo} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\log D_{oo}(x) \f}

The derivative is given by

\f{equation}{ \dot{S} = -N_f\sum_{x~\mathrm{odd}}\mathrm{tr}\left[D_{oo}^{-1}(x)\dot{D}_{oo}(x)\right] \f}

with \f$ D_{oo}(x)\f$ given by

\f{equation}{ D_{oo}(x) = 4+m_0-\frac{c_{sw}}{2}\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x)\,.\f}

Both the determinant and the inverse of \f$ D_{oo}(x)\f$ can be calculated from an LDL decomposition.
If we insert the above definition we obtain

\f{equation}{ \dot{S} = \frac{N_fc_{sw}}{2}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{F}_{\mu\nu}(x)) \f}

\f{equation}{ \dot{S} = \frac{N_fc_{sw}}{2\cdot8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu} \mathrm{tr}(D_{oo}^{-1}(x)\sigma_{\mu\nu}\dot{Q}_{\mu\nu}(x) + D_{oo}^{-1}(x)\sigma_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger(x)) \f}

Since \f$ D_{oo}^{-1}\f$ is Hermitian we can write the result as two times the real part.
To simplify the result we define \f$ X_{\mu\nu}(x)=D_{oo}^{-1}(x)\sigma_{\mu\nu}\f$ such that

\f{equation}{ \dot{S} = \frac{N_fc_{sw}}{8}\sum_{x~\mathrm{odd}}\sum_{\mu<\nu}\mathrm{Re}~\mathrm{tr}[X_{\mu\nu}(x)\dot{Q}_{\mu\nu}(x)] \,.\f}

This is equivalent to eq.\f$(\ref{eq:force})\f$ except from the factor \f$ N_f\f$ and the definition of \f$ X_{\mu\nu}(x)\f$.
Notice that we still have the identity \f$ X_{\mu\nu}^\dagger=-X_{\mu\nu}\f$.
The sum over \f$ x\f$ can be extended to all sites by setting \f$ X_{\mu\nu}\f$ to zero on the even sites.
To calculate the inverse \f$ D_{oo}^{-1}\f$ we introduce the definitions:

\f{equation}{ D_{oo} = D_{oo}^\dagger = \begin{pmatrix}
 D_+ & 0 \\
 0 & D_-
 \end{pmatrix} \f}

\f{equation}{ D_{oo}^{-1} = \begin{pmatrix}
 D_{+}^{-1} & 0 \\
 0 & D_{-}^{-1}
 \end{pmatrix} \f}

\f{equation}{ D_{+}^{-1} =
 \begin{pmatrix}
  D_{11} & D_{12} \\
  D_{21} & D_{22} \\
 \end{pmatrix} \f}

\f{equation}{ D_{-}^{-1} =
 \begin{pmatrix}
  D_{33} & D_{34} \\
  D_{43} & D_{44} \\
 \end{pmatrix} \f}

Because of hermiticity we know that \f$ D_{12} = D_{21}^\dagger\f$ and \f$ D_{34} = D_{43}^\dagger\f$.
The six independent elements of \f$ X_{\mu\nu}\f$ can now be written as

\f{equation}{ \begin{aligned}
 X_{01} &= i(D_{34}+D_{43}) - i(D_{12}+D_{21}) \\
 X_{02} &= ~(D_{12}-D_{21}) + ~(D_{43}-D_{34}) \\
 X_{03} &= i(D_{22}-D_{11}) + i(D_{33}-D_{44}) \\
 X_{12} &= i(D_{11}-D_{22}) + i(D_{33}-D_{44}) \\
 X_{13} &= ~(D_{12}-D_{21}) + ~(D_{34}-D_{43}) \\
 X_{23} &= i(D_{12}+D_{21}) + i(D_{34}+D_{43})
\end{aligned} \f}


## Even-odd Preconditioning

### Method 1

We can write the determinant as

\f{equation}{ \det D = \det(D_{oo})\det(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe}) \f}

Use the notation

\f{equation}{ \begin{aligned}
 Q_{oo} &= \gamma_5D_{oo} \\
 Q &= \gamma_5(D_{ee}-D_{eo}D_{oo}^{-1}D_{oe})
\end{aligned} \f}

The action is

\f{equation}{ S = S_1 + S_2 = \phi_1^\dagger Q_{oo}^{-2}\phi_1 + \phi_2^\dagger Q^{-2}\phi_2 \f}

#### Forces for \f$ \phi_1\f$-term

The derivative is

\f{equation}{ \dot{S}_1 = -2\mathrm{Re}\left[\phi_1^\dagger(Q_{oo}^{-2}Q_{oo}\dot{Q}_{oo}Q_{oo}^{-2})\phi_1\right] \f}

and we can write it as

\f{equation}{ \dot{S}_1 = -2\mathrm{Re}\left[\xi^\dagger\dot{Q}_{oo}\eta\right] \f}

with

\f{equation}{ \begin{aligned}
 \eta &= Q_{oo}^{-2}\phi_1, \\
 \xi &= Q_{oo}\eta.
\end{aligned} \f}

#### Forces for \f$ \phi_2\f$-term

The derivative is

\f{equation}{ \dot{S}_2 = -2\mathrm{Re}\left[\phi_2^\dagger(Q^{-2}Q\dot{Q}Q^{-2})\phi_2\right] \f}

and we can write it as

\f{equation}{ \dot{S}_2 = -2\mathrm{Re}\left[\xi^\dagger\dot{Q}\eta\right] \f}

with

\f{equation}{ \begin{aligned}
 \eta &= Q^{-2}\phi_2, \\
 \xi &= Q\eta.
\end{aligned} \f}

The explicit expression for \f$ \xi^\dagger\dot{Q}\eta\f$ is given by

\f{equation}{ \xi^\dagger \dot{Q}\eta
 = \xi^\dagger\gamma_5\dot{D}_{ee}\eta
 - \xi^\dagger\gamma_5\dot{D}_{eo}D_{oo}^{-1}D_{oe}\eta
 + \xi^\dagger\gamma_5D_{eo}D_{oo}^{-1}\dot{D}_{oo}D_{oo}^{-1}D_{oe}\eta
 - \xi^\dagger\gamma_5D_{eo}D_{oo}^{-1}\dot{D}_{oe}\eta \f}

and it can be written as

\f{equation}{ \xi^\dagger \dot{Q}\eta
 = \xi^\dagger\gamma_5\dot{D}_{ee}\eta
 - \xi^\dagger\gamma_5\dot{D}_{eo}\eta_1
 + \xi_1^\dagger\gamma_5\dot{D}_{oo}\eta_1
 - \xi_1^\dagger\gamma_5\dot{D}_{oe}\eta \f}

with

\f{equation}{ \begin{aligned}
 \eta_1 &= D_{oo}^{-1}D_{oe}\eta \\
  \xi_1 &= D_{oo}^{-1}D_{oe}\xi
\end{aligned} \f}

### Method 2

The action of \f$ D_{oo}\f$ can also be expressed directly as the logarithm of the determinant.

\f{equation}{ S = -2\log\det D_{oo} + \phi^\dagger Q^{-2}\phi \f}

This is the approach implemented in the code.

## LDL factorization

With even-odd preconditioning we need to calculate the inverse \f$ D_{oo}^{-1}\f$ when applying the dirac operator and when calculating the forces.
Because this matrix is Hermitian and block diagonal it can be inverted locally with an exact solver.
The most practical solver is via an LDL decomposition.

\f{equation}{ A = LDL^\dagger \f}

Sum over \f$ j\f$

\f{equation}{ D_j = A_{jj} - \sum_{k=1}^{j-1}L_{jk}L_{jk}^*D_k \f}

Sum over \f$ i>j\f$.

\f{equation}{ L_{ij} = \frac{1}{D_j}\left(A_{ij}-\sum_{k=1}^{j-1}L_{ik}L_{jk}^*D_k\right) \f}

The determinant is given by

\f{equation}{ \det(A) = \prod_k D_k \f}

### LDL Decomposition

Calculates the LDL decomposition \f$ A=LDL^\dagger\f$ in-place.
After the decomposition, the lower triangular part of \f$ A\f$ is \f$ L\f$ and the diagonal is \f$ D\f$.

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

### Forward substitution

Calculates \f$ x=L^{-1}b\f$.

```
do i=0, N-1
  x_i = b_i
  do k=0, i-1
    x_i = x_i - A_ik * x_k
  enddo
enddo
```

### Backward substitution with diagonal

Calculates \f$ x=(L^\dagger)^{-1}D^{-1}x\f$.

```
do i=N-1, 0
  x_i = x_i/A_ii
  do k=i+1, N-1
    x_i = x_i - conj(A_ki) * x_k
  enddo
enddo
```

### Full inversion

This algorithm calculates the inverse \f$ B=A^{-1}\f$ from the LDL decomposition.
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


# Exponential Clover Term

The exponential version of the clover term (including mass term) can be
written as

\f{equation}{D_{sw} = \sum_x  (4+m_0) \exp ( A(x) ), \text{ with } A(x) =  -\frac{c_{sw}}{4(4+m_0)}\sum_{\mu,\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) \f} (eq:clover_exp)

where \f$ \sigma_{\mu\nu}\f$ is again defined by

\f{equation}{\sigma_{\mu\nu} = \frac{1}{2}[\gamma_\mu,\gamma_\nu].\f}

As for the clover term above, we can simplify the sum over \f$ \mu\nu\f$ to a sum over \f$ \mu < \nu\f$ and introduce a factor of two. We define

\f{equation}{A(x) = -\frac{c_{sw}}{2(4+m_0)}\sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu}(x) \f} (eq:clover2)

The quantity \f$ \sigma_{\mu\nu}F_{\mu\nu}\f$ is Hermitian and block
diagonal. It can be written as 

\f{equation}{\begin{aligned} A(x)=-\frac{c_{sw}}{2(4+m_0)} \sum_{\mu<\nu}\sigma_{\mu\nu}F_{\mu\nu} = 
 \begin{pmatrix}
 a & b & 0 & 0 \\
 b^\dagger & -a & 0 & 0 \\
 0  & 0 & c & d \\
 0 & 0 & d^\dagger & -c
 \end{pmatrix} \equiv \begin{pmatrix}
 A^+ & 0 \\ 
 0 & A^- 
 \end{pmatrix}, \end{aligned}\f} (eq:blocks)

where \f$ A^\pm\f$ are \f$ 2 \times 2\f$ matrices in spin space and \f$ a,b,c,d\f$ are \f$ N_F \times N_F\f$.

This formulation of \f$ A(x)\f$ as a block matrix will be useful for the
exponentiation.

## Evaluation of the operator

The evaluation of the exponential of \f$ A(x)\f$ can be split as:

\f{equation}{\exp A(x)  = \begin{pmatrix}
\exp A^+ & 0 \\
0& \exp A^-
\end{pmatrix}\f} 

and so, the problem is reduced to the exponential of two
\f$ (2 N_F) \times  (2 N_F)\f$ matrices. The evaluation can be performed in
two ways.

1.  Using the Taylor expansion:

    \f{equation}{\exp(A^\pm) = \sum_{k=0}^N  \frac{1}{k!} (A^\pm)^k.\f}

2.  Using the Horner scheme:

\f{equation}{\exp(A^\pm) = \sum_{k=0}^{\dim A^\pm -1} b_k(A^\pm)  (A^\pm)^k,  \label{eq:exphorner} \f}

where \f$ b_k\f$ are computed recursively as follows. We start with

\f{equation}{\begin{aligned} q_{N,0} = 1/N!, q_{N,1 \cdots (\dim A^\pm)-1} = 0.\end{aligned}\f}

Then, the recursion proceeds: 

\f{equation}{q_{n,0} = - p_0 q_{n+1, \dim A^\pm-1}  + 1/n!,\f}
\f{equation}{q_{n,i} = - p_i q_{n+1, \dim A^\pm-1}  + q_{n+1,i-1},\f}
\f{equation}{\text{ with } i=1 \cdots (\dim A^\pm) -1, \label{eq:horner} \f} 

where \f$ p_i\f$ represent the coefficients of the characteristic polynomial of the matrix \f$ A^\pm\f$

\f{equation}{P(A^\pm) = \sum_{n=0}^{\dim A\pm} p_n (A^\pm)^n.\f} 

For instance, the characteristic polynomial of a \f$ 4 \times 4\f$ traceless matrix has the following coefficients:

\f{equation}{p_0=\frac{1}{8 } \left(\mathrm{tr}A^2\right)^2 - \frac{1}{4} \mathrm{tr}A^4 ,\ p_1 = -\frac{1}{3}\mathrm{tr}A^3, \ p_2= -\frac{1}{2}\mathrm{tr}A^2, \ p_3=0, \ p_3=1.\f}

Finally, the coefficients of eq.\f$(\ref{eq:exphorner})\f$ are \f$ b_k  =q_{0,k}\f$.

The Horner scheme method is currently implemented only for \f$ SU(2)\f$ and
\f$ SU(3)\f$ with fundamental fermions.

## Pseudofermion Forces

For the forces we use the following shorthand notation for the
derivative with respect to the link variables.

\f{equation}{\dot{S} = \partial_{x,\mu}^a S\f} 

To calculate the pseudofermion
forces let us write down the action as 

\f{equation}{S = \phi^\dagger(H^{-2})\phi,\f}

where \f$ H=\gamma^5D\f$ is the Hermitian Dirac operator. When differentiating the action we obtain

\f{equation}{\dot{S} = -2\mathrm{Re}~\xi^\dagger\dot{H}\eta, \f} (eq:dotS)

with the definitions 

\f{equation}{\begin{aligned}
 \eta &= H^{-2}\phi, \\
 \xi &= H\eta.\end{aligned}\f}
 
 #### Forces

Here we will only consider the forces from the clover term. For the
exponential version of the clover term, the implementation is very
similar to the traditional clover term.

The clover part of the Dirac operator is given by

\f{equation}{H_{sw} = (4+m_0) \gamma_5 \exp \left(- \frac{c_{sw}}{2(4+m_0)}  \sum_{\mu<\nu}{\sigma}_{\mu\nu}F_{\mu\nu}(x)  \right) = (4+m) \gamma_5 \exp A(x).\label{eq:expHsw}\f} 

An optimized way of calculating the derivative is provided by the double
Horner scheme. The basic idea is that the derivative of a matrix can be
expressed as:

\f{equation}{d e^A = \sum_{k=0}^{\dim A -1} \sum_{l=0}^{\dim A-1} C_{kl}(A) A^l dA A^k, \label{eq:dexphorner}\f} 

where the \f$ C_{kl}\f$ coefficients depend on the matrix \f$ A\f$, similarly to
the ones eq.\f$(\ref{eq:exphorner})\f$. They are calculated performing first the
iteration in eq.\f$(\ref{horner})\f$, and then repeating the iteration process on the result of the first iteration. For compactness, we shall omit the limits of the sum henceforth. When inserting eq.\f$(\ref{eq:expHsw})\f$ in eq.\f$(\ref{eq:dotS})\f$, and using eq.\f$(\ref{eq:dexphorner})\f$ we obtain 

\f{equation}{\dot{S} = c_{sw} \sum_{k}  \sum_{\mu<\nu}\mathrm{Re}(  \xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{F}_{\mu\nu}\eta_k),\f}

with 

\f{equation}{\xi_k = \sum_l \begin{pmatrix}
C^+_{kl} (A^+)^l \xi^+ \\
C^-_{kl} (A^-)^l \xi^-
\end{pmatrix}, \text{ and } 
\eta_k =  \begin{pmatrix}
 (A^+)^k \eta^+ \\
(A^-)^k \eta^-
\end{pmatrix}.\f}

From the definition of \f$ F_{\mu\nu}\f$ it follows that 

\f{equation}{\begin{aligned}
 \dot{S} &=  \frac{1}{8}c_{sw}\sum_k \sum_{\mu<\nu}\mathrm{Re}(\xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta_k + \xi_k^\dagger\bar{\sigma}_{\mu\nu}^\dagger\dot{Q}_{\mu\nu}^\dagger\eta_k), \\
 &=
 \frac{1}{8}c_{sw}\sum_k \sum_{\mu<\nu}\mathrm{Re}(\xi_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\eta_k + \eta_k^\dagger\bar{\sigma}_{\mu\nu}\dot{Q}_{\mu\nu}\xi_k).\end{aligned}\f}

This can in be written as

\f{equation}{\dot{S} = \frac{1}{8}c_{sw}\sum_{\mu<\nu} \mathrm{Re}~\mathrm{tr}\left[\dot{Q}_{\mu\nu} \sum_k\left\{\bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x)\right\}\right]\f}

As for the clover term above we need to calculate now

\f{equation}{\dot{S} = \frac{1}{8}c_{sw}\mathrm{Re}~\mathrm{tr}[\dot{Q}_{\mu\nu}(x)X_{\mu\nu}(x)]  \label{eq:force2} \f}

now with

\f{equation}{X_{\mu\nu}(x) = \sum_k \bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x).\f}

The total force can now be expressed as in the clover term above. 

## Even-odd Preconditioning

Even-odd preconditioning is particularly simple for the exponential
case, since the force coming from the little determinant vanished. This
can be seen because of the fact that:

\f{equation}{\det D_{oo}  = \exp(\log \det D_{oo} ) =\exp( \mathrm{tr}\log D_{oo}) =  1,\f}

and so it is a constant term in the action that does not contribute to
the force.

## Implementation using Taylor Series 

In the current version of the code, the Horner scheme is only implemented
for \f$ SU(2)\f$ and \f$ SU(3)\f$ with fundamental fermions. For other theories, a
less efficient, but more flexible, alternative is used. For this,
we use the Taylor series

\f{equation}{dA = \sum_{k=0}^N \sum_{l=0}^{N-k} \frac{1}{(k+l+1)!} A^{k} dA A^{l},\f}

with \f$ N\f$ sufficiently large. The implementation changes only in the definition of \f$ X_{\mu\nu}\f$:

\f{equation}{X_{\mu\nu}(x) = \sum_{k=0}^N \bar{\sigma}_{\mu\nu}\eta_k(x)\otimes\xi_k^\dagger(x) + \bar{\sigma}_{\mu\nu}\xi_k(x)\otimes\eta_k^\dagger(x),\f}

where now 

\f{equation}{\xi_k = \sum_l \frac{1}{(k+l+1)!} \begin{pmatrix}
(A^+)^l \xi^+ \\
 (A^-)^l \xi^-
\end{pmatrix}, \text{ and } 
\eta_k =  \begin{pmatrix}
 (A^+)^k \eta^+ \\
(A^-)^k \eta^-
\end{pmatrix}.\f}

# Stout smearing 

The implementation follows @cite Morningstar:2003gk closely. We define the
smeared links as 

\f{equation}{U'_\mu(x) = e^{Q_\mu(x)}U_\mu(x)\f}

where \f$ \Sigma_\mu(x)\f$ is an element of the Lie algebra, defined via the
projection 

\f{equation}{Q_\mu(x) = \mathcal{P}\{\Omega_\mu(x)\}.\f}

The projection operator is not unique, but the most common choice is

\f{equation}{\mathcal{P}(\Omega) = \frac{1}{2}(\Omega-\Omega^\dagger) - \frac{1}{2N}\mathrm{tr}(\Omega-\Omega^\dagger).\f}

However, in a setup with mixed representations, it is convenient to use
the following two-step procedure for the projection. This allows us to
project matrices from different representations onto the generators of
the fundamental representation. 

\f{equation}{\begin{aligned} A_\mu^a(x) &= -\frac{1}{T_f}\mathrm{tr}[iT^a_R\Omega_\mu(x)] \\ Q_\mu(x) &= iT^a_F A_\mu^a(x)\end{aligned}\f}

The matrix \f$ \Omega_\mu(x)\f$ is defined as

\f{equation}{\Omega_\mu(x) = U_\mu(x)V_\mu(x)\f} 

\f{equation}{V_\mu(x) = \sum_{\nu\neq\mu}
 \rho_{\mu\nu}\left(
 U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x) +
 U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})
 \right)\f}

For the force calculation we use the chain rule.

\f{equation}{\frac{dS}{dU} = \frac{dS}{dU'}\frac{dU'}{dU}\f}

The first derivative on the right-hand side is the usual force
\f$ \Sigma_\mu'(x)\f$ evaluated using the smeared links. The second term is
the derivative of the smeared links with respect to the fundamental
links. This can be written in the following way, because the derivative
of the action is surrounded by a trace.

\f{equation}{\frac{dS}{dU} = e^Q\Sigma' + \frac{d\Omega}{dU}\mathcal{P}(X)\f}

When using a Taylor expansion to define the exponential function, we can
use the following definition of \f$ X\f$.

\f{equation}{X = \sum_{n=0}^\infty\sum_{k=0}^n Q^kU\Sigma' Q^{n-k}\f}

The derivative of the \f$ \Omega\f$ matrix is the last missing piece. Define
\f$ \Lambda=\mathcal{P}(X)\f$ and consider

\f{equation}{\frac{d}{d U_\mu(x)} U_\mu(x)V_\mu(x)\Lambda_\mu(x)\f}

Here we have a sum over \f$ \nu\neq\mu\f$. There are eight contributions to
the above derivative. 

\f{equation}{\rho_{\mu\nu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\Lambda_\mu(x)\f} 

\f{equation}{\rho_{\nu\mu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\Lambda_\nu(x-\hat{\nu})U_\nu(x-\hat{\nu})\f}

\f{equation}{\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\Lambda_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})\f}

\f{equation}{\rho_{\nu\mu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\Lambda_\nu^\dagger(x)\f}


\f{equation}{\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})\Lambda_\mu(x)\f}

\f{equation}{\rho_{\nu\mu}U_\nu^\dagger(x-\hat{\nu}+\hat{\mu})\Lambda_\nu^\dagger(x-\hat{\nu}+\hat{\mu})U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu})\f}

\f{equation}{\rho_{\mu\nu}U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})\Lambda_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\f}

\f{equation}{\rho_{\nu\mu}\Lambda_\nu(x+\hat{\mu})U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})U_\nu^\dagger(x)\f}

This can be simplified because several products appear more than once
and we can use \f$ \Lambda^\dagger=-\Lambda\f$ to remove some of the
Hermitian conjugates. In the following we also assume that
\f$ \rho_{\mu\nu}=\rho_{\nu\mu}\f$. 

\f{equation}{\rho_{\mu\nu}W_2U_\nu^\dagger(x)\{\Lambda_\mu(x)-\Lambda_\nu(x)\}\f}

\f{equation}{\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})U_\mu^\dagger(x-\hat{\nu})\{\Lambda_\nu(x-\hat{\nu})-\Lambda_\mu(x-\hat{\nu})\}U_\nu(x-\hat{\nu})\f}

\f{equation}{\rho_{\mu\nu}U_\nu^\dagger(x+\hat{\mu}-\hat{\nu})\{W_1\Lambda_\mu(x)-\Lambda_\nu(x+\hat{\mu}-\hat{\nu})W_1\}\f}

\f{equation}{\rho_{\mu\nu}\{\Lambda_\nu(x+\hat{\mu})W_2-W_2\Lambda_\mu(x+\hat{\nu})\}U_\nu^\dagger(x)\f}

Here 

\f{equation}{\begin{aligned}
 W_1 &= U_\mu^\dagger(x-\hat{\nu})U_\nu(x-\hat{\nu}) \\
 W_2 &= U_\nu(x+\hat{\mu})U_\mu^\dagger(x+\hat{\nu})\end{aligned}\f}

This brings us down to 13 multiplications.