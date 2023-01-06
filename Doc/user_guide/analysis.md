@page analysis Analysis

## Contractions

### Connected Two-Point Correlation Functions

First we choose interpolating operators with the quantum numbers of the
meson we would like to study:

\f{equation}{O_M = \bar{\psi}^{(1)}(x) \Gamma \psi^{(2)}(x)\f}

\f{equation}{\bar{O}_M = \bar{\psi}^{(2)}(x) \bar{ \Gamma } \psi^{(1)}(x)\f}

where  \f$ \bar{ \Gamma } = \gamma_0 \Gamma^\dagger \gamma_0 \f$. The gamma matrix is
chosen to have the same \f$ J^{PC} \f$ quantum numbers as the meson we want to
create. We then calculate the expectation value to create a meson at  \f$ y \f$ 
and destroy it again at  \f$ x \f$ : 

\f{equation}{
    \begin{aligned}
\langle O_M(x) \bar{O}_M'(y) \rangle =& \langle \bar{\psi}^{(1)}(x) \Gamma \psi^{(2)}(x) \bar{\psi}^{(2)}(y) \bar{ \Gamma }' \psi^{(1)}(y) \rangle 
\\
=& \Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \bar{\psi}^{(1)}_\alpha(x)  \psi^{(2)}_\beta(y) \bar{\psi}^{(2)}_\gamma(y)  \psi^{(1)}_\delta(x) \rangle \\
=& -\Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \psi^{(2)}_\beta(x) \bar{\psi}^{(2)}_\gamma(y)  \psi^{(1)}_\delta(y) 
\bar{\psi}^{(1)}_\alpha (x) \rangle
\end{aligned}
\label{eq:2ptcorr}
\f}

where the sign in the last line comes from exchanging the anticommuting \f$ \bar{\psi}_\alpha \f$ three times. Wick contracting where \f$ \langle \psi_\alpha(x) \bar{\psi}_\beta(y) \rangle = S_{\alpha \beta}(x,y) \f$ gives 

\f{equation}{\begin{aligned}
\langle O_M(x) \bar{O}_M'(y) \rangle =& -\Gamma_{\alpha \beta} S^{(2)}_{\beta \gamma} (x,y) \bar{ \Gamma }'_{\gamma \delta} S^{(1)}_{\delta \alpha} (y,x) \\
=& -\text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]\end{aligned}
\f}

To get correlation functions we Fourier transform, ie. sum over all
source and sink seperations while projecting onto the momentum we want
to give to the particle: 

\f{equation}{\begin{aligned}
C(t - \tau, \vec{p}) =& \sum_{\vec{x} \vec{y}} e^{-i \vec{p}(\vec{x} - \vec{y}) } \langle O_M(\vec{x}, t) O_M'(\vec{y}, \tau)\rangle \\
=& -\sum_{\vec{x} \vec{y}} e^{-i \vec{p}(\vec{x} - \vec{y}) } \text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]\end{aligned}\f}

here \f$ x = (\vec{x}, t) \f$ and \f$ y = (\vec{y}, \tau) \f$ and the zero momentum
correlator is:

\f{equation}{C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]\f}

We now use \f$ \gamma_5 \f$  Hermiticity: \f$ \gamma_5 S^\dagger(x,y) \gamma_5 = S(y,x) \f$, 

\f{equation}{C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \gamma_5 \Gamma S^{(2)} (x,y) \bar{ \Gamma }' \gamma_5 S^{\dagger (1)} (x,y) \right]\label{eqn:corr}\f}

### Point Sources

Using a delta function source and solving the Dirac equation gives a
point propagator,

\f{equation}{D_{\alpha a, \beta b} (x,y) S_{\beta b, \gamma c} (y, z) = \delta(x,z) \delta_{ac} \delta_{\alpha \gamma}\f}

usually \f$ z = (\vec{0}, 0) \f$  so we get  \f$ S(y,0) = \gamma_5 S^{\dagger} (0,y) \gamma_5 \f$. Then we use these to
calculate correlation functions, 

\f{equation}{C(t, 0) = -\sum_{\vec{x}} \text{Tr} e^{-i \vec{p} \vec{x}} \left[ \gamma_5 \Gamma S (x,0) \bar{ \Gamma }' \gamma_5 S (x,0) \right]\label{eqn:pointcorr}\f}

Translational invariance in the limit of infinitely many gauge
configurations implies  \f$ S(x,y) = S(|x - y|) \f$, so the sum over \f$ \vec{y} \f$ 
in equation eq.\f$(\ref{eqn:corr})\f$ just gives \f$ V \f$  times equation eq.\f$(\ref{eqn:pointcorr})\f$. We place the source at the time origin so  \f$ \tau = 0 \f$.

### One-end Trick

For this method it helps to write all the indices out,

\f{equation}{C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,bc} (x,y) (\bar{ \Gamma }' \gamma_5)_{\gamma \delta} S^{\dagger (1)}_{\delta \alpha, cb} (x,y)\f}

Greek indices \f$ \alpha, \beta, \gamma, \delta, ... \f$ are spinor indices
and Latin indices \f$ a, b, c, d... \f$ are colour indices. The one-end trick
involves inserting a delta function in colour, spin and space.

\f{equation}{C(t - \tau, 0) = -\sum_{\vec{x} \vec{y} \vec{z}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) \delta_{\gamma \lambda} \delta_{cd} \delta(\vec{y}, \vec{z}) (\bar{ \Gamma }' \gamma_5)_{\lambda \delta} S^{\dagger (1)}_{\delta \alpha, d b} (\vec{x}, t;\vec{z}, \tau)\f}

The delta function is aproximated with a \f$ Z(2) \times Z(2) \f$ noise source on timeslice \f$ \tau \f$ 

\f{equation}{\delta_{\gamma \lambda} \delta_{cd} \delta(\vec{y}, \vec{z}) \approx \frac{1}{K} \sum_{k = 0}^{K} | \eta^{(k)}_{\gamma c }(\vec{y})\rangle \langle \eta^{(k)}_{\lambda d }(\vec{z}) |\f}

which is exact in the limit \f$ K \rightarrow \infty \f$. 

\f{equation}{
C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\vec{x} \vec{y} \vec{z}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle \langle \eta^{(k)}_{\lambda d }(\vec{z}) | (\bar{ \Gamma }' \gamma_5)_{\lambda \delta} S^{\dagger (1)}_{\delta \alpha, d b} (\vec{x}, t;\vec{z}, \tau)\label{eq:oet}\f}

Defining

\f{equation}{\phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) = \sum_{\vec{y}} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle\f}

and

\f{equation}{\phi^{\Gamma (k)}_{\alpha, b}(\vec{x}, t; \tau) = \sum_{\vec{z}} S^{(1)}_{\alpha \delta , b d} (\vec{x}, t;\vec{z}, \tau) (\bar{ \Gamma }' \gamma_5)^{\dagger}_{\delta \lambda} | \eta^{(k)}_{\lambda d }(\vec{z}) \rangle\f}

the correlator can be evaluated as,

\f{equation}{C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) \phi^{\Gamma \dagger (k)}_{\alpha, b}(\vec{x}, t; \tau)\f}

#### Implementation in HiRep

In HiRep the code ```Spectrum/mk_mesons_with_z2semwall.c``` does two solves to
calculate  \f$ S^{(1)} | \eta \rangle \f$ and \f$ S^{(2)} (\bar{ \Gamma }' \gamma_5)^{\dagger} | \eta \rangle \f$. HiRep has

\f{equation}{\rho = \rho_{ c }(\vec{y})\f}

a \f$ Z(2) \times Z(2) \f$ colour vector at all (even) spatial sites \f$ \vec{y} \f$ and non-zero only on timeslice \f$ \tau \f$.

\f{equation}{\rho^{\alpha}_\beta = \delta_{\alpha \beta} \rho^\alpha\f}

eg.

\f{equation}{\rho^1_\beta = \left( \begin{matrix}
  0 \\
  \rho^1 \\
  0 \\
  0
 \end{matrix} \right)\f}
 
It then solves for the four objects
 
\f{equation}{\chi^\alpha_\beta = S_{\beta \gamma} \rho^{\alpha}_\gamma\f}

eg.

\f{equation}{\chi^0_\beta = \left( \begin{matrix}
  S_{00} \rho^0 \\
  S_{10} \rho^0 \\
  S_{20} \rho^0 \\
  S_{30} \rho^0
 \end{matrix} \right)\f}
 
For every different \f$ \Gamma \f$ that is required it does four more inversions,

\f{equation}{\chi^{\Gamma \alpha}_\beta = S_{\beta \gamma} (\Gamma \gamma_5)^\dagger_{\gamma \delta} \rho^{\alpha}_\delta\f}

before calculating the correlator as,

\f{equation}{C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\lambda = 0}^{3}\sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \chi^{\lambda}_{\beta, b}(\vec{x}, t; \tau) \chi^{\Gamma \lambda \dagger}_{\alpha, b}(\vec{x}, t; \tau)\f}

where the  \f$ \lambda \f$ sum is over the \f$ 4 \f$  spinor components.

We should be able to improve the signal and reduce the number of
inversions with two modifications. First, instead of having a different
noise vector for every spin component we reuse the same noise, i.e.

\f{equation}{\rho^{\alpha}_\beta = \delta_{\alpha \beta} \rho\f}

for fixed \f$ \rho \f$.
Using less noise seems to be generally preferred.

Secondly there is no need to invert for every different \f$ \Gamma \f$. Let,

\f{equation}{\chi^{\Gamma \alpha}_\beta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha} \chi^{\gamma}_\beta\f}

This is true because,

\f{equation}{(\Gamma \gamma_5)^\dagger_{\gamma \alpha} \chi^{\gamma}_\beta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \delta} \rho^{\gamma}_\delta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \delta} \delta_{\gamma \delta} \rho = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \gamma} \rho = S_{\beta \gamma} (\Gamma \gamma_5)^\dagger_{\gamma \alpha} \rho\f}

then the correlation function is

\f{equation}{C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\lambda = 0}^{3}\sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \chi^{\lambda}_{\beta, b}(\vec{x}, t; \tau) \chi^{\Gamma \lambda \dagger}_{\alpha, b}(\vec{x}, t; \tau)\f}

as before. By using the spin_matrix object in HiRep to construct the
objects \f$ \chi^{\lambda}_{\beta, b} \f$ the correlators can be calculated
with only \f$ 4 N_F \f$ inversions.

### Disconnected

The disconnected contributions occur when we have fermion species of the
same type in the hadron interpolator \f$ O_M \f$ :

\f{equation}{O_M(x) = \bar{ \psi } (x) \Gamma \psi(x)\f}

The same manipulations that lead to equation eq.\f$(\ref{eq:2ptcorr})\f$ give, 

\f{equation}{\langle O_M(x) \bar{O}_M'(y) \rangle = \langle \bar{\psi}(x) \Gamma \psi(x) \bar{\psi}(y) \bar{ \Gamma }' \psi(y) \rangle = \Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \bar{\psi}_\alpha(x)  \psi_\beta(y) \bar{\psi}_\gamma(y)  \psi_\delta(x) \rangle\f}

There are two allowed Wick contractions,

\f{equation}{\langle O_M(x) \bar{O}_M'(y) \rangle = -\text{Tr}\left[ \Gamma S (x,y) \bar{ \Gamma }' S (y,x) \right] + \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]\f}

the connected and disconnected contributions. Fourier transforming the
first term gives us the same result as before. For the disconnected part, 

\f{equation}{\begin{aligned}
D(t - \tau, \vec{p}) = \sum_{\vec{x} \vec{y}} e^{-i \vec{p} (\vec{x} - \vec{y} )} \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]\end{aligned}\f}

again \f$ x = (\vec{x}, t) \f$ and \f$ y = (\vec{y}, \tau) \f$, the zero-momentum
correlator is 

\f{equation}{D(t - \tau, 0) = \sum_{\vec{x} \vec{y}} \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]
 = \sum_{\vec{x}} \text{Tr}\left[ \Gamma S (x,x) \right] \sum_{\vec{y}} \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]\,.\f}
 
This means we have to evaluate objects like

\f{equation}{\begin{aligned}
d(t) = \sum_{\vec{x}} \text{Tr}\left[ \Gamma S (x,x) \right] = \sum_{\vec{x}} \Gamma_{\alpha \beta} S_{\beta \alpha}(x,x)\,,\end{aligned}\f}

using

\f{equation}{\phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) = \sum_{\vec{y}} S_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle\,,\f}

which implies 

\f{equation}{\begin{aligned}
\frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \text{Tr} \left[ \langle \eta^{(k)}(\vec{x}) | \Gamma | \phi^{(k)} (\vec{x}, t; \tau) \rangle \right] =& \frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \Gamma_{\beta \gamma} | \phi^{(k)}_{\gamma, b}(\vec{x}, t; \tau) \rangle \langle \eta^{(k)}_{\beta b }(\vec{x}) | \\
 =& \frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \Gamma_{\beta \gamma} \sum_{\vec{y}} S_{\gamma \alpha,b c} (\vec{x}, t; \vec{y}, \tau) 
| \eta^{(k)}_{\alpha c }(\vec{y})\rangle \langle \eta^{(k)}_{\beta b }(\vec{x}) |\end{aligned}
\f}

Using the limit \f$ K \rightarrow \infty \f$ this becomes, 

\f{equation}{\sum_{\vec{x}} \Gamma_{\beta \gamma} \sum_{\vec{y}} S_{\gamma \alpha,b c} (\vec{x}, t; \vec{y}, \tau) \delta_{bc} \delta_{\alpha \beta} \delta(\vec{y}, \vec{x}) =
\sum_{\vec{x}} \Gamma_{\beta \gamma} S_{\gamma \beta,b b} (\vec{x}, t; \vec{x}, \tau)\,,\f}

\f{equation}{d(t, \tau) = \text{Tr} \left[ \Gamma S(\vec{x}, t; \vec{x}, \tau) \right]\,.\f}

We want only cases where \f$ t = \tau \f$ so we need either four noise
vectors on every timeslice or noise vectors that are nonzero on all
timeslices. In the latter case we would evaluate, 

\f{equation}{\begin{aligned}\frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \text{Tr} \left[ \langle \eta^{(k)}(\vec{x}, t) | \Gamma | \phi^{(k)} (\vec{x}, t) \rangle \right] \end{aligned}\f}

with

\f{equation}{\phi^{(k)}_{\beta, b}(\vec{x}, t) = \sum_{\vec{y}} S_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y}, \tau)\,.\rangle\f}

### Cancelling Backwards Propagation

The two-point function evaluated in the center of the lattice is
(including the backward propagating part to give the extra factor of
$2$),

\f{equation}{C(T/2, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   2 e^{- E_\pi(\vec{p}) (T/2) }\f}

therefore

\f{equation}{\frac{1}{2} C(T/2, \vec{p}) e^{ - E_\pi(\vec{p}) (T/2 - t)  } = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) (T - t) }\f}

then

\f{equation}{C_{\rightarrow}(t, \vec{p}) = C(t, \vec{p}) - \frac{1}{2} C(T/2, \vec{p}) e^{ - E_\pi(\vec{p}) (T/2 - t)  }\f}

is the forward propagating part only.  \f$ E_\pi(\vec{p}) \f$  is obtained by
fitting the zero momentum correlator and using
 \f$ E(\vec{p}) = \sqrt{m_\pi^2 + \vec{p}^2} \f$. The factor  \f$ C(T/2, \vec{p}) \f$ 
can be obtained also from the zero momentum correlator, by fitting to
obtain  \f$ |Z_\pi|^2 \f$  and using the fact that this is momentum independent.
Since  \f$ 0 \f$  momentum results are used this might not be too noisy.

Alternatively the Wilson action is invariant under 

\f{equation}{\begin{aligned}\psi(x) \rightarrow {\cal P_\mu}[ \psi(x) ] = \gamma_\mu \psi( P_\mu[x]) \\
\bar{ \psi } (x) \rightarrow {\cal P_\mu}[ \bar{ \psi } (x) ] = \bar{ \psi }( P_\mu[x])  \gamma_\mu \end{aligned}\f}

where \f$ P_\mu[x] \f$ reverses the sign of all the components of \f$ x \f$ except
the \f$ \mu \f$ one. Time reversal corresponds to \f$ {\cal T} = {\cal P}_1{\cal P}_2{\cal P}_3 \f$. 

\f{equation}{\begin{aligned}\psi(x) \rightarrow {\cal T}[ \psi(x) ] = \gamma_0 \gamma_5 \psi( T[x]) \\
\bar{ \psi } (x) \rightarrow {\cal T}[ \bar{ \psi } (x) ] = \bar{ \psi }( T[x])  \gamma_5 \gamma_0 \end{aligned}\f}

Using this the T symmetry of operators used to construct the correlators
can be calculated to calculate the sign on the backwards propagating
part.

\f{equation}{\langle O_1(t) O_2(0) \rangle = C(t) = A \left( e^{ -Et} + \tau_1 \tau_2 e^{-E(T - t)}\right)\f}

Here \f$ \tau_i = \pm 1 \f$ is the \f$ {\cal T} \f$ eigenvalue of \f$ O_i \f$. We mostly use correlators where \f$ O_1 = O_2 \f$ so \f$ \tau_1 \tau_2 = 1 \f$ then the correlator is

\f{equation}{C_{pp}(t, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } + e^{- E_\pi(\vec{p})(T - t) } + e^{- E_\pi(\vec{p})(2T - t) } + \ldots \right)\f}

The subscript on \f$ C_{pp} \f$ refers to the fact that both propagators used
periodic boundary conditions. We want to cancel the backwards
propagating part which can be done by solving the forward propagator
 \f$ S(0,x) \f$  using antiperiodic time bc's and the backward  \f$ S(x,0) \f$  with
periodic time bc's to give an extra minus sign,

\f{equation}{C_{ap}(t, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } - e^{- E_\pi(\vec{p})(T - t) } + e^{- E_\pi(\vec{p})(2T - t) } - \ldots \right)\f}

so,

\f{equation}{C_{ap}(t, \vec{p}) + C_{p}(t, \vec{p}) = \frac{ 2 |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } + e^{- E_\pi(\vec{p})(2T - t) } + \ldots \right)\f}

cancelling the subleading exponential. This method requires two
inversions and the calculation of

\f{equation}{S_{A \pm P} (x,y) = S_A(x,y) \pm S_P(x,y)\f}

Where the subscript refers to (A)ntiperiodic/(P)eriodic boundary conditions. Then,

\f{equation}{C_{\pm}(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \gamma_5 \Gamma S_{A \pm P} (x,y) \bar{ \Gamma }' \gamma_5 S_{A \pm P} (x,y) \right]\f}

where  \f$ C_{+}(t - \tau, 0) \f$  gives the forward propagating part from  \f$ 0 \f$ 
to  \f$ T \f$  and  \f$ C_{-}(t - \tau, 0) \f$  gives the backwards propagating part
from  \f$ 2T \f$  to  \f$ T \f$.

### Form Factors and Sequential Sources

The electromagnetic form factor of a 'pion' requires the evaluation of
the matrix element

\f{equation}{\langle \pi(p_f) | V_\mu | \pi(p_i) \rangle = (p_i + p_f)_\mu f(q^2)\f}

where \f$ q^2 = (p_i - p_f)^2 \f$ and

\f{equation}{V_\mu = q_u \bar{u} \gamma_\mu u + q_d \bar{d} \gamma_\mu d\f}

is the electromagnetic current and \f$ q_i \f$ is the charge of the fermion \f$ i \f$. This
is the local (not conserved) current, so there will be a factor  \f$ Z_V \f$ 
for renormalization. The matrix elements required look like:

\f{equation}{\begin{aligned}
C_3(t_f,t,t_i,\vec{p_i}, \vec{p_f}) = Z_V \sum_{\vec{x} \vec{y} \vec{z}} e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} (\vec{y} - \vec{z)}}\langle 0 | \bar{u} \gamma_5 d ( \vec{x}, t_f)  V_0(\vec{y}, t) \bar{d} \gamma_5 u (\vec{z}, t_i) | 0 \rangle\end{aligned}\f}

We take the \f$ \mu = 0 \f$ component since this is statistically cleaner and
also nonzero independent of the momentum direction. The contractions
give three propagators eg. taking the  \f$ \bar{d} \gamma_\mu d \f$  part of
\f$V_\mu \f$, 

\f{equation}{\begin{aligned}
Z_V \sum_{\vec{x} \vec{y} \vec{z}} e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{-i\vec{p_i} (\vec{y} - \vec{z} ) } Tr \left[ S_u(\vec{z},t_i;\vec{x},t_f) \gamma_5 S_d(\vec{x},t_f;\vec{y},t)  \gamma_0 S_d(\vec{y},t;\vec{z},t_i) \gamma_5  \right]\end{aligned}\f}

There are also disconnected contributions from contracting the two
fermions in the current together but we ignore those. The usual
sequential source trick consists of solving

\f{equation}{S_u(\vec{x}, t; \vec{z}, t_i) = \sum_{\vec{y},\tau} D_u ( \vec{x}, t; \vec{y},\tau)  \delta(\vec{y}, \tau; \vec{z}, t_i)\f}

to get the point-to-all propagator (for a specific  \f$ \vec{z} \f$  and  \f$ t_i \f$ 
as well as dropping the sum over  \f$ \vec{z} \f$  and using translational
invariance). Then taking a single timeslice of the propagator
 \f$ S_u(\vec{x}, t_f; \vec{z}, t_i) \f$  and solving, 

\f{equation}{D_d ( \vec{x},t_f; \vec{y}, t ) G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i) = e^{i\vec{p_f} \vec{x}} \gamma_5 S_u(\vec{x}, t_f; \vec{z}, t_i)\f}


\f{equation}{G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i) = \sum_{ \vec{x} } e^{i\vec{p_f} \vec{x}} S_d(\vec{y},t; \vec{x}, t_f ) \gamma_5 S_u(\vec{x}, t_f; \vec{z}, t_i)\f}

to get the all-to-all-to-point contribution. Then 

\f{equation}{\begin{aligned}\gamma_5 \left[ G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i)  \right]^\dagger \gamma_5 =& \sum_{ \vec{x} } e^{-i\vec{p_f} \vec{x}} \gamma_5 S_u^\dagger (\vec{x}, t_f;\vec{y},t ) \gamma_5 \gamma_5 \gamma_5 \gamma_5 \gamma_5 S_d^\dagger (\vec{z}, t_i;\vec{x}, t_f) \gamma_5 \\
=& \sum_{ \vec{x} } e^{-i\vec{p_f} \vec{x}} S_u (\vec{z}, t_i;\vec{x}, t_f) \gamma_5 S_d (\vec{x}, t_f;\vec{y},t ) \\
=& G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t )\end{aligned}\f}

\f{equation}{\begin{aligned}
C_3(t_f,t,t_i, \vec{p_i}, \vec{p_f}) = Z_V \sum_{\vec{y} } e^{-i(\vec{p_i} - \vec{p_f})\vec{y}} Tr \left[ G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t )  \gamma_0 S_d(\vec{y},t;\vec{z},t_i) \gamma_5 \right]\end{aligned}\label{eq:3ptfn}\f}

This \f$ Z_V \f$ factor is unknown. We show how to calculate it later, or
cancel it, but an alternative is to use the conserved vector current in
place of the local current

\f{equation}{V_\mu = \frac{1}{2} \left[ \bar{\psi}(x + \mu)(1 + \gamma_\mu)U_\mu^\dagger(x) \psi(x) - \bar{\psi}(x)(1 - \gamma_\mu)U_\mu^\dagger(x) \psi(x + \mu) \right]\,.\f}

The trace in eq.\f$(\ref{eq:3ptfn})\f$ becomes 

\f{equation}{\begin{aligned}Tr[ S_d(\vec{y},t+1;\vec{z},t_i) &\gamma_5 (1 + \gamma_0)U_0^\dagger(\vec{y},t) G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t ) -\\ &S_d(\vec{y},t;\vec{z},t_i) \gamma_5 (1 - \gamma_0)U_0(\vec{y},t) G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t+1 ) ]\end{aligned}\f}

If we use this then all the following formulas are the same except \f$ Z_V \rightarrow 1 \f$.

There is an alternative that doesn't require the sequential source
trick. Using the properties of our noise sources

\f{equation}{S_d(\vec{x},t_f;\vec{y},t) \approx \frac{1}{K} \sum_{i=0}^K | \psi^{(i)}(\vec{x},t_f) \rangle \langle \eta^{(i)}(\vec{y},t) |\f}

the three point correlation function becomes

\f{equation}{\begin{aligned}
Z_V \frac{1}{K} \sum_{i=0}^K \sum_{\vec{x} \vec{y} } e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} \vec{y}} Tr \left[  \langle \eta^{(i)}(\vec{y},t) | \gamma_0 S_d(\vec{y},t;\vec{0},0) \gamma_5  S_u(\vec{0},0;\vec{x},t_f) \gamma_5 | \psi^{(i)}(\vec{x},t_f) \rangle \right] \\
Z_V \frac{1}{K} \sum_{i=0}^K \sum_{\vec{x} \vec{y} } e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} \vec{y}} Tr \left[  \langle \eta^{(i)}(\vec{y},t) | \gamma_0 S_d(\vec{y},t;\vec{0},0) S_u^{\dagger}(\vec{x},t_f;\vec{0},0) | \psi^{(i)}(\vec{x},t_f) \rangle \right]\end{aligned}\f}

Using this method we can inject arbitrary momentum at the source without
the need for extra inversions.

#### Two-Point Functions

A complete set of hadrons is given by,

\f{equation}{\sum_n \frac{ | n \rangle \langle n |}{ 2 E_n V}\f}

the first term is the pion. The two-point function (from point sources) is,

\f{equation}{C(t, \vec{p}) = \sum_{\vec{x}} e^{-i \vec{p} \vec{x}} \langle O_\pi (\vec{x}, t) O^\dagger_\pi (\vec{0}, 0)\rangle 
= \sum_{\vec{x}} \sum_n e^{i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, t)  | n \rangle \langle n | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_n }\f}

Use \f$ \sum_{ \vec{y} } e^{-i\vec{p} \vec{y}} O^\dagger_n (\vec{y}, 0) | 0 \rangle  = | n(\vec{p}) \rangle \f$, the time evolution operator \f$ e^{-Ht} \f$ and also the fact that the lightest meson dominates the sum to get

\f{equation}{\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} \vec{z}} \sum_{\vec{p'} } e^{-i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, 0)  O^\dagger_\pi (\vec{y}, 0) e^{i \vec{p'} \vec{y} }   | 0 \rangle \langle 0 | e^{i \vec{p'} \vec{z} } O_\pi (\vec{z}, 0) O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p'}) } e^{- E_\pi(\vec{p'}) t }\end{aligned}\f}

The sum over \f$ \vec{p'} \f$ gives a delta function leaving

\f{equation}{\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} } e^{-i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, 0)  O^\dagger_\pi (\vec{y}, 0) | 0 \rangle \langle 0 | O_\pi (\vec{y}, 0) O^\dagger_\pi (\vec{0}, 0) | 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}\f}

Translational invariance lets us write

\f{equation}{\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} } e^{-i \vec{p} ( \vec{x} - \vec{y} ) } e^{ -i \vec{p} \vec{y} }  
\frac{ \langle 0 | O_\pi (\vec{0}, 0)  O^\dagger_\pi (\vec{x}-\vec{y}, 0) | 0 \rangle \langle 0 | O_\pi (\vec{y}, 0) O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}\f}

Now changing variables gives us two Fourier transforms

\f{equation}{\begin{aligned}
C(t, \vec{p}) =& \frac{ \langle 0 |  O_\pi (\vec{0}, 0) | \pi(p) \rangle \langle \pi(p) | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}\f}

and finally using the time evolution operator we get

\f{equation}{\begin{aligned} C(t, \vec{p}) =& \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}\f}

where 

\f{equation}{Z_\pi = \langle \pi(p) | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle\f}

#### Three-Point Functions

In less detail we insert two complete sets of states into the correlator
( point sources so \f$ (\vec{x_i}, t_i) = (\vec{0}, 0) \f$ ) 

\f{equation}{\begin{aligned}\langle \pi(p_f) | V_\mu | \pi(p_i) \rangle &=\, \langle 0| O(\vec{p_f}, t_f) V_\mu(\vec{x}, t) O^\dagger(\vec{p_i}, t_i) |0\rangle \\
 &=\, \langle 0| O(\vec{0}, 0) | \pi(\vec{p_f}) \rangle \frac{e^{-(t_f - t) E_\pi(\vec{p_f}) }}{2 E_\pi(\vec{p_f}) } \langle \pi(\vec{p_f}) | V_\mu(\vec{0}, 0) | \pi(\vec{p_i}) \rangle \times \frac{e^{-(t - t_i) E_\pi(\vec{p_i}) }} {2 E_\pi(\vec{p_i}) } \langle \pi(\vec{p_i}) | O^\dagger(\vec{0}, 0) |0 \rangle \\ \nonumber
&=\,  \frac{ Z_{\pi, f}^\dagger Z_{\pi, i} }{4 E(\vec{p_f}) E(\vec{p_i}) } \langle \pi(\vec{p_f}) | V_\mu(\vec{0}, 0) | \pi(\vec{p_i}) \rangle  e^{-(t_f - t) E_\pi(\vec{p_f}) -(t-t_i) E_\pi(\vec{p_i}) }\end{aligned}\f}

if \f$ t < t_f \f$ we have the backwards contribution and the exponential
changes to

\f{equation}{-e^{-(t - t_f) E_\pi(\vec{p_f}) -(T - t + t_i) E_\pi(\vec{p_i}) }\f}

#### Correlator Ratios: \f$Z_V\f$

 \f$ Z_V \f$  can be obtained as follows: The ratio, 

\f{equation}{\begin{aligned}\frac{ C_{\rightarrow}(t_f, \vec{0}) }{ C_3(t_f, t, \vec{p_i}, \vec{p_f}) } =&  
\frac{ \frac{ |Z_\pi( \vec{0} )|^2 }{ 2 m_\pi }   e^{- m_\pi t_f }  }{ \frac{ |Z_\pi( \vec{0} )|^2 }{4 m_\pi^2 } \langle \pi(\vec{0}) | V_\mu | \pi(\vec{0}) \rangle e^{-(t_f - t) m_\pi -t m_\pi } } \\
=& \frac{ 1  }{ \frac{ 1 }{2 m_\pi } \langle \pi(\vec{0}) | V_\mu | \pi(\vec{0}) \rangle } \\
=& \frac{ 1  }{ \frac{ 1 }{2 m_\pi } 2 m_pi f(0)/Z_V } = Z_V.\end{aligned}\f}

Where we used that the renormalized form factor \f$ f(0) = 1 \f$ 

#### Correlator Ratios: \f$ f(q) \f$ 

There are various ways to cancel the unwanted terms and get \f$ f(q) \f$.

##### RBC-UKQCD Ratio

We examine the ratio,

\f{equation}{2 m_\pi \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{3} (t, t_f, \vec{0}, \vec{0} )  C_{\rightarrow}(t, \vec{p}) }\f}

Assuming \f$ Z_\pi \f$ is momentum independent. This also works if
 \f$ Z_\pi = E(\vec{p}) f_\pi \f$, which is the case for
 \f$ O = \bar{u} \gamma_0 \gamma_5 d \f$  type interpolators. The numerator is,

\f{equation}{\frac{ Z_V |Z_\pi|^2 }{ 4 E(\vec{p}) E(\vec{0})} f(q^2) ( E(\vec{p}) + m_\pi ) \frac{|Z_\pi|^2}{2 E(\vec{0})} e^{ -E(\vec{p})t - E(\vec{0})(t_f - t) -E(\vec{0})t }\f}

and the denominator is,

\f{equation}{\frac{ Z_V |Z_\pi|^2 }{ 4 E(\vec{0}) E(\vec{0})} f(0) ( m_\pi + m_\pi ) \frac{|Z_\pi|^2}{2 E(\vec{p})} e^{ -E(\vec{0})t - E(\vec{0})(t_f - t) -E(\vec{p})t }\f}

Cancelling leaves,

\f{equation}{2 m_\pi \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{3} (t, t_f, \vec{0}, \vec{0} )  C_{\rightarrow}(t, \vec{p}) } = f(q^2) ( E(\vec{p}) + m_\pi )\f}

note there is no \f$ Z_V \f$  here.

##### Bonnet et. al. Ratio

\f{equation}{\frac{2 Z_V m_\pi}{E(\vec{p}) + m_\pi} \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{\rightarrow} (t, \vec{p} )  C_{\rightarrow}(t_f, \vec{0}) }\f}

the numerator of the right term is,

\f{equation}{\frac{ |Z_\pi|^2 }{ 4 E(\vec{p}) m_\pi} f_B(q^2) ( E(\vec{p}) + m_\pi ) \frac{|Z_\pi|^2}{2 m_\pi} e^{ -E(\vec{p})t - m_\pi(t_f - t) -m_\pi t }\f}

the denominator of the right term is,

\f{equation}{\frac{ |Z_\pi|^2 }{ 2 m_\pi } \frac{|Z_\pi|^2}{2 E(\vec{p})} e^{ -m_\pi t - m_\pi (t_f - t) -E(\vec{p})t }\f}

Cancelling leaves,

\f{equation}{\frac{ f_B(q^2) ( E(\vec{p}) + m_\pi ) }{ 2 E(\vec{p}) }\f}

the kinematic factors are cancelled

\f{equation}{\frac{2 Z_V m_\pi}{E(\vec{p}) + m_\pi} \frac{ f_B(q^2) ( E(\vec{p}) + m_\pi ) }{ 2 E(\vec{p}) } = Z_V f_B(q^2) = f(q^2)\f}

you need to actually know  \f$ Z_V \f$  or use the conserved current.


## Estimation of Disconnected Contributions

### Conventions

We choose the hermitian basis of gamma matrices given in Tab. 1. Each element of the basis is referred by an index in \[0,15\] shown in the following table

| No | Matrix			  					|
|:---|:---------------------------------------|
| 0  | \f$ \gamma_5 \f$ 		  					|
| 1  | \f$ \gamma_1 \f$ 		  					|
| 2  | \f$ \gamma_2 \f$ 		  					|  
| 3  | \f$ \gamma_3 \f$ 		  					|
| 4  | \f$ -\mathrm{i}\gamma_0\gamma_5 \f$ 			|
| 5  | \f$ -\mathrm{i}\gamma_0\gamma_1 \f$ 			|
| 6  | \f$ -\mathrm{i}\gamma_0\gamma_2 \f$ 			|
| 7  | \f$ -\mathrm{i}\gamma_0\gamma_3 \f$ 	    		|
| 8  | 1										|
| 9  | \f$ -\mathrm{i}\gamma_5\gamma_1 \f$ 			|
| 10 | \f$ -\mathrm{i}\gamma_5\gamma_2 \f$ 			|
| 11 | \f$ -\mathrm{i}\gamma_5\gamma_3 \f$ 			|
| 12 | \f$ \gamma_0 \f$ 							|
| 13 | \f$ -\mathrm{i}\gamma_5\gamma_0\gamma_1 \f$ 	|
| 14 | \f$ -\mathrm{i}\gamma_5\gamma_0\gamma_2 \f$  |
| 15 | \f$ -\mathrm{i}\gamma_5\gamma_0\gamma_3 \f$  |

### Singlet Two-Point Functions

Consider a gauge theory on a group G coupled to  \f$ N_f \f$  fermions in an arbitrary representation  \f$ R \f$. Let us denote:

\f{equation}{C(t, x_0) = \dfrac{1}{N_f}\sum_{\vec{x}}\langle \bar{q}\Gamma q(x)\bar{q}\Gamma q(x_0)\rangle\f}

where  \f$ q \f$, \f$\bar{q} \f$  are the  \f$ N_f \f$  quark fields and  \f$ \Gamma \f$  denotes as arbitrary Dirac structure. The  \f$ 1/N_f \f$  factor is only there for convenience. The Wick contractions read:

\f{equation}{C(t, x_0) = \sum_{\vec{x}} \langle -\mathrm{tr}\left(\Gamma S(x,x_0)\Gamma S(x_0,x)\right) + N_f\,\mathrm{tr}\left(\Gamma S(x,x)\right)\mathrm{tr}\left(\Gamma S(x_0,x_0)\right)\rangle\f}

### Stochastic Evaluation of Disconnected Loops

The simple one consist to evaluate stochastically the disconnected contribution without any variance reduction techniques. Considering a general volume source  \f$ \xi \f$, we define  \f$ \phi \f$  using the Dirac operator  \f$ D \f$ :

\f{equation}{\phi = D^{-1}\xi\f}

For a given element X of the basis defined in the previous section, we then have

\f{equation}{\sum \left(\xi^{*}X\phi\right)_{R} = \sum XM^{-1} + \text{noise}\f}

where the symbol  \f$ (\ldots)_{R} \f$  refers to the average over R samples of the stochastic source. 

It should be observed that in evaluating the disconnected contributions to the neutral meson correlators each one of the two quark loops arising from Wick contractions must be averaged over completely independent samples of stochastic sources for the purpose of avoiding unwanted biases.

#### Implemented Source Types

TODO: XX
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

TODO: Notation

`#T #iGamma #iSrc #\[color and/or e/o \] #Re #Im`

where iGamma refers to the index of the Gamma matrix defined in Table 1.

### Debugging Options

If the code is executed with the following additional arguments

```bash
-p  <propagator_name> -s <source_name>
```

This will read the two files and perform the contraction accordingly computing \f$ \chi^{\dagger}\Gamma \psi \f$.

## Mesonic Correlators of the Isotriplet

The two fermionic flavors are denoted by \f$ u \f$ and \f$ d \f$. We are interested in the mesonic correlators

\f{equation}{\begin{aligned}
&& C_{\Gamma_1,\Gamma_2}(x-y) = \left<
\left( \bar{u} \Gamma_1 d \right)^\dagger(x)
\left( \bar{u} \Gamma_2 d \right)(y)
\right> = \nonumber \\
&& \quad = \left<
\left( \bar{d} \gamma_0 \Gamma_1^\dagger \gamma_0 u \right)(x)
\left( \bar{u} \Gamma_2 d \right)(y)
\right>\end{aligned}\f}

where \f$ \Gamma_i \f$ are the generic products of the \f$ \gamma\f$-matrices.

We can integrate out the fermionic fields explicitly. Here we use the definition \f$ H(x,y) = G(x,y) \gamma_5 \f$  for the hermitian Dirac operator, with \f$ G(x,y) \f$ defined as its inverse.

\f{equation}{\begin{aligned}
&& C_{\Gamma_1,\Gamma_2}(x-y) = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,y) \Gamma_2 G(y,x) \right]
\right> = \nonumber \\
&& \quad = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,y) \Gamma_2 \gamma_5 G(x,y)^\dagger \gamma_5 \right]
\right> = \nonumber \\
&& \quad = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 H(x,y) \gamma_5 \Gamma_2 H(x,y)^\dagger \gamma_5 \right]
\right> \; .\end{aligned}\f}

Since the \f$ \gamma\f$-matrices commute, we can conclude that the matrix \f$ \gamma_0 \Gamma^\dagger \gamma_0 \f$ is equal to \f$ \Gamma \f$ up to a sign

\f{equation}{\label{eq:gamma0_adj}\gamma_0 \Gamma^\dagger \gamma_0 = s(\Gamma) \Gamma \qquad \text{with } s(\Gamma) = \pm 1 \; .\f}

In addition, a generic matrix \f$ \Gamma \f$ has the following properties:

1.	Its matrix elements can be \f$ 0 \f$, \f$ \pm 1 \f$, \f$ \pm i \f$ 
2.	Its entries are either all real or all imaginary
3.	In each row and correspondingly each column, there is only one non-zero element

Consequently, we can write

\f{equation}{\label{gamma_ab}
\Gamma_{\alpha\beta} = t_\alpha(\Gamma) \delta_{\sigma_\alpha(\Gamma), \beta}\f}

where \f$ \sigma(\Gamma) \f$ constitutes a permutation of four elements. Putting this together we find

\f{equation}{\begin{aligned}
& C_{\Gamma_1,\Gamma_2}(x-y) = - s(\Gamma_1) \left< \mathrm{tr}
\left[ \gamma_5 \Gamma_1 H(x,y) \gamma_5 \Gamma_2 H(x,y)^\dagger \right]
\right> = \\
& \quad = - s(\Gamma_1) \sum_{\alpha\beta} t_\alpha(\gamma_5 \Gamma_1) t_\beta(\gamma_5 \Gamma_2) \times \left< \mathrm{tr}
\left[ H_{\sigma_\alpha(\gamma_5 \Gamma_1), \beta}(x,y) H_{\alpha, \sigma_\beta(\gamma_5 \Gamma_2)}(x,y)^\dagger \right]
\right> \; .\label{triplet_corr}\end{aligned}\f}

### Implementation of the Point-To-All Propagator

In order to calculate mesonic masses we are interested in correlators satisfying \f$ \Gamma_1=\Gamma_2 \f$. Using translational invariance, we can set \f$ y=0 \f$. In this case the formula simplifies to

\f{equation}{C_{\Gamma}(x) = - s(\Gamma) \sum_{\alpha\beta} t_\alpha(\gamma_5 \Gamma)t_\beta(\gamma_5 \Gamma) \times \left< \mathrm{tr}\left[ H_{\sigma_\alpha(\gamma_5 \Gamma), \beta}(x,0) H_{\alpha, \sigma_\beta(\gamma_5 \Gamma)}(x,0)^\dagger \right]\right> \; .\label{eq:triplet_point_to_all_corr}\f}

This is implemented into HiRep in the following way

* The data of the point-like source \f$ \xi^{(\alpha,a)} \f$ defined by

   \f{equation}{\xi^{(\alpha,a)}_{\beta b}(x) = \delta_{\alpha,\beta} \delta_{a,b} \delta_{x,0} \; ,\f}
	
    The function `quark_propagator` applies the inverse of the hermitian Dirac operator to the source
	
   \f{equation}{\eta^{(\alpha,a)} = H \xi^{(\alpha,a)} \qquad \eta^{(\alpha,a)}_{\beta b}(x) = H_{\beta \alpha}^{b a}(x,0) \; .\f}

-   The functions `void *_correlator(float *out, suNf_spinor **qp)`in `Observables/mesons.c` implement the formulae eq.\f$(\ref{eq:triplet_point_to_all_corr})\f$, where `out` stands for the correlator and `qp` for the spinor array. The functions \f$ s(\Gamma) \f$, \f$ t_\alpha(\gamma_5 \Gamma) \f$ and \f$ \sigma_\alpha(\gamma_5 \Gamma) \f$ where calculated using `Mathematica`, see file `mesons.nd` and implemented in the code using macros, defined as follows\
    
    `      `\
    `      _C1_ = `\f$\sigma_1(\gamma_5 \Gamma)\f$\
    `      _C2_ = `\f$\sigma_2(\gamma_5 \Gamma)\f$\
    `      _C3_ = `\f$\sigma_3(\gamma_5 \Gamma)\f$\
    `      _C4_ = `\f$\sigma_4(\gamma_5 \Gamma)\f$\
    `      `
    
    If \f$ t_\alpha(\gamma_5 \Gamma) \f$ are read:
    
    `      `\
    `      _S0_ = `\f$-s(\Gamma)\f$\
    `      _S1_ = `\f$t_1(\gamma_5 \Gamma)\f$\
    `      _S2_ = `\f$t_2(\gamma_5 \Gamma)\f$\
    `      _S3_ = `\f$t_3(\gamma_5 \Gamma)\f$\
    `      _S4_ = `\f$t_4(\gamma_5 \Gamma)\f$\
    `      `
    
    whereas if \f$ t_\alpha(\gamma_5 \Gamma) \f$ are imaginary:
    
    `      `\
    `      _S0_ = `\f$s(\Gamma)\f$\
    `      _S1_ = `\f$-i t_1(\gamma_5 \Gamma)\f$\
    `      _S2_ = `\f$-i t_2(\gamma_5 \Gamma)\f$\
    `      _S3_ = `\f$-i t_3(\gamma_5 \Gamma)\f$\
    `      _S4_ = `\f$-i t_4(\gamma_5 \Gamma)\f$\
    `      `

## Mesonic Correlators of the Isosinglet

We are now concerned with the genertic mesonic correlator given by

\f{equation}{\begin{aligned}
C_{\Gamma_1,\Gamma_2}(x-y) &=& \left<
\left( \bar{u} \gamma_0 \Gamma_1^\dagger \gamma_0 u \right)(x)
\left( \bar{u} \Gamma_2 u \right)(y)
\right>\end{aligned}\f}

considering only a single flavor \f$ u \f$.

Integration of the fermionic fields now yiels one additional term, the hairpin diagram

\f{equation}{\begin{aligned}
C_{\Gamma_1,\Gamma_2}(x-y) =&- \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,y) \Gamma_2 G(y,x) \right]
\right>  + \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,x) \right]
\mathrm{tr}\left[ \Gamma_2 G(y,y) \right]
\right>\\
=& - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 H(x,y) \gamma_5 \Gamma_2 H(x,y)^\dagger \gamma_5 \right]
\right> + \left< \mathrm{tr}
\left[ \gamma_5 \gamma_0 \Gamma_1^\dagger \gamma_0 H(x,x) \right]
\mathrm{tr}\left[ \gamma_5 \Gamma_2 H(y,y) \right]
\right>
\; .\end{aligned}\f}

All other contributions are identical to the contributions to the isotriplet correlator. We can, therefore, focus on the contribution through the hairpin diagram. Using the formulae eq.\f$(\ref{eq:gamma0_adj})\f$ we can write

\f{equation}{\left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,x) \right]
\mathrm{tr}\left[ \Gamma_2 G(y,y) \right] \right> = \quad = s(\Gamma_1) \sum_{\alpha \beta} t_\alpha(\Gamma_1) t_\beta(\Gamma_2) \times \left<
\mathrm{tr}G_{\sigma_\alpha(\gamma_5 \Gamma_1),\alpha}(x,x) \;
\mathrm{tr}G_{\sigma_\beta(\gamma_5 \Gamma_2),\beta}(y,y)
\right>
\; .\label{eq:hairpin}\f}

### All-to-all Propagator

It is clear from eq.\f$(\ref{eq:hairpin})\f$, from the fact that we are employing point source, that one must compute the entire inverse matrix of the Dirac operator. The alternative is to use a statistic estimate for \f$ H \f$ followed by variance reduction procedures.

Suppose there are \f$ N_s \f$ available random fermion sources \f$ \xi^{(i)} \f$ such that the only non-zero correlators are 

\f{equation}{\left< \xi^{(i)}_{\alpha a}(x)^\dagger \xi^{(j)}_{\beta b}(y) \right> = \delta_{\alpha,\beta} \delta_{a,b} \delta_{x,y} \delta_{i,j} \; .\f}

Current literature proposes mainly either Gaussian noise or \f$ Z_2 \f$ noise. In the following, we will choose \f$ Z_2 \f$ noise, following @cite Foster_1999. Each component of the spinor is randomly chosen from the values \f$ \pm 1/\sqrt{2} \f$.

Then the matrix \f$ H \f$ can be estimated as follows:

\f{equation}{H_{\alpha \beta}^{a b}(x,y) \simeq \sum_{i=1}^{N_s} \eta^{(i)}_{\alpha a}(x) \xi^{(i)}_{\beta b}(y)^\dagger \eta^{(i)} \equiv H \xi^{(i)}\label{eq:naive_noisy_estimate}\f}

Stochastic estimation can then be used to calculate the relevant tracks for correlators

\f{equation}{\mathrm{tr}\left[ \Gamma_1 G(x,y) \Gamma_2 G(y,x) \right] = \sum_{ij} \xi^{(i)}(x)^\dagger \gamma_5 \Gamma_1 \eta^{(j)}(x) \times \xi^{(j)}(y)^\dagger \gamma_5 \Gamma_2 \eta^{(i)}(y)\f}

\f{equation}{\mathrm{tr}\left[ \Gamma G(x,x) \right] = \sum_i \xi^{(i)}(x)^\dagger \gamma_5 \Gamma \eta^{(i)}(x)\, .\f}

### Variance reduction

The noise obtained from stochastic estimation of the matrix \f$ G \f$ in the formula eq.`naive_noisy_estimate` can be reduced using the trick from @cite McNeile_2001 for Wilson fermions. Here, the Dirac operator has the form \f$ D = 1 - K \f$. As a result, for the matrix \f$ G \f$ the following formula applies 

\f{equation}{\begin{aligned}
&& G = D^{-1} = \left( 1 - K \right)^{-1} = \\
&& \quad = 1 + K + \dots + K^m + K^{n_1} G K^{n_2}\end{aligned}\f}

with \f$ n_1 + n_2 = n = m+1 \f$. In particular, for the evaluation of the hairpin diagram

\f{equation}{\begin{aligned}
\mathrm{tr}\left[ \Gamma G(x,x) \right] &=& \mathrm{tr}\Gamma + \mathrm{tr}\left[ \Gamma K^4 \right](x,x) + \mathrm{tr}\left[ \Gamma K^6 \right](x,x) + \dots + \\
&& + \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) + \mathrm{tr}\left[ \Gamma K^{n_1} G K^{n_2} \right](x,x)\end{aligned}\f}

with \f$ m=2k \f$. (TODO: fix this sentence) Here, we can use the fact that the matrix \f$ K \f$ connects first neighboring sites as thus \f$ K^p(x,x) \neq 0 \f$ only when \f$ p \f$ is even. Further, \f$ r_0=1 \f$ and consequently \f$ K^2(x,x) = 0 \f$ The first \f$ k+1 \f$ terms can be calculated explicitly and we can estimate the last term stochastically.

\f{equation}{\begin{aligned}
\mathrm{tr}\left[ \Gamma G(x,x) \right] =& \mathrm{tr}\Gamma + \mathrm{tr}\left[ \Gamma K^4 \right](x,x) + \mathrm{tr}\left[ \Gamma K^6 \right](x,x) + \dots + \nonumber \\
& + \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) + \sum_{iy} \xi^{(i)}(y)^\dagger \gamma_5 K^{n_2}(y,x) \Gamma K^{n_1}(x,y) \eta^{(i)}(y) \nonumber \\
\end{aligned}\label{eq:hairpin_with_variance_reduction}\f}

@cite McNeile_2001 use this trick only for the calculation of the hairpin diagram. It might be possible to generalize it to the the isotriplet part as well, as an alternative to the point-to-all propagator. 

### Time dilution

This is a trick introduced in @cite Foley_2005 for noise reduction in the computation of null-moment propagators. Whenever stochastic estimation of the \f$ H \f$ matrix is required, such as in eq.\f$(\ref{eq:naive_noisy_estimate})\f$, it is possible to replace each stochastic source \f$ \xi^{(i)} \f$ with a set of sources each with support on a different time slice.

\f{equation}{\begin{aligned}
 \label{time_dilution}
&& \xi^{(i)} \rightarrow \xi^{(i,\tau)}(\mathbf{x},t) = \xi^{(i)}(\mathbf{x},t) \delta_{t,\tau}\end{aligned}\f}

Stochastic estimation is now obtained similarly to the naive case:

\f{equation}{H_{\alpha \beta}^{a b}(x,y) \simeq \sum_{i=1}^{N_s} \sum_{\tau=1}^{N_t} \eta^{(i,\tau)}_{\alpha a}(x) \xi^{(i,\tau)}_{\beta b}(y)^\dagger  \eta^{(i,\tau)} \equiv H \xi^{(i,\tau)}\label{eq:diluted_noisy_estimate}\f}

### Implementation Scheme 

TODO: Add this to function reference instead if this is still implemented this way

The following functions will be implemented


-   Calculation of the exact terms of the formula eq.`hairpin_with_variance_reduction`
-   
    ` `\
    ` void GAMMA_variance_reduction_exact_terms(float *out, int k)`
    
    Here, ` out ` is a real vector with its components equal to the volume and \f$ k \f$ corresponds to the index in the \f$ \gamma\f$-matrix.
    
    This function evaluates
    
   \f{equation}{h_k(x) = \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) = \left( \frac{\kappa}{2} \right)^{2k} \sum_{\mathcal{C}_x} \mathrm{tr}\left( \Gamma \tilde{\gamma}(\mathcal{C}_x) \right) \mathcal{W}(\mathcal{C}_x)\f}
    
    where \f$ \mathcal{C}_x = (x,\hat{\mu}_1,\hat{\mu}_2,\dots,\hat{\mu}_{2k}) \f$ is the generic closed path of \f$ x \f$ of length \f$ 2k \f$ obtained by moving from \f$ x \f$ in the directions \f$ \hat{\mu}_i \f$. \f$ \mathcal{W}(\mathcal{C}_x) \f$ is the trace of the parallel transport through \f$ \mathcal{C}_x \f$ in the corresponding fermionic representation. \f$ \tilde{\gamma}(\mathcal{C}_x) \f$. is the matrix defined as
    
   \f{equation}{\tilde{\gamma}(\mathcal{C}_x) = (1-\gamma_{\hat{\mu}_1})(1-\gamma_{\hat{\mu}_2})\cdots(1-\gamma_{\hat{\mu}_{2k}})\f}
    
    having defined \f$ \gamma_{-\hat{\mu}_i} = - \gamma_{\hat{\mu}_i} \f$.

	It should be noted that since \f$ (1-\gamma_{\hat{\mu}_i})(1-\gamma_{-\hat{\mu}_i}) = 0 \f$, one can exclude the paths in which a pair of subsequences \f$ (\dots, \hat{\mu}_i,-\hat{\mu}_i, \dots) \f$ appears from the sum. In addition, the matrix \f$ \tilde{\gamma}(\mathcal{C}_x) \f$ does not depend on \f$ x \f$. There is is convenient to compute the list of paths and the matrix \f$ \tilde{\gamma}(\mathcal{C}_x) \f$ only once.

-   It is convenient to have a function that calculates the traces of the parallel transports:
-   
    ` void tr_r_pexp(complex *out, int *path, int length)`
    
    Here again ` out ` is the complex vector with number of components according to the volume and \f$ \phi(x) \f$ \ is the number of directions (`length`) of which the path is composed.
    
    This returns
    
   \f{equation}{\phi(x) = \mathrm{tr}\mathbf{R} \left[ U_{x, \hat{\mu}_1} U_{x+\hat{\mu}_1, \hat{\mu}_2} \cdots \right]\f}


-   Calculation of the time-diluted estimators
    
    ` void time_diluted_stochastic_estimate`\
    `      (suNf_spinor *csi, suNf_spinor **eta,`\
    `      int nm, float *mass, double ac)`\
    
    With the parameters:
    ` csi `  is the spinor \f$ \xi\f$\
    ` path ` is the list of \f$ N_t \times \textrm{nm} \f$ spinors along the \f$ \eta^{(\tau,m)}\f$
    ` `\
    
    This function generates the spinor \f$ \xi \f$ with \f$ Z_2 \f$ noise. Here the spinors are defined as 
    
   \f{equation}{\xi^{(\tau)}(\mathbf{x},t) = \xi(\mathbf{x},t) \delta_{t,\tau}\f}
    
   and then returned as 
    
   \f{equation}{\eta^{(\tau,m)} \equiv H_m \xi^{(\tau)}\f}
    
    using the inverter of the Dirac operator with parameters `nm`, `mass` and `acc`.
    
    The calculation of the hairpin term

    ` `\
    ` void GAMMA_hairpin_VR_TD(float **out,`\
    `      int n1, int n2, int nrs,`\
    `      int nm, float *mass, double acc)`\
    
    Takes the parameters
    		`out` = list of real `nm` vectors with the number of components corresponding to the volume 
    		` n1, n2 ` = parameters of the variance reduction algorithm
    		` nrs ` = number of random sources for the statistical estimation
    		` nm, mass ` = number and list of masses
    		`acc` = inverter parameter
	    ` `\
	    

The functions above implement the formula eq.\f$(\ref{eq:hairpin_with_variance_reduction})\f$, summing exact terms and the statistical term, generated with `nrs` to dilute, for a total of \f$ \times N_t \f$ matrix invertions for each mass value and returns the result as `out`.


