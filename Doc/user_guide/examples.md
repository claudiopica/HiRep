# Examples

## Contractions

### Connected two point correlation functions

First we choose interpolating operators with the quantum numbers of the
meson we would like to study:

$$O_M = \bar{\psi}^{(1)}(x) \Gamma \psi^{(2)}(x)$$

$$\bar{O}_M = \bar{\psi}^{(2)}(x) \bar{ \Gamma } \psi^{(1)}(x)$$ 

where $\bar{ \Gamma } = \gamma_0 \Gamma^\dagger \gamma_0$. The gamma matrix is
chosen to have the same $J^{PC}$ quantum numbers as the meson we want to
create. We then calculate the expectation value to create a meson at $y$
and destroy it again at $x$: 

$$\begin{aligned}
\langle O_M(x) \bar{O}_M'(y) \rangle =& \langle \bar{\psi}^{(1)}(x) \Gamma \psi^{(2)}(x) \bar{\psi}^{(2)}(y) \bar{ \Gamma }' \psi^{(1)}(y) \rangle 
\\
=& \Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \bar{\psi}^{(1)}_\alpha(x)  \psi^{(2)}_\beta(y) \bar{\psi}^{(2)}_\gamma(y)  \psi^{(1)}_\delta(x) \rangle \\
=& -\Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \psi^{(2)}_\beta(x) \bar{\psi}^{(2)}_\gamma(y)  \psi^{(1)}_\delta(y) 
\bar{\psi}^{(1)}_\alpha (x) \rangle\end{aligned}$$ 

where the sign in the last line comes from exchanging the anticommuting $\bar{\psi}_\alpha$ three times. Wick contracting where $\langle \psi_\alpha(x) \bar{\psi}_\beta(y) \rangle = S_{\alpha \beta}(x,y)$ gives 

$$\begin{aligned}
\langle O_M(x) \bar{O}_M'(y) \rangle =& -\Gamma_{\alpha \beta} S^{(2)}_{\beta \gamma} (x,y) \bar{ \Gamma }'_{\gamma \delta} S^{(1)}_{\delta \alpha} (y,x) \\
=& -\text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]\end{aligned}$$

To get correlation functions we Fourier transform, ie. sum over all
source and sink seperations while projecting onto the momentum we want
to give to the particle: 

$$\begin{aligned}
C(t - \tau, \vec{p}) =& \sum_{\vec{x} \vec{y}} e^{-i \vec{p}(\vec{x} - \vec{y}) } \langle O_M(\vec{x}, t) O_M'(\vec{y}, \tau)\rangle \\
=& -\sum_{\vec{x} \vec{y}} e^{-i \vec{p}(\vec{x} - \vec{y}) } \text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]\end{aligned}$$

here $x = (\vec{x}, t)$ and $y = (\vec{y}, \tau)$ and the zero momentum
correlator is:

$$C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \Gamma S^{(2)} (x,y) \bar{ \Gamma }' S^{(1)} (y,x) \right]$$

We now use $\gamma_5$ Hermiticity: $\gamma_5 S^\dagger(x,y) \gamma_5 = S(y,x)$, 

$$C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \gamma_5 \Gamma S^{(2)} (x,y) \bar{ \Gamma }' \gamma_5 S^{\dagger (1)} (x,y) \right]$$(eqn:corr)

### Point Sources

Using a delta function source and solving the Dirac equation gives a
point propagator,

$$D_{\alpha a, \beta b} (x,y) S_{\beta b, \gamma c} (y, z) = \delta(x,z) \delta_{ac} \delta_{\alpha \gamma}$$

usually $z = (\vec{0}, 0)$ so we get
$S(y,0) = \gamma_5 S^{\dagger} (0,y) \gamma_5$. Then we use these to
calculate correlation functions, 

$$C(t, 0) = -\sum_{\vec{x}} \text{Tr} e^{-i \vec{p} \vec{x}} \left[ \gamma_5 \Gamma S (x,0) \bar{ \Gamma }' \gamma_5 S (x,0) \right]$$(eqn:pointcorr)

Translational invariance in the limit of infinitely many gauge
configurations implies $S(x,y) = S(|x - y|)$, so the sum over $\vec{y}$
in equation {eq}`eq:corr` just gives $V$ times equation {eq}`eq:pointcorr`. We place the source at the time origin so $\tau = 0$.

### One-end Trick

For this method it helps to write all the indices out,

$$C(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,bc} (x,y) (\bar{ \Gamma }' \gamma_5)_{\gamma \delta} S^{\dagger (1)}_{\delta \alpha, cb} (x,y)$$

Greek indices $\alpha, \beta, \gamma, \delta, ...$ are spinor indices
and Latin indices $a, b, c, d...$ are colour indices. The one-end trick
involves inserting a delta function in colour, spin and space.

$$C(t - \tau, 0) = -\sum_{\vec{x} \vec{y} \vec{z}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) \delta_{\gamma \lambda} \delta_{cd} \delta(\vec{y}, \vec{z}) (\bar{ \Gamma }' \gamma_5)_{\lambda \delta} S^{\dagger (1)}_{\delta \alpha, d b} (\vec{x}, t;\vec{z}, \tau)$$

The delta function is aproximated with a $Z(2) \times Z(2)$ noise source on timeslice $\tau$

$$\delta_{\gamma \lambda} \delta_{cd} \delta(\vec{y}, \vec{z}) \approx \frac{1}{K} \sum_{k = 0}^{K} | \eta^{(k)}_{\gamma c }(\vec{y})\rangle \langle \eta^{(k)}_{\lambda d }(\vec{z}) |$$

which is exact in the limit $K \rightarrow \infty$. 

$$
C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\vec{x} \vec{y} \vec{z}} (\gamma_5 \Gamma)_{\alpha \beta} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle \langle \eta^{(k)}_{\lambda d }(\vec{z}) | (\bar{ \Gamma }' \gamma_5)_{\lambda \delta} S^{\dagger (1)}_{\delta \alpha, d b} (\vec{x}, t;\vec{z}, \tau)$$(eq:oet)

Defining

$$\phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) = \sum_{\vec{y}} S^{(2)}_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle$$

and

$$\phi^{\Gamma (k)}_{\alpha, b}(\vec{x}, t; \tau) = \sum_{\vec{z}} S^{(1)}_{\alpha \delta , b d} (\vec{x}, t;\vec{z}, \tau) (\bar{ \Gamma }' \gamma_5)^{\dagger}_{\delta \lambda} | \eta^{(k)}_{\lambda d }(\vec{z}) \rangle$$

the correlator can be evaluated as,

$$C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) \phi^{\Gamma \dagger (k)}_{\alpha, b}(\vec{x}, t; \tau)$$

#### Implementation in HiRep

In HiRep the code ```Spectrum/mk_mesons_with_z2semwall.c``` does two solves to
calculate $S^{(1)} | \eta \rangle$ and $S^{(2)} (\bar{ \Gamma }' \gamma_5)^{\dagger} | \eta \rangle$. HiRep has

$$\rho = \rho_{ c }(\vec{y})$$ 

a $Z(2) \times Z(2)$ colour vector at all (even) spatial sites $\vec{y}$ and non-zero only on timeslice $\tau$.

$$\rho^{\alpha}_\beta = \delta_{\alpha \beta} \rho^\alpha$$ 

eg.

$$\rho^1_\beta = \left( \begin{matrix}
  0 \\
  \rho^1 \\
  0 \\
  0
 \end{matrix} \right)$$ 
 
It then solves for the four objects
 
$$\chi^\alpha_\beta = S_{\beta \gamma} \rho^{\alpha}_\gamma$$ 

eg.

$$\chi^0_\beta = \left( \begin{matrix}
  S_{00} \rho^0 \\
  S_{10} \rho^0 \\
  S_{20} \rho^0 \\
  S_{30} \rho^0
 \end{matrix} \right)$$ 
 
For every different $\Gamma$ that is required it does four more inversions,

$$\chi^{\Gamma \alpha}_\beta = S_{\beta \gamma} (\Gamma \gamma_5)^\dagger_{\gamma \delta} \rho^{\alpha}_\delta$$

before calculating the correlator as,

$$C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\lambda = 0}^{3}\sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \chi^{\lambda}_{\beta, b}(\vec{x}, t; \tau) \chi^{\Gamma \lambda \dagger}_{\alpha, b}(\vec{x}, t; \tau)$$

where the $\lambda$ sum is over the $4$ spinor components.

We should be able to improve the signal and reduce the number of
inversions with two modifications. First, instead of having a different
noise vector for every spin component we reuse the same noise, i.e.

$$\rho^{\alpha}_\beta = \delta_{\alpha \beta} \rho$$ 

for fixed $\rho$.
Using less noise seems to be generally preferred.

Secondly there is no need to invert for every different $\Gamma$. Let,

$$\chi^{\Gamma \alpha}_\beta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha} \chi^{\gamma}_\beta$$

This is true because,

$$(\Gamma \gamma_5)^\dagger_{\gamma \alpha} \chi^{\gamma}_\beta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \delta} \rho^{\gamma}_\delta = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \delta} \delta_{\gamma \delta} \rho = (\Gamma \gamma_5)^\dagger_{\gamma \alpha}
S_{\beta \gamma} \rho = S_{\beta \gamma} (\Gamma \gamma_5)^\dagger_{\gamma \alpha} \rho$$

then the correlation function is

$$C(t - \tau, 0) = -\frac{1}{K} \sum_{k = 0}^{K} \sum_{\lambda = 0}^{3}\sum_{\vec{x} } (\gamma_5 \Gamma)_{\alpha \beta} \chi^{\lambda}_{\beta, b}(\vec{x}, t; \tau) \chi^{\Gamma \lambda \dagger}_{\alpha, b}(\vec{x}, t; \tau)$$

as before. By using the spin_matrix object in HiRep to construct the
objects $\chi^{\lambda}_{\beta, b}$ the correlators can be calculated
with only $4 N_F$ inversions.

### Disconnected

The Disconnected contributions occur when we have fermion species of the
same type in the hadron interpolator $O_M$:

$$O_M(x) = \bar{ \psi } (x) \Gamma \psi(x)$$ 

The same manipulations that lead to equation (5) give, 

$$\begin{aligned}
\langle O_M(x) \bar{O}_M'(y) \rangle &=& \langle \bar{\psi}(x) \Gamma \psi(x) \bar{\psi}(y) \bar{ \Gamma }' \psi(y) \rangle \\ 
=& \Gamma_{\alpha \beta} \bar{ \Gamma }'_{\gamma \delta} \langle \bar{\psi}_\alpha(x)  \psi_\beta(y) \bar{\psi}_\gamma(y)  \psi_\delta(x) \rangle \end{aligned}$$

There are two allowed Wick contractions,

$$\langle O_M(x) \bar{O}_M'(y) \rangle = -\text{Tr}\left[ \Gamma S (x,y) \bar{ \Gamma }' S (y,x) \right] + \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]$$

the connected and disconnected contributions. Fourier transforming the
first term gives us the same result as before. For the disconnected part, 

$$\begin{aligned}
D(t - \tau, \vec{p}) = \sum_{\vec{x} \vec{y}} e^{-i \vec{p} (\vec{x} - \vec{y} )} \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right]\end{aligned}$$

again $x = (\vec{x}, t)$ and $y = (\vec{y}, \tau)$, the zero momentum
correlator is 

$$\begin{aligned}
D(t - \tau, 0) &=& \sum_{\vec{x} \vec{y}} \text{Tr}\left[ \Gamma S (x,x) \right] \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right] \\
 =& \sum_{\vec{x}} \text{Tr}\left[ \Gamma S (x,x) \right] \sum_{\vec{y}} \text{Tr} \left[ \bar{ \Gamma }' S (y,y) \right] \end{aligned}$$
 
This means we have to evaluate objects like, 

$$\begin{aligned}
d(t) = \sum_{\vec{x}} \text{Tr}\left[ \Gamma S (x,x) \right] = \sum_{\vec{x}} \Gamma_{\alpha \beta} S_{\beta \alpha}(x,x)\end{aligned}$$

Using

$$\phi^{(k)}_{\beta, b}(\vec{x}, t; \tau) = \sum_{\vec{y}} S_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y})\rangle$$

which implies 

$$\begin{aligned}
\frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \text{Tr} \left[ \langle \eta^{(k)}(\vec{x}) | \Gamma | \phi^{(k)} (\vec{x}, t; \tau) \rangle \right] =& \frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \Gamma_{\beta \gamma} | \phi^{(k)}_{\gamma, b}(\vec{x}, t; \tau) \rangle \langle \eta^{(k)}_{\beta b }(\vec{x}) | \\
 =& \frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \Gamma_{\beta \gamma} \sum_{\vec{y}} S_{\gamma \alpha,b c} (\vec{x}, t; \vec{y}, \tau) 
| \eta^{(k)}_{\alpha c }(\vec{y})\rangle \langle \eta^{(k)}_{\beta b }(\vec{x}) |\end{aligned}$$

Using the limit $K \rightarrow \infty$ this becomes, 

$$\begin{aligned}\sum_{\vec{x}} \Gamma_{\beta \gamma} \sum_{\vec{y}} S_{\gamma \alpha,b c} (\vec{x}, t; \vec{y}, \tau) \delta_{bc} \delta_{\alpha \beta} \delta(\vec{y}, \vec{x}) =& 
\sum_{\vec{x}} \Gamma_{\beta \gamma} S_{\gamma \beta,b b} (\vec{x}, t; \vec{x}, \tau) \\
d(t, \tau) =& \text{Tr} \left[ \Gamma S(\vec{x}, t; \vec{x}, \tau) \right]\end{aligned}$$

We want only cases where $t = \tau$ so we need either (a) four noise
vectors on every timeslice or (b) noise vectors that are nonzero on all
timeslices. In case (b) we would evaluate, 

$$\begin{aligned}\frac{1}{K} \sum_{k}^K \sum_{\vec{x}} \text{Tr} \left[ \langle \eta^{(k)}(\vec{x}, t) | \Gamma | \phi^{(k)} (\vec{x}, t) \rangle \right] \end{aligned}$$

with

$$\phi^{(k)}_{\beta, b}(\vec{x}, t) = \sum_{\vec{y}} S_{\beta\gamma,b c} (\vec{x}, t; \vec{y}, \tau) | \eta^{(k)}_{\gamma c }(\vec{y}, \tau)\rangle$$

### Cancelling Backwards Propagation

The two point function evaluated in the centre of the lattice is
(including the backward propagating part to give the extra factor of
$2$),

$$C(T/2, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   2 e^{- E_\pi(\vec{p}) (T/2) }$$

therefore

$$\frac{1}{2} C(T/2, \vec{p}) e^{ - E_\pi(\vec{p}) (T/2 - t)  } = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) (T - t) }$$

then

$$C_{\rightarrow}(t, \vec{p}) = C(t, \vec{p}) - \frac{1}{2} C(T/2, \vec{p}) e^{ - E_\pi(\vec{p}) (T/2 - t)  }$$

is the forward propagating part only. $E_\pi(\vec{p})$ is obtained by
fitting the zero momentum correlator and using
$E(\vec{p}) = \sqrt{m_\pi^2 + \vec{p}^2}$. The factor $C(T/2, \vec{p})$
can be obtained also from the zero momentum correlator, by fitting to
obtain $|Z_\pi|^2$ and using the fact that this is momentum independent.
Since $0$ momentum results are used this might not be too noisy.

Alternatively the Wilson action is invariant under 

$$\begin{aligned}\psi(x) \rightarrow {\cal P_\mu}[ \psi(x) ] = \gamma_\mu \psi( P_\mu[x]) \\
\bar{ \psi } (x) \rightarrow {\cal P_\mu}[ \bar{ \psi } (x) ] = \bar{ \psi }( P_\mu[x])  \gamma_\mu \end{aligned}$$

where $P_\mu[x]$ reverses the sign of all the components of $x$ except
the $\mu$ one. Time reversal corresponds to ${\cal T} = {\cal P}_1{\cal P}_2{\cal P}_3$. 

$$\begin{aligned}\psi(x) \rightarrow {\cal T}[ \psi(x) ] = \gamma_0 \gamma_5 \psi( T[x]) \\
\bar{ \psi } (x) \rightarrow {\cal T}[ \bar{ \psi } (x) ] = \bar{ \psi }( T[x])  \gamma_5 \gamma_0 \end{aligned}$$

Using this the T symmetry of operators used to construct the correlators
can be calculated to calculate the sign on the backwards propagating
part.

$$\langle O_1(t) O_2(0) \rangle = C(t) = A \left( e^{ -Et} + \tau_1 \tau_2 e^{-E(T - t)}\right)$$

where $\tau_i = \pm 1$ is the ${\cal T}$ eigenvalue of $O_i$.

We mostly use correlators where $O_1 = O_2$ so $\tau_1 \tau_2 = 1$ then
the correlator is

$$C_{pp}(t, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } + e^{- E_\pi(\vec{p})(T - t) } + e^{- E_\pi(\vec{p})(2T - t) } + \ldots \right)$$

The subscript on $C_{pp}$ refers to the fact that both propagators used
periodic boundary conditions. We want to cancel the backwards
propagating part which can be done by solving the forward propagator
$S(0,x)$ using antiperiodic time bc's and the backward $S(x,0)$ with
periodic time bc's to give an extra minus sign,

$$C_{ap}(t, \vec{p}) = \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } - e^{- E_\pi(\vec{p})(T - t) } + e^{- E_\pi(\vec{p})(2T - t) } - \ldots \right)$$

so,

$$C_{ap}(t, \vec{p}) + C_{p}(t, \vec{p}) = \frac{ 2 |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) }   \left( e^{- E_\pi(\vec{p})t } + e^{- E_\pi(\vec{p})(2T - t) } + \ldots \right)$$

cancelling the subleading exponential. This method requires two
inversions and the calculation of

$$S_{A \pm P} (x,y) = S_A(x,y) \pm S_P(x,y)$$ 

Where the subscript refers to (A)ntiperiodic/(P)eriodic boundary conditions. Then,

$$C_{\pm}(t - \tau, 0) = -\sum_{\vec{x} \vec{y}} \text{Tr}\left[ \gamma_5 \Gamma S_{A \pm P} (x,y) \bar{ \Gamma }' \gamma_5 S_{A \pm P} (x,y) \right]$$

where $C_{+}(t - \tau, 0)$ gives the forward propagating part from $0$
to $T$ and $C_{-}(t - \tau, 0)$ gives the backwards propagating part
from $2T$ to $T$.

### Form Factors and Sequential Sources

The electromagnetic form factor of a 'pion' requires the evaluation of
the matrix element

$$\langle \pi(p_f) | V_\mu | \pi(p_i) \rangle = (p_i + p_f)_\mu f(q^2)$$

where $q^2 = (p_i - p_f)^2$ and

$$V_\mu = q_u \bar{u} \gamma_\mu u + q_d \bar{d} \gamma_\mu d$$ 

is the electromagnetic current and $q_i$ is the charge of the fermion $i$. This
is the local (not conserved) current, so there will be a factor $Z_V$
for renormalization. The matrix elements required look like:

$$\begin{aligned}
C_3(t_f,t,t_i,\vec{p_i}, \vec{p_f}) = Z_V \sum_{\vec{x} \vec{y} \vec{z}} e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} (\vec{y} - \vec{z)}}\langle 0 | \bar{u} \gamma_5 d ( \vec{x}, t_f)  V_0(\vec{y}, t) \bar{d} \gamma_5 u (\vec{z}, t_i) | 0 \rangle\end{aligned}$$

We take the $\mu = 0$ component since this is statistically cleaner and
also nonzero independant of the momentum direction. The contractions
give three propagators eg. taking the $\bar{d} \gamma_\mu d$ part of
$V_\mu$, 

$$\begin{aligned}
Z_V \sum_{\vec{x} \vec{y} \vec{z}} e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{-i\vec{p_i} (\vec{y} - \vec{z} ) } Tr \left[ S_u(\vec{z},t_i;\vec{x},t_f) \gamma_5 S_d(\vec{x},t_f;\vec{y},t)  \gamma_0 S_d(\vec{y},t;\vec{z},t_i) \gamma_5  \right]\end{aligned}$$

There are also disconnected contributions from contracting the two
fermions in the current together but we ignore those. The usual
sequential source trick consists of solving

$$S_u(\vec{x}, t; \vec{z}, t_i) = \sum_{\vec{y},\tau} D_u ( \vec{x}, t; \vec{y},\tau)  \delta(\vec{y}, \tau; \vec{z}, t_i)$$

to get the point-to-all propagator (for a specific $\vec{z}$ and $t_i$
as well as dropping the sum over $\vec{z}$ and using translational
invarience). Then taking a single timeslice of the propagator
$S_u(\vec{x}, t_f; \vec{z}, t_i)$ and solving, 

$$\begin{aligned}D_d ( \vec{x},t_f; \vec{y}, t ) G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i) =& e^{i\vec{p_f} \vec{x}} \gamma_5 S_u(\vec{x}, t_f; \vec{z}, t_i) \\
G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i) =& \sum_{ \vec{x} } e^{i\vec{p_f} \vec{x}} S_d(\vec{y},t; \vec{x}, t_f ) \gamma_5 S_u(\vec{x}, t_f; \vec{z}, t_i)\end{aligned}$$

to get the all-to-all-to-point contribution. Then 

$$\begin{aligned}\gamma_5 \left[ G_{du}(\vec{y}, t; \vec{p_f}; t_f; \vec{z}, t_i)  \right]^\dagger \gamma_5 &=& \sum_{ \vec{x} } e^{-i\vec{p_f} \vec{x}} \gamma_5 S_u^\dagger (\vec{x}, t_f;\vec{y},t ) \gamma_5 \gamma_5 \gamma_5 \gamma_5 \gamma_5 S_d^\dagger (\vec{z}, t_i;\vec{x}, t_f) \gamma_5 \\
=& \sum_{ \vec{x} } e^{-i\vec{p_f} \vec{x}} S_u (\vec{z}, t_i;\vec{x}, t_f) \gamma_5 S_d (\vec{x}, t_f;\vec{y},t ) \\
=& G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t )\end{aligned}$$

$$\begin{aligned}
C_3(t_f,t,t_i, \vec{p_i}, \vec{p_f}) = Z_V \sum_{\vec{y} } e^{-i(\vec{p_i} - \vec{p_f})\vec{y}} Tr \left[ G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t )  \gamma_0 S_d(\vec{y},t;\vec{z},t_i) \gamma_5 \right]\end{aligned}$$(eq:3ptfn)

This $Z_V$ factor is unknown. We show how to calculate it later, or
cancel it, but an alternative is to use the conserved vector current in
place of the local current:

$$V_\mu = \frac{1}{2} \left[ \bar{\psi}(x + \mu)(1 + \gamma_\mu)U_\mu^\dagger(x) \psi(x) - \bar{\psi}(x)(1 - \gamma_\mu)U_\mu^\dagger(x) \psi(x + \mu) \right]$$

The trace in {eq}`eq:3ptfn` becomes 

$$\begin{aligned}
Tr[ S_d(\vec{y},t+1;\vec{z},t_i) \gamma_5 (1 + \gamma_0)U_0^\dagger(\vec{y},t) G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t )   - \\ \nonumber
S_d(\vec{y},t;\vec{z},t_i) \gamma_5 (1 - \gamma_0)U_0(\vec{y},t) G_{ud}(\vec{z}, t_i; t_f; \vec{p_f}; \vec{y}, t+1 ) ]\end{aligned}$$

If we use this then all the following formulas are the same except $Z_V \rightarrow 1$.

There is an alternative that doesn't require the sequential source
trick. Using the properties of our noise sources,

$$S_d(\vec{x},t_f;\vec{y},t) \approx \frac{1}{K} \sum_{i=0}^K | \psi^{(i)}(\vec{x},t_f) \rangle \langle \eta^{(i)}(\vec{y},t) |$$

the three point correlation function becomes, 

$$\begin{aligned}
Z_V \frac{1}{K} \sum_{i=0}^K \sum_{\vec{x} \vec{y} } e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} \vec{y}} Tr \left[  \langle \eta^{(i)}(\vec{y},t) | \gamma_0 S_d(\vec{y},t;\vec{0},0) \gamma_5  S_u(\vec{0},0;\vec{x},t_f) \gamma_5 | \psi^{(i)}(\vec{x},t_f) \rangle \right] \\
Z_V \frac{1}{K} \sum_{i=0}^K \sum_{\vec{x} \vec{y} } e^{-i\vec{p_f}(\vec{x} - \vec{y}) } e^{i\vec{p_i} \vec{y}} Tr \left[  \langle \eta^{(i)}(\vec{y},t) | \gamma_0 S_d(\vec{y},t;\vec{0},0) S_u^{\dagger}(\vec{x},t_f;\vec{0},0) | \psi^{(i)}(\vec{x},t_f) \rangle \right]\end{aligned}$$

Using this method we can inject arbitrary momentum at the source without
the need for extra inversions.

#### Two Point Function

A complete set of hadrons is given by,

$$\sum_n \frac{ | n \rangle \langle n |}{ 2 E_n V}$$ 

the first term is the pion. The two point function (from point sources) is,

$$\begin{aligned}
C(t, \vec{p}) &=& \sum_{\vec{x}} e^{-i \vec{p} \vec{x}} \langle O_\pi (\vec{x}, t) O^\dagger_\pi (\vec{0}, 0)\rangle \\
&=& \sum_{\vec{x}} \sum_n e^{i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, t)  | n \rangle \langle n | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_n } \end{aligned}$$

Use $\sum_{ \vec{y} } e^{-i\vec{p} \vec{y}} O^\dagger_n (\vec{y}, 0) | 0 \rangle  = | n(\vec{p}) \rangle$, the time evolution operator $e^{-Ht}$ and also the fact that the lightest meson dominates the sum to get, 

$$\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} \vec{z}} \sum_{\vec{p'} } e^{-i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, 0)  O^\dagger_\pi (\vec{y}, 0) e^{i \vec{p'} \vec{y} }   | 0 \rangle \langle 0 | e^{i \vec{p'} \vec{z} } O_\pi (\vec{z}, 0) O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p'}) } e^{- E_\pi(\vec{p'}) t }\end{aligned}$$

The sum over $\vec{p'}$ gives a delta function leaving,

$$\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} } e^{-i \vec{p} \vec{x}}  
\frac{ \langle 0 | O_\pi (\vec{x}, 0)  O^\dagger_\pi (\vec{y}, 0) | 0 \rangle \langle 0 | O_\pi (\vec{y}, 0) O^\dagger_\pi (\vec{0}, 0) | 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}$$

Translational invarience lets us write, 

$$\begin{aligned}
C(t, \vec{p}) 
&=& \sum_{\vec{x} \vec{y} } e^{-i \vec{p} ( \vec{x} - \vec{y} ) } e^{ -i \vec{p} \vec{y} }  
\frac{ \langle 0 | O_\pi (\vec{0}, 0)  O^\dagger_\pi (\vec{x}-\vec{y}, 0) | 0 \rangle \langle 0 | O_\pi (\vec{y}, 0) O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}$$

Now changing variables gives us two Fourier transforms,

$$\begin{aligned}
C(t, \vec{p}) =& \frac{ \langle 0 |  O_\pi (\vec{0}, 0) | \pi(p) \rangle \langle \pi(p) | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}$$

and finally using the time evolution operator we get, 

$$\begin{aligned} C(t, \vec{p}) =& \frac{ |Z_\pi|^2 }{ 2 E_\pi(\vec{p}) } e^{- E_\pi(\vec{p}) t }\end{aligned}$$

where 

$$Z_\pi = \langle \pi(p) | O^\dagger_\pi (\vec{0}, 0)| 0 \rangle$$

#### Three Point Function

In less detail we insert two complete sets of states into the correlator
( point sources so $(\vec{x_i}, t_i) = (\vec{0}, 0)$ ) 

$$\begin{aligned}\langle \pi(p_f) | V_\mu | \pi(p_i) \rangle =& \langle 0| O(\vec{p_f}, t_f) V_\mu(\vec{x}, t) O^\dagger(\vec{p_i}, t_i) |0\rangle \\
 =& \langle 0| O(\vec{0}, 0) | \pi(\vec{p_f}) \rangle \frac{e^{-(t_f - t) E_\pi(\vec{p_f}) }}{2 E_\pi(\vec{p_f}) } \langle \pi(\vec{p_f}) | V_\mu(\vec{0}, 0) | \pi(\vec{p_i}) \rangle \\ \times& \frac{e^{-(t - t_i) E_\pi(\vec{p_i}) }} {2 E_\pi(\vec{p_i}) } \langle \pi(\vec{p_i}) | O^\dagger(\vec{0}, 0) |0 \rangle \\ \nonumber
=&  \frac{ Z_{\pi, f}^\dagger Z_{\pi, i} }{4 E(\vec{p_f}) E(\vec{p_i}) } \langle \pi(\vec{p_f}) | V_\mu(\vec{0}, 0) | \pi(\vec{p_i}) \rangle  e^{-(t_f - t) E_\pi(\vec{p_f}) -(t-t_i) E_\pi(\vec{p_i}) }\end{aligned}$$

if $t < t_f$ we have the backwards contribution and the exponential
changes to

$$-e^{-(t - t_f) E_\pi(\vec{p_f}) -(T - t + t_i) E_\pi(\vec{p_i}) }$$

#### Correlator Ratios: $Z_V$

$Z_V$ can be obtained as follows: The ratio, 

$$\begin{aligned}\frac{ C_{\rightarrow}(t_f, \vec{0}) }{ C_3(t_f, t, \vec{p_i}, \vec{p_f}) } =&  
\frac{ \frac{ |Z_\pi( \vec{0} )|^2 }{ 2 m_\pi }   e^{- m_\pi t_f }  }{ \frac{ |Z_\pi( \vec{0} )|^2 }{4 m_\pi^2 } \langle \pi(\vec{0}) | V_\mu | \pi(\vec{0}) \rangle e^{-(t_f - t) m_\pi -t m_\pi } } \\
=& \frac{ 1  }{ \frac{ 1 }{2 m_\pi } \langle \pi(\vec{0}) | V_\mu | \pi(\vec{0}) \rangle } \\
=& \frac{ 1  }{ \frac{ 1 }{2 m_\pi } 2 m_pi f(0)/Z_V } = Z_V.\end{aligned}$$

Where we used that the renormalized form factor $f(0) = 1$

#### Correlator Ratios: $f(q)$

There are various ways to cancel the unwanted terms and get $f(q)$,

##### RBC-UKQCD Ratio

We examine the ratio,

$$2 m_\pi \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{3} (t, t_f, \vec{0}, \vec{0} )  C_{\rightarrow}(t, \vec{p}) }$$

Assuming $Z_\pi$ is momentum independant (this also works (probably) if
$Z_\pi = E(\vec{p}) f_\pi$ , which is the case for
$O = \bar{u} \gamma_0 \gamma_5 d$ type interpolators) the numerator is,

$$\frac{ Z_V |Z_\pi|^2 }{ 4 E(\vec{p}) E(\vec{0})} f(q^2) ( E(\vec{p}) + m_\pi ) \frac{|Z_\pi|^2}{2 E(\vec{0})} e^{ -E(\vec{p})t - E(\vec{0})(t_f - t) -E(\vec{0})t }$$

and the denominator is,

$$\frac{ Z_V |Z_\pi|^2 }{ 4 E(\vec{0}) E(\vec{0})} f(0) ( m_\pi + m_\pi ) \frac{|Z_\pi|^2}{2 E(\vec{p})} e^{ -E(\vec{0})t - E(\vec{0})(t_f - t) -E(\vec{p})t }$$

Cancelling leaves,

$$2 m_\pi \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{3} (t, t_f, \vec{0}, \vec{0} )  C_{\rightarrow}(t, \vec{p}) } = f(q^2) ( E(\vec{p}) + m_\pi )$$

note there is no $Z_V$ here.

##### Bonnet et. al. Ratio

$$\frac{2 Z_V m_\pi}{E(\vec{p}) + m_\pi} \frac{ C_{3} (t, t_f, \vec{p}, \vec{0} )  C_{\rightarrow}(t, \vec{0}) }{ C_{\rightarrow} (t, \vec{p} )  C_{\rightarrow}(t_f, \vec{0}) }$$

the numerator of the right term is,

$$\frac{ |Z_\pi|^2 }{ 4 E(\vec{p}) m_\pi} f_B(q^2) ( E(\vec{p}) + m_\pi ) \frac{|Z_\pi|^2}{2 m_\pi} e^{ -E(\vec{p})t - m_\pi(t_f - t) -m_\pi t }$$

the denominator of the right term is,

$$\frac{ |Z_\pi|^2 }{ 2 m_\pi } \frac{|Z_\pi|^2}{2 E(\vec{p})} e^{ -m_\pi t - m_\pi (t_f - t) -E(\vec{p})t }$$

Cancelling leaves,

$$\frac{ f_B(q^2) ( E(\vec{p}) + m_\pi ) }{ 2 E(\vec{p}) }$$ 

the kinematic factors are cancelled

$$\frac{2 Z_V m_\pi}{E(\vec{p}) + m_\pi} \frac{ f_B(q^2) ( E(\vec{p}) + m_\pi ) }{ 2 E(\vec{p}) } = Z_V f_B(q^2) = f(q^2)$$

you need to actually know $Z_V$ or use the conserved current.


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

# Mesonic Correlators of the Isotriplet

The two fermionic flavors are denoted by $u$ and $d$. We are interested in the mesonic correlators

$$\begin{aligned}
&& C_{\Gamma_1,\Gamma_2}(x-y) = \left<
\left( \bar{u} \Gamma_1 d \right)^\dagger(x)
\left( \bar{u} \Gamma_2 d \right)(y)
\right> = \nonumber \\
&& \quad = \left<
\left( \bar{d} \gamma_0 \Gamma_1^\dagger \gamma_0 u \right)(x)
\left( \bar{u} \Gamma_2 d \right)(y)
\right>\end{aligned}$$ 

where $\Gamma_i$ are the generic producs of the $\gamma$-matrices.

We can integrate out the fermionic fields explicitly. Here we use the definition $H(x,y) = G(x,y) \gamma_5$ for the hermitian Dirac operator, with $G(x,y)$ defined as its inverse.

$$\begin{aligned}
&& C_{\Gamma_1,\Gamma_2}(x-y) = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,y) \Gamma_2 G(y,x) \right]
\right> = \nonumber \\
&& \quad = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,y) \Gamma_2 \gamma_5 G(x,y)^\dagger \gamma_5 \right]
\right> = \nonumber \\
&& \quad = - \left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 H(x,y) \gamma_5 \Gamma_2 H(x,y)^\dagger \gamma_5 \right]
\right> \; .\end{aligned}$$

Since the $\gamma$-matrices commute, we can conclude that the matrix $\gamma_0 \Gamma^\dagger \gamma_0$ is equal to $\Gamma$ up to a sign

$$\label{gamma0_adj}\gamma_0 \Gamma^\dagger \gamma_0 = s(\Gamma) \Gamma \qquad \text{with } s(\Gamma) = \pm 1 \; .$$

In addition, a generic matrix $\Gamma$ has the following properties:

1.	Its matrix elements can be $0$, $\pm 1$, $\pm i$
2.	Its entries are either all real or all imaginary
3.	In each row and correspondingly each column, there is only one non-zero element

Consequently, we can write

$$\label{gamma_ab}
\Gamma_{\alpha\beta} = t_\alpha(\Gamma) \delta_{\sigma_\alpha(\Gamma), \beta}$$

where $\sigma(\Gamma)$ constitutes a permutation of four elements. Putting this together we find

$$\begin{aligned}
& C_{\Gamma_1,\Gamma_2}(x-y) = - s(\Gamma_1) \left< \mathrm{tr}
\left[ \gamma_5 \Gamma_1 H(x,y) \gamma_5 \Gamma_2 H(x,y)^\dagger \right]
\right> = \\
& \quad = - s(\Gamma_1) \sum_{\alpha\beta} t_\alpha(\gamma_5 \Gamma_1) t_\beta(\gamma_5 \Gamma_2) \times \left< \mathrm{tr}
\left[ H_{\sigma_\alpha(\gamma_5 \Gamma_1), \beta}(x,y) H_{\alpha, \sigma_\beta(\gamma_5 \Gamma_2)}(x,y)^\dagger \right]
\right> \; .\label{triplet_corr}\end{aligned}$$

## Implementation of the Point-To-All Propagator

In order to calculate mesonic masses we are interested in correlators satisfying $\Gamma_1=\Gamma_2$. Using translational invariance, we can set $y=0$. In this case the formula simplifies to

$$C_{\Gamma}(x) = - s(\Gamma) \sum_{\alpha\beta} t_\alpha(\gamma_5 \Gamma)t_\beta(\gamma_5 \Gamma) \times \left< \mathrm{tr}\left[ H_{\sigma_\alpha(\gamma_5 \Gamma), \beta}(x,0) H_{\alpha, \sigma_\beta(\gamma_5 \Gamma)}(x,0)^\dagger \right]\right> \; .$$(eq:triplet_point_to_all_corr)

This is implemented into HiRep in the following way

* The data of the point-like source $\xi^{(\alpha,a)}$ defined by

    $$\xi^{(\alpha,a)}_{\beta b}(x) = \delta_{\alpha,\beta} \delta_{a,b} \delta_{x,0} \; ,$$
	
    The function `quark_propagator` applies the inverse of the hermitian Dirac operator to the source
	
    $$\eta^{(\alpha,a)} = H \xi^{(\alpha,a)} \qquad \eta^{(\alpha,a)}_{\beta b}(x) = H_{\beta \alpha}^{b a}(x,0) \; .$$

-   The functions `void *_correlator(float *out, suNf_spinor **qp)`in `Observables/mesons.c` implement the formulae {eq}`eq:triplet_point_to_all_corr`, where `out` stands for the correlator and `qp` for the spinor array. The functions $s(\Gamma)$, $t_\alpha(\gamma_5 \Gamma)$ and $\sigma_\alpha(\gamma_5 \Gamma)$ where calculated using `Mathematica`, see file `mesons.nd` and implemented in the code using macros, defined as follows\
    
    `      `\
    `      _C1_ = `$\sigma_1(\gamma_5 \Gamma)$\
    `      _C2_ = `$\sigma_2(\gamma_5 \Gamma)$\
    `      _C3_ = `$\sigma_3(\gamma_5 \Gamma)$\
    `      _C4_ = `$\sigma_4(\gamma_5 \Gamma)$\
    `      `
    
    If $t_\alpha(\gamma_5 \Gamma)$ are read:
    
    `      `\
    `      _S0_ = `$-s(\Gamma)$\
    `      _S1_ = `$t_1(\gamma_5 \Gamma)$\
    `      _S2_ = `$t_2(\gamma_5 \Gamma)$\
    `      _S3_ = `$t_3(\gamma_5 \Gamma)$\
    `      _S4_ = `$t_4(\gamma_5 \Gamma)$\
    `      `
    
    whereas if $t_\alpha(\gamma_5 \Gamma)$ are imaginary:
    
    `      `\
    `      _S0_ = `$s(\Gamma)$\
    `      _S1_ = `$-i t_1(\gamma_5 \Gamma)$\
    `      _S2_ = `$-i t_2(\gamma_5 \Gamma)$\
    `      _S3_ = `$-i t_3(\gamma_5 \Gamma)$\
    `      _S4_ = `$-i t_4(\gamma_5 \Gamma)$\
    `      `

# Mesonic Correlators of the Isosinglet

We are now concerned with the genertic mesonic correlator given by

$$\begin{aligned}
C_{\Gamma_1,\Gamma_2}(x-y) &=& \left<
\left( \bar{u} \gamma_0 \Gamma_1^\dagger \gamma_0 u \right)(x)
\left( \bar{u} \Gamma_2 u \right)(y)
\right>\end{aligned}$$

considering only a single flavor $u$.

Integration of the fermionic fields now yiels one additional term, the hairpin diagram

$$\begin{aligned}
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
\; .\end{aligned}$$

All other contributions are identical to the contributions to the isotriplet correlator. We can, therefore, focus on the contribution through the hairpin diagram. Using the formulae {eq}`eq:gamma0_adj` we can write

$$\left< \mathrm{tr}
\left[ \gamma_0 \Gamma_1^\dagger \gamma_0 G(x,x) \right]
\mathrm{tr}\left[ \Gamma_2 G(y,y) \right] \right> = \quad = s(\Gamma_1) \sum_{\alpha \beta} t_\alpha(\Gamma_1) t_\beta(\Gamma_2) \times \left<
\mathrm{tr}G_{\sigma_\alpha(\gamma_5 \Gamma_1),\alpha}(x,x) \;
\mathrm{tr}G_{\sigma_\beta(\gamma_5 \Gamma_2),\beta}(y,y)
\right>
\; .$$(eq:hairpin)

## All-to-all Propagator

It is clear from {eq}`eq:hairpin`, from the fact that we are employing point source, that one must compute the entire inverse matrix of the Dirac operator. The alternative is to use a statistic estimate for $H$ followed by variance reduction procedures.

Suppose there are $N_s$ available random fermion sources $\xi^{(i)}$ such that the only non-zero correlators are 

$$\left< \xi^{(i)}_{\alpha a}(x)^\dagger \xi^{(j)}_{\beta b}(y) \right> = \delta_{\alpha,\beta} \delta_{a,b} \delta_{x,y} \delta_{i,j} \; .$$

Current literature proposes mainly either Gaussian noise or $Z_2$ noise. In the following, we will choose $Z_2$ noise, following [@Foster:1998vw]. Each component of the spinor is randomly chosen from the values $\pm 1/\sqrt{2}$.

Then the matrix $H$ can be estimated as follows:

$$H_{\alpha \beta}^{a b}(x,y) \simeq \sum_{i=1}^{N_s} \eta^{(i)}_{\alpha a}(x) \xi^{(i)}_{\beta b}(y)^\dagger \eta^{(i)} \equiv H \xi^{(i)}$$(eq:naive_noisy_estimate)

Stochastic estimation can then be used to calculate the relevant tracks for correlators

$$\mathrm{tr}\left[ \Gamma_1 G(x,y) \Gamma_2 G(y,x) \right] = \sum_{ij} \xi^{(i)}(x)^\dagger \gamma_5 \Gamma_1 \eta^{(j)}(x) \times \xi^{(j)}(y)^\dagger \gamma_5 \Gamma_2 \eta^{(i)}(y)$$

$$\mathrm{tr}\left[ \Gamma G(x,x) \right] = \sum_i \xi^{(i)}(x)^\dagger \gamma_5 \Gamma \eta^{(i)}(x)\, .$$

## Variance reduction

The noise obtained from stochastic estimation of the matrix $G$ in the formula {eq}`naive_noisy_estimate` can be reduced using the trick from [@McNeile:2000xx] for Wilson fermions. Here, the Dirac operator has the form $D = 1 - K$. As a result, for the matrix $G$ the following formula applies 

$$\begin{aligned}
&& G = D^{-1} = \left( 1 - K \right)^{-1} = \\
&& \quad = 1 + K + \dots + K^m + K^{n_1} G K^{n_2}\end{aligned}$$ 

with $n_1 + n_2 = n = m+1$. In particular, for the evaluation of the hairpin diagram

$$\begin{aligned}
\mathrm{tr}\left[ \Gamma G(x,x) \right] &=& \mathrm{tr}\Gamma + \mathrm{tr}\left[ \Gamma K^4 \right](x,x) + \mathrm{tr}\left[ \Gamma K^6 \right](x,x) + \dots + \\
&& + \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) + \mathrm{tr}\left[ \Gamma K^{n_1} G K^{n_2} \right](x,x)\end{aligned}$$

with $m=2k$. (TODO: fix this sentence) Here, we can use the fact that the matrix $K$ connects first neighboring sites as thus $K^p(x,x) \neq 0$ only when $p$ is even. Further, $r_0=1$ and consequently $K^2(x,x) = 0$ The first $k+1$ terms can be calculated explicitly and we can estimate the last term stochastically.

$$\begin{aligned}
\mathrm{tr}\left[ \Gamma G(x,x) \right] =& \mathrm{tr}\Gamma + \mathrm{tr}\left[ \Gamma K^4 \right](x,x) + \mathrm{tr}\left[ \Gamma K^6 \right](x,x) + \dots + \nonumber \\
& + \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) + \sum_{iy} \xi^{(i)}(y)^\dagger \gamma_5 K^{n_2}(y,x) \Gamma K^{n_1}(x,y) \eta^{(i)}(y) \nonumber \\
\end{aligned}$$(eq:hairpin_with_variance_reduction)

[@McNeile:2000xx] use this trick only for the calculation of the hairpin diagram. It might be possible to generalize it to the the isotriplet part as well, as an alternative to the point-to-all propagator. 

## Time dilution

This is a trick introduced in [@Foley:2005ac] for noise reduction in the computation of null-moment propagators. Whenever stochastic estimation of the $H$ matrix is required, such as in {eq}`eq:naive_noisy_estimate`, it is possible to replace each stochastic source $\xi^{(i)}$ with a set of sources each with support on a different time slice.

$$\begin{aligned}
 \label{time_dilution}
&& \xi^{(i)} \rightarrow \xi^{(i,\tau)}(\mathbf{x},t) = \xi^{(i)}(\mathbf{x},t) \delta_{t,\tau}\end{aligned}$$

Stochastic estimation is now obtained similarly to the naive case:

$$H_{\alpha \beta}^{a b}(x,y) \simeq \sum_{i=1}^{N_s} \sum_{\tau=1}^{N_t} \eta^{(i,\tau)}_{\alpha a}(x) \xi^{(i,\tau)}_{\beta b}(y)^\dagger  \eta^{(i,\tau)} \equiv H \xi^{(i,\tau)}$$(eq:diluted_noisy_estimate)

## Implementation Scheme 

TODO: Add this to function reference instead if this is still implemented this way

The following functions will be implemented


-   Calculation of the exact terms of the formula {eq}`hairpin_with_variance_reduction`
-   
    ` `\
    ` void GAMMA_variance_reduction_exact_terms(float *out, int k)`
    
    Here, ` out ` is a real vector with its components equal to the volume and $k$ corresponds to the index in the $\gamma$-matrix.
    
    This function evaluates
    
    $$h_k(x) =& \mathrm{tr}\left[ \Gamma K^{2k} \right](x,x) = \left( \frac{\kappa}{2} \right)^{2k} \sum_{\mathcal{C}_x} \mathrm{tr}\left( \Gamma \tilde{\gamma}(\mathcal{C}_x) \right) \mathcal{W}(\mathcal{C}_x)$$
    
    where $\mathcal{C}_x = (x,\hat{\mu}_1,\hat{\mu}_2,\dots,\hat{\mu}_{2k})$ is the generic closed path of $x$ of length $2k$ obtained by moving from $x$ in the directions $\hat{\mu}_i$. $\mathcal{W}(\mathcal{C}_x)$ is the trace of the parallel transport through $\mathcal{C}_x$ in the corresponding fermionic representation. $\tilde{\gamma}(\mathcal{C}_x)$. is the matrix defined as
    
    $$\tilde{\gamma}(\mathcal{C}_x) = (1-\gamma_{\hat{\mu}_1})(1-\gamma_{\hat{\mu}_2})\cdots(1-\gamma_{\hat{\mu}_{2k}})$$
    
    having defined $\gamma_{-\hat{\mu}_i} = - \gamma_{\hat{\mu}_i}$.

	It should be noted that since $(1-\gamma_{\hat{\mu}_i})(1-\gamma_{-\hat{\mu}_i}) = 0$, one can exclude the paths in which a pair of subsequences $(\dots, \hat{\mu}_i,-\hat{\mu}_i, \dots)$ appears from the sum. In addition, the matrix $\tilde{\gamma}(\mathcal{C}_x)$ does not depend on $x$. There is is convenient to compute the list of paths and the matrix $\tilde{\gamma}(\mathcal{C}_x)$ only once.

-   It is convenient to have a function that calculates the traces of the parallel transports:
-   
    ` void tr_r_pexp(complex *out, int *path, int length)`
    
    Here again ` out ` is the complex vector with number of components according to the volume and $\phi(x)$\ is the number of directions (`length`) of which the path is composed.
    
    This returns
    
    $$\phi(x) = \mathrm{tr}\mathbf{R} \left[ U_{x, \hat{\mu}_1} U_{x+\hat{\mu}_1, \hat{\mu}_2} \cdots \right]$$


-   Calculation of the time-diluted estimators
    
    ` void time_diluted_stochastic_estimate`\
    `      (suNf_spinor *csi, suNf_spinor **eta,`\
    `      int nm, float *mass, double ac)`\
    
    With the parameters:
    ` csi `  is the spinor $\xi$\
    ` path ` is the list of $N_t \times \textrm{nm}$ spinors along the $\eta^{(\tau,m)}$
    ` `\
    
    This function generates the spinor $\xi$ with $Z_2$ noise. Here the spinors are defined as 
    
    $$\xi^{(\tau)}(\mathbf{x},t) = \xi(\mathbf{x},t) \delta_{t,\tau}$$ 
    
   and then returned as 
    
    $$\eta^{(\tau,m)} \equiv H_m \xi^{(\tau)}$$
    
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
	    

The functions above implement the formula {eq}`eq:hairpin_with_variance_reduction`, summing exact terms and the statistical term, generated with `nrs` to dilute, for a total of $\times N_t$ matrix invertions for each mass value and returns the result as `out`.

TODO: Add this to bibliography

J. Foley, K. Jimmy Juge, A. O'Cais, M. Peardon, S. M. Ryan and
J. I. Skullerud, Comput. Phys. Commun. **172**, 145 (2005)
\[arXiv:hep-lat/0505023\].

C. McNeile and C. Michael \[UKQCD Collaboration\], Phys. Rev. D **63**,
114503 (2001) \[arXiv:hep-lat/0010019\].

M. Foster and C. Michael \[UKQCD Collaboration\], Phys. Rev. D **59**,
074503 (1999) \[arXiv:hep-lat/9810021\].

A. Hart, C. McNeile, C. Michael and J. Pickavance \[UKQCD
Collaboration\], Phys. Rev. D **74**, 114504 (2006)
\[arXiv:hep-lat/0608026\].
:::


