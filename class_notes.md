# SIR model

The SIR model is the most simple and widely used model of infectious disease. It consists of three compartments: susceptible (S), infected (I), and recovered (R), with simple flow through these compartments. 

![[./figures/lecture1/SIR_dag.png]]
This model can be formulated with discrete or continuous time. 
## Discrete time SIR model
$$
\begin{eqnarray}
S_{t+1} &=& S_t - \lambda_tS_t\\
I_{t+1} &=& I_t + \lambda_tS_t - rI_t\\
R_{t+1} &=& R_t + rI_t\\
N &=& S_{t+1} + I_{t+1} + R_{t+1}
\end{eqnarray}
$$
$\lambda_t$ is the proportion of susceptible individuals infected in the interval $[t, t+1]$. $r$ is the proportion of infected individuals who recover in $[t, t+1]$. This can also be considered the **risk of infection** or **course of infection**. Note that the population is unchanging (the terms with $\lambda_t$ and $r$ cancel with each other, so $S_t + I_t + R_t = S_{t+1} + I_{t+1} + R_{t+1}$). *(In fact, a common bug when writing simulations is forgetting a term so that the population is allowed to change. This will cause really weird behavior.)* 

We can define $\lambda_t = \beta I_t$ , reflecting that the risk of infection depends on the number of infected individuals. $\beta$ is the per-capita rate at which two individuals come into **effective contact**, i.e. a contact that transmits the infection. This itself can be broken down into two components: the overall contact rate and the effectivity of that contact. Mathematically, $\beta = cq$, where $c$ is the contact rate in units of $\frac{[contacts]}{[time]}$ and $q$ is the probability of transmission given contact, in units of $\frac{[infections]}{[contact]}$. So $\beta \sim \frac{[infections]}{time}$. 
- Contacts $c$ have high heterogeneity within populations (think of the different contacts you have during the day, with different durations, closeness, etc.), and a current move in the infectious disease modeling field is to incorporate heterogeneity into models and see how they influence predictions. 
- Effectivity $q$ also has heterogeneity by pathogen (infectivity, pathogenicity, virulence) and host (age, health status, etc.) characteristics.

We might define a new parameter 
$$ R_0 = cqD = \beta D$$
	, where $D \sim \frac{[time]}{[infection]}$ is the duration of infectiveness. This new parameter is unitless. 

Consider a worked example with the following values, to get a feel for the numbers:
$$
\begin{matrix}
\beta = .02\\
r = 2\\
S = 200, I = 1, R = 99, N = 300
\end{matrix}
$$
The results will be: 
$$
\begin{array}{|T|N|S|I|R|} \hline \textbf{Time} & \textbf{N} & \textbf{S} & \textbf{I} & \textbf{R} \\ 

\hline t & 300 & 200 & \hline 1 & 99 \\ 
t+1 & 300 & 196 & 3 & 101 \\ 
t + 2 & 300 & 184.24 & 8.76 & 107 \\\hline 
\end{array}
$$
A major issue with using the discrete time model is defining what $t$ is -- changing the size of the time step can change the results of the model. 

## Continuous time SIR model
We can do a bit of calculus to convert the discrete-time SIR model to a continuous time version. All we do is compress the time step $t$ to be infinitely small. 

$$
	\lim_{\Delta t \rightarrow 0} \frac{I(t+\Delta t)-I(t)}{\Delta t} = \frac{\partial I}{\partial t}
$$
![[derivative.png]]


From the discrete time model,
$$
S_{t+1}-S_t = S_t - \beta I_tS_t - S_t = -\beta I_tS_t
$$
Change variables from the discrete time model: $t+1 \rightarrow t+ \delta$ so $\beta \rightarrow \delta \beta$. 
$$
dS = -\beta I_t S_t dt \rightarrow \frac{dS}{dt} = -\beta I_tS_t
$$
We can do  a similar calculation for $I$ and $R$.
$$
\begin{matrix}
dI = (\beta SI -rI)dt\\
dR = rI
\end{matrix}
$$
To differentiate from the discrete time model, we replace $r$ with $\gamma$ and we get the classic SIR model. 
$$
\begin{eqnarray}
\frac{dS}{dt} &=& -\beta SI\\
\frac{dI}{dt} &=& I(\beta S -\gamma )\\
\frac{dR}{dt} &=& \gamma I
\end{eqnarray}
$$
Some things to notice with this model:
1. There is no way to leave the recovered compartment (no reinfections).
2. The interaction between variables in the $\beta S I$ terms make this a *nonlinear* model. 
3. In the continuous model we often set $N =1$ and define $S, I, R \in [0, 1]$ as proportions. 

### Defining $R_0$
To analyze this model further, let's investigate what is needed for $\frac{dI}{dt} > 0$ (i.e. for the infection to spread.) We need  $I>0$ (trivial; of course there need to be infections to seed other infections) and $\beta S - \gamma > 0$.
- Rewriting this second condition as $S > \frac{\gamma}{\beta}$, we can recognize that this defines a *critical proportion* of susceptibles required for the infection to spread. ==
	- In a fully susceptible population (novel disease), $S = N =1$ and this critical proportion becomes.
$$R_0 \equiv \frac{\beta}{\gamma} > 1$$
		- So $R_0$ has a clear interpretation: if it is greater than 1, the disease will spread. We can also understand it as the average number of secondary infections arising from a single infection in a fully susceptible population (the multiplicative growth rate) or the **basic reproduction number**. 

$\gamma$ has units of $\frac{1}{[time]}$ (because $\frac{dR}{dt}$ is defined in $\frac{[infections]}{[time]}$ and $I$ in $[infections]$), so we can define a unit of $[time]$ as 
$$
	D \equiv \frac{1}{\gamma}
$$
, so $D$ is the duration of infectiousness. Substituting, $R_0 = \beta D$, matching the definition of $R_0$ for the discrete time model (with $N =1$). 

### Defining $R_{e}$
In a endemic/non-pandemic context where less than the whole population is susceptible, the **effective reproduction number** rather than the basic reproduction number is the threshold for invasion:
$$
\begin{eqnarray}
\frac{dI}{dt} > 0 \rightarrow S &>& \frac{\gamma}{\beta}\\
	S &>& \frac{1}{R_0}\\
	R_{e} &\equiv& R_0S > 1
\end{eqnarray}
$$
We can use $R_e$ to plot the course of an epidemic:
![[./figures/lecture1/pandemic_timecourse.png]]
We can re-write this equation in terms of $S$, because the number of susceptibles is the main policy lever we can control in a pandemic (i.e. via vaccination, etc.). To prevent invasion, we need $S < \frac{1}{R_0}$.
For an infection to spread, 
$$
\begin{eqnarray}
-\frac{1}{R_0} &<& -S\\
1-\frac{1}{R_0} &<& 1 -S \\
1-S &>& 1-\frac{1}{R_0} \equiv P_c
\end{eqnarray}
$$
$1-S$ is the proportion of the population with immunity. (There are two ways to get immunity: vaccination and infection.) So $P_c$ is the **critical proportion** or **herd immunity threshold** that represents the minimum immune population fraction required to prevent invasion. 

### Estimating $R_0$ in early pandemic
In the early pandemic $S = N = 1$. This removes the non-linearity from the SIR model, leaving 
$$
\frac{dI}{dt} = I(\beta -\gamma)
$$
We can easily solve this equation: 
$$
I(t) = I_0 e^{kt}
$$
, where $k = \beta - \gamma$. Taking a log to get a linear equation, 
$$
ln(I(t)) = ln(I_0) + kt
$$
So, if we have data on the number of infections per time, we can fit a linear model to estimate $k$. 

We can also rewrite $k$:
$$
\begin{eqnarray}
k &=& \beta - \gamma\\
\frac{k}{\gamma} &=& \frac{\beta}{\gamma} - 1\\
&=& R_0 -1 \\
R_0 &=& \frac{k}{\gamma} +1 = kD +1
\end{eqnarray}
$$
(remember that $D = \frac{1}{\gamma}$.) So if we have duration data, we can combine that with our linear model to estimate $R_0$ (and then we will know if the disease is spreading, what the critical proportion is, etc.)
![[./figures/lecture1/estimating_k.png]]