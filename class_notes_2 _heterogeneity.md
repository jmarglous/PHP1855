## Lecture 10/14: Heterogeneity
To date we have assumed that everyone has an equal chance of getting infected -- now we expand our modeling framework. 
	Heterogeneity developed to model STIs. 
		- More likely to have asymptomatic/latent periods and not produce immunity.
		- Gonorrhea produced a key model for the field: non-immunizing bacteiral infection with or without immunity

#### Model without heterogeneity

S -> I (lambda S)
I <- S (gamma I)

$$
\begin{eqnarray}
\frac{dS}{dt} &=& -\lambda S +\gamma I \\
\frac{dI}{dt} &=& \lambda S -\gamma I)\\
\end{eqnarray}
$$
$\lambda$ is the force of infection (what we called $\beta I$ before). We define it as $\lambda = c\beta_p I/N$. 
- $c$ is the number of partners per year (sometimes called partner change rate).
- $\beta_p$  is transmission probability per partnership.
- $I/N$ is the proportion of infected individuals in the population.

#### R0 Calculation
In a fully susceptible population where $S = N$:
$$
\frac{dI}{dt} = \frac{c\beta_p I}{N} S - \gamma I > 0
$$
$$
R_0 \equiv \frac{c\beta_p}{\gamma} > 1
$$
or, 
$$
R_0 = c\beta_pD 
$$
where $D = 1/\gamma$ is the duration of infectiousness. 

Reasonable values for these parameters are $c = 2, \beta_p = .75, D = 2/12$ (in units of years). These give $R_0 = .25$.
	So even with a high probability of infection the infection wouldn't invade -- so why do we see gonorrhea in the population?

#### With heterogeneity
Two groups: high and low risk, with the same recovery rate but varying force of infection
S_H -> I_H (lambda_H S_H)
I_H -> S_H (gamma I_H)

S_L -> I_L (lambda_L S_L)
I_L -> S_L (gamma I_L)

$H = S_H + I_H + S_L + I_L$. 
Now, $\lambda_j(t)  = c_j \beta_p p (t)$ where $p(t)$ is the probability of meeting an infected partner
	$c(j)$ is the partner change rate
	 $p(t)$ is the population prevalence of infectious individuals in $H$.
		 $p(t) = \frac{g_H I_H}{M_H} + \frac{g_L I_L}{M_L}$
	and $g_H$ and $g_L$ are the probabilities of selecting a partner from group $H$ and$L$, respectively:
				$g_H + g_L = 1$.
Let's set some parameter values: assume that $N_H$ and $N_L$ are fixed as $\frac{N_H}{N}=.2$, $\frac{N_L}{N} = .98$. (so only 2% of population is in high-risk group)
	$$g_H = \frac{\frac{c_h N_H}{N}}{\frac{c_HN_H}{N} + \frac{c_LN_L}{N}}$$ (i..e partnerships in H per year over partnerships in whole population per year)
$$g_L = \frac{\frac{c_L N_L}{N}}{\frac{c_LN_L}{N} + \frac{c_HN_H}{N}}$$
	There are two aspects being woven in: 
		1. Ratio of individuals in high-risk and low-risk groups.
		2. With assumption of random choice of partners, there is a higher chance of drawing a high-risk individual.

Assume $c_{mean} = 2$. $c_L = 1.4$, then $c_H = \frac{c_{mean} - \frac{c_LN_L}{N}}{\frac{N_H}{N}} = 31.4$
	The cut-off between high- and low-risk groups is a bit arbitrary so this is potentially not super realistic but was a useful start for the field. 
		$g_h = .314$, $g_l = .686$ ($g_h$ is much higher than the population proportion of high-risk individuals, $.02$). 

### Estimating R0
Estimating $R0$ relies on a new method. This is the gold-standard method for estimating R0. 

We would also need to use this next-generation approach for SEIR model as well because there are two infected classes. This approach will work for most compartmental models of any complexity. This is a matrix-based linear algebra approach. 