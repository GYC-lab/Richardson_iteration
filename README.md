**1. Introduction**

We are concerned with solving a given linear system
$$Au=b, \quad A \in {R}^{N \times N}, \quad b \in {R}^N$$ by an
iteration scheme of the form
$${u}_{k+1}= {u}_k+\alpha\left( {b}-A  {u}_k\right), \quad k=0,1,2, \ldots
    \label{Richardson}$$ which is called as **Richardson iteration**.
The basic idea of Richardson iteration is to start with a guess and
improve it. With different relaxation parameters $\alpha$, the iteration
may diverge or converge with different speed. We can also convert
([[\[Richardson\]](#Richardson){reference-type="ref"
reference="Richardson"}]{style="color: blue"}) into the following form
$${u}_{k+1}=(I-\alpha A ) {u}_k+\alpha {b}$$ where $I$ is the unit
matrix. The definations of the residual and error are listed below:
$${r}={b}-A {v}, \quad r \in {R}^N$$ $${e}={v}-{u}, \quad e \in {R}^N$$

As required by Homework 2, the linear system to be solved is
$$\left[\begin{array}{ll}
        6 & 3 \\
        3 & 4
        \end{array}\right]\left[\begin{array}{l}
        x_1 \\
        x_2
        \end{array}\right]=\left[\begin{array}{l}
        -3 \\
        -9
        \end{array}
    \right]$$

In $\S$ [2]{style="color: blue"}, the numerical results obtained by
Richardson iteration with a initial solution $x_0=[0.0,0.0]^T$ and
different relaxation parameters $\alpha$ are prresented. In addition, a
discussion on the effects of relaxation parameters and preconditioning
on the convergence behaviours is provided.

**2. Solutions**

::: center
2.1. *Numerical results of Richardson iteration*
:::

With the aid of MATLAB, the numerical results of different cases have
been calculated. The errors and residuals of the first 10 steps in the
form of 2-norm are listed in Table [[1](#errors){reference-type="ref"
reference="errors"}]{style="color: blue"} and Table
[[2](#residuals){reference-type="ref"
reference="residuals"}]{style="color: blue"} respectively. At the same
time, the convergence curve of the first 50 steps are shown in Fig.
[[1](#velocity){reference-type="ref"
reference="velocity"}]{style="color: blue"}.

::: {#errors}
  steps   $\alpha=0.06$   $\alpha=0.1$   $\alpha=0.2$   $\alpha=0.22$   $\alpha=0.24$   $\alpha=0.4$
  ------- --------------- -------------- -------------- --------------- --------------- --------------
  1       2.728369        2.469818       2.000000       1.948333        1.914158        2.280351
  2       2.402067        2.011219       1.264911       1.226250        1.281050        4.841487
  3       2.129703        1.641417       0.800000       0.795471        0.984463        10.955181
  4       1.892674        1.339762       0.505964       0.536039        0.850022        24.812191
  5       1.683337        1.093551       0.320000       0.376458        0.782130        56.197395
  6       1.497538        0.892586       0.202386       0.274704        0.738963        127.282101
  7       1.332359        0.728554       0.128000       0.206677        0.704987        288.282639
  8       1.185433        0.594666       0.080954       0.158897        0.674855        652.934540
  9       1.054718        0.485383       0.051200       0.123888        0.646762        1478.838663
  10      0.938420        0.396183       0.032382       0.097423        0.620084        3349.438050

  : $||e||_2$ of the first 10 steps
:::

[]{#errors label="errors"}

::: {#residuals}
  steps   $\alpha=0.06$   $\alpha=0.1$   $\alpha=0.2$   $\alpha=0.22$   $\alpha=0.24$   $\alpha=0.4$
  ------- --------------- -------------- -------------- --------------- --------------- --------------
  1       6.307139        4.743416       6.000000       6.958448        8.004998        17.492856
  2       4.826813        3.704727       3.794733       5.255294        7.286283        39.481641
  3       4.038605        3.016828       2.400000       4.050665        6.855425        89.418119
  4       3.515160        2.462126       1.517893       3.163205        6.530513        202.523954
  5       3.104361        2.009643       0.960000       2.489847        6.248187        458.698740
  6       2.755240        1.640326       0.607157       1.968976        5.987036        1038.911852
  7       2.449437        1.338880       0.384000       1.561250        5.739734        2353.042948
  8       2.178769        1.092831       0.242863       1.239844        5.503606        5329.433007
  9       1.938359        0.891999       0.153600       0.985453        5.277504        12070.691784
  10      1.724579        0.728074       0.097145       0.783639        5.060792        27339.043372

  : $||r||_2$ of the first 10 steps
:::

[]{#residuals label="residuals"}

![The convergence curves with different relaxation
parameters](../images/1_convergence_curve_for_every_alpha.png){#velocity
width="80%"}

::: center
2.2. *The condition of convergence*
:::

The constrant for the Richardson iteration to convergence is
$$0<\alpha<\frac{2}{\lambda_{max}}$$ and the optimal relaxation
parameter $\alpha_{opt}$ and the corresponding spectral radius
$\rho_{opt}$ for the iteration to convergence are
$$\alpha_{opt}=\frac{2}{\lambda_{max}+\lambda_{min}}$$
$$\rho_{opt}=\frac{\lambda_{max}-\lambda_{min}}{\lambda_{max}+\lambda_{min}}$$

The calculated results are $\alpha_{opt}=0.200000$,
$\rho_{opt}=0.632456$. And the range of $\rho$ for the Richardson
iteration to convergence is $(0.000000,0.245030)$. Since $\alpha_6$ is
out of the convergence range($\alpha_6=0.4>0.245030$), the case of
$\alpha_6=0.4$ will diverge. And for the given linear system, the
Richardson iteration converge with $\alpha = 0.06,0.1,0.2,0.22,0.24$. As
the current relaxation parameter is closer to the $\alpha_{opt}$, the
iteration converges faster. Hence, from Fig.
[[1](#velocity){reference-type="ref"
reference="velocity"}]{style="color: blue"}, we can observe that with
$\alpha=0.2$ the iteration converges fastest. In addition, if $\alpha$
is close to the boundary of convergence range, e.g. $\alpha_5=0.24$,
which is close to $0.245030$, the iteration will converge at a slow
speed.

::: center
2.3. *The case of $\alpha_{opt}$*
:::

The convergence curve with the optimal relaxation parameter
$\alpha_{opt}$ has been shown in Fig.
[[2](#optimal){reference-type="ref"
reference="optimal"}]{style="color: blue"}, which is close to a straight
line. In addition, we are concerned with the rate of convergence $\mu$.
In this work, the rate of convergence is calculated with errors. We
assume that the error $e$ finally converges to $e^*=0$, then the rate of
convergence is

![The convergence curve with the optimal relaxation parameter
$\alpha_{opt}$](../images/2_convergence_curve_for_optimal_alpha.png){#optimal
width="80%"}

$$\mu = \lim_{n\rightarrow\infty} \frac{\left|e_{n+1}-e^*\right|}{\left|e_n-e^*\right|}$$

After calculating, the rate of convergence is $0.632456$, which is
consistent with $\rho_{opt}$.

::: center
2.4. *Preconditioning*
:::

In this section, the diagonal elements of $A$ are selected to construct
the left preprocessor to acceleration the iteration, that is $$M=
    \left[\begin{array}{ll}
        6 & 0 \\
        0 & 4
    \end{array}\right]$$

Then the linear system is converted to $A'u=b'$ and the iteration
equation is converted into the following form:
$${u}_{k+1}= {u}_k+\alpha  {b'}-\alpha A'  {u}_k$$ where $A'=M^{-1}A$,
$b'=M^{-1}b$.

The numerical results of different cases have also been calculated with
the aid of MATLAB. The errors and residuals of the first 10 steps in the
form of 2-norm are listed in Table
[[3](#errors_preconditioning){reference-type="ref"
reference="errors_preconditioning"}]{style="color: blue"} and Table
[[4](#residuals_preconditioning){reference-type="ref"
reference="residuals_preconditioning"}]{style="color: blue"}
respectively. At the same time, the convergence curve of the first 50
steps are shown in Fig. [[3](#preconditioning){reference-type="ref"
reference="preconditioning"}]{style="color: blue"}.

::: {#errors_preconditioning}
  steps   $\alpha=0.06$   $\alpha=0.1$   $\alpha=0.2$   $\alpha=0.22$   $\alpha=0.24$   $\alpha=0.4$
  ------- --------------- -------------- -------------- --------------- --------------- --------------
  1       3.044524        2.967006       2.777139       2.739913        2.702961        2.418677
  2       2.935271        2.794920       2.479970       2.422739        2.367299        1.980909
  3       2.833570        2.641738       2.239383       2.169841        2.103336        1.656098
  4       2.738586        2.504058       2.036831       1.958631        1.884345        1.393763
  5       2.649587        2.379185       1.861292       1.776473        1.696255        1.175795
  6       2.565935        2.264987       1.706061       1.616037        1.531250        0.992841
  7       2.487068        2.159777       1.566904       1.472835        1.384651        0.838673
  8       2.412499        2.062216       1.441018       1.343932        1.253408        0.708556
  9       2.341799        1.971239       1.326444       1.227273        1.135360        0.598665
  10      2.274596        1.885992       1.221740       1.121327        1.028872        0.505832

  : $||e||_2$ of the first 10 steps after preconditioning
:::

[]{#errors_preconditioning label="errors_preconditioning"}

::: {#residuals_preconditioning}
  steps   $\alpha=0.06$   $\alpha=0.1$   $\alpha=0.2$   $\alpha=0.22$   $\alpha=0.24$   $\alpha=0.4$
  ------- --------------- -------------- -------------- --------------- --------------- --------------
  1       2.130860        2.015952       1.733854       1.678560        1.623730        1.209339
  2       1.974223        1.775241       1.354140       1.283106        1.216457        0.833142
  3       1.833244        1.574680       1.099779       1.030373        0.968219        0.659215
  4       1.706351        1.407427       0.925598       0.863209        0.808942        0.545436
  5       1.592114        1.267689       0.801686       0.746296        0.698547        0.457462
  6       1.489236        1.150576       0.709088       0.658984        0.615548        0.385431
  7       1.396536        1.051985       0.636260       0.589688        0.548845        0.325295
  8       1.312947        0.968492       0.576340       0.532008        0.492701        0.274727
  9       1.237498        0.897264       0.525275       0.482374        0.444049        0.232084
  10      1.169312        0.835981       0.480637       0.438722        0.401140        0.196083

  : $||r||_2$ of the first 10 steps after preconditioning
:::

[]{#residuals_preconditioning label="residuals_preconditioning"}

![The convergence curves after
preconditioning](../images/3_convergence_curve_for_every_alpha_after_preconditioning.png){#preconditioning
width="80%"}

::: center
2.5. *Discussion*
:::

![The comparision of preconditioning and
$\alpha_{opt}$](../images/4_convergence_curve_for_every_alpha_after_preconditioning.png){#compare
width="80%"}

Figure [[4](#compare){reference-type="ref"
reference="compare"}]{style="color: blue"} shows that the acceleration
effect of preconditioning is not obvious compared to the case of
$\alpha_{opt}$ discussed in "2.3. *The situation of $\alpha_{opt}$*". To
dig in the reason for this result, the preconditioned linear system
$A'u=b'$ is rewritten as $${u}_{k+1}=(I-\alpha A')  {u}_k+\alpha  {b'}$$

If ${b'}= {0}$, we get $${u}_k=(I-\alpha A')^k  {u}_0$$ And the
convergence of the Richardson iteration is determined by the property of
$(I-\alpha A')$. The spectral radius $\rho'=I-\alpha A'$ is $0.612372$,
which is close to $\rho_{opt}=0.632456$. This result indicates that the
rates of convergence of two cases are close, which explains the
acceleration effect of preconditioning is not obvious.

Besides, if we select $M^{-1}=A^{-1}$ as the left preconditioner, then
the exact solution can be calculated with only one step, that is
$${u}_{k+1} \equiv A^{-1}b ,\quad k=0,1,2, \ldots$$

**3. Readme of codes**

The codes of the corresponding solutions are contained in `Homework2.m`
and `Richardson_iteration.m`, which are developed by MATLAB. And the Git
repository is used during the development of codes. By typing *git log*
in the terminal of the current folder, the development process will be
displayed.

In addition, the program written by C language is also developed to
realise Richardson iteration (see `Richardson_iteration_C_version.c`),
with the aim of reviewing C language and make better use of C+PETSc.

::: center
REFERENCES
:::

::: thebibliography
10 Michelle Schatzman (2002), Numerical analysis: a mathematical
introduction, Clarendon Press, Oxford.
:::
