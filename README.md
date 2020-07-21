# Dispersion_MTSI

In this project we will try to write some Python code to solve the dispersion relation and study the MTSI mode in different situations.


## The MTSI Dispersion relation

From S. Janhunen et. al, Physics of Plasmas 25, 082308 (2018);  doi: 10.1063/1.5033896

We have to solve :
<!-- $$1 -  \frac{\omega_{pi}^2}{\omega^2} - \frac{\omega_{pe}^2 k_z^2}{ (\omega - k_y v_0)^2 k^2} - \frac{\omega_{pe}^2 k_y^2}{ ((\omega - k_y v_0)^2 - \Omega_{ce}^2 ) k^2} = 0$$ -->
![equation](https://latex.codecogs.com/gif.latex?1%20-%20%5Cfrac%7B%5Comega_%7Bpi%7D%5E2%7D%7B%5Comega%5E2%7D%20-%20%5Cfrac%7B%5Comega_%7Bpe%7D%5E2%20k_z%5E2%7D%7B%20%28%5Comega%20-%20k_y%20v_0%29%5E2%20k%5E2%7D%20-%20%5Cfrac%7B%5Comega_%7Bpe%7D%5E2%20k_y%5E2%7D%7B%20%28%28%5Comega%20-%20k_y%20v_0%29%5E2%20-%20%5COmega_%7Bce%7D%5E2%20%29%20k%5E2%7D%20%3D%200)

with $\omega_{pi}$ and $\omega_{pe}$ the ion and electron plasma frequencies, respectively, $\omega = \omega_r + i \gamma$ the complex pulsation, $k$ the norm of the wave vector $\vec{k} = k_x \vec{e_x} + k_y \vec{e_y} + k_z \vec{e_z}$.
The direction $z$ is parallel to the magnetic field $B$, $x$ is in the direction of the electric field $E$, and $y$ is in the direction of the $E x B$ drift.
The drift velocity $\vec{v_0} = \frac{E x B}{B^2} = \frac{E}{B} \vec{e_y}$.
Lastly, $\Omega_{ce} = \frac{q B}{m_e}$ is the electron cyclotron frequency.

## Approche

The usual approache is to

   1. Normalise the relation, using the usual $\lambda_{De}$ Debye lengh and the Ions plasma pulsation
   2. Fixe the wave vector $\vec{k}$
   3. Solve for the complex pulsation  $\omega = \omega_r + i \gamma$
   4. I iterate step 2. and 3. over the parameter space needed.

## Normalised Relation

The wave vectors are normalized by the debye length $$\lambda_{De} = \sqrt{ \frac{\epsilon_0 k_B T_e}{n_e q_e^2}}$$
with $\epsilon_0$ the vacuum permitivity, $k_B$ the Boltzman constant, $T_e$ the electron temperature (in Kelvin), $n_e$ the electron density, and $q_e$ the elementary charge.

The pulsation are normalized by the ion plasma pulsation
$$\omega_{pi} = \sqrt{  \frac{n_e q_e^2}{m_i \epsilon_0}}$$
with $m_i$ the ion mass.

The velocities are normalized by
$$\lambda_{De} \omega_{pi} = \sqrt{\frac{k_B T_e}{m_i}} = u_B$$ the Bohm speed.

The normalized equation reads
<!-- $$ 1 - \frac{1}{\tilde{\omega}^2} - \frac{m_i}{m_e} \frac{ \tilde{k_z}}{ (\tilde{\omega} - \tilde{k} \tilde{v_0})^2 \tilde{k}^2} - \frac{m_i}{m_e} \frac{ \tilde{k_y}}{ ((\tilde{\omega} - \tilde{k} \tilde{v_0})^2 - \frac{\Omega_{ce}^2}{\omega_{pi}^2}) \tilde{k}^2} =0$$ -->
![equation](https://latex.codecogs.com/gif.latex?1%20-%20%5Cfrac%7B1%7D%7B%5Ctilde%7B%5Comega%7D%7D%20-%20%5Cfrac%7Bm_i%7D%7Bm_e%7D%20%5Cfrac%7B%20%5Ctilde%7Bk_z%7D%7D%7B%20%28%5Ctilde%7B%5Comega%7D%20-%20%5Ctilde%7Bk%7D%20%5Ctilde%7Bv_0%7D%29%5E2%20%5Ctilde%7Bk%7D%5E2%7D%20-%20%5Cfrac%7Bm_i%7D%7Bm_e%7D%20%5Cfrac%7B%20%5Ctilde%7Bk_y%7D%7D%7B%20%28%28%5Ctilde%7B%5Comega%7D%20-%20%5Ctilde%7Bk%7D%20%5Ctilde%7Bv_0%7D%29%5E2%20-%20%5Cfrac%7B%5COmega_%7Bce%7D%5E2%7D%7B%5Comega_%7Bpi%7D%5E2%7D%29%20%5Ctilde%7Bk%7D%5E2%7D%20%3D0)

with $\tilde{\omega} = \frac{\omega}{\omega_{pi}}$, $\tilde{k} = k \lambda_{De}$, and $\tilde{v_0} = \frac{v_0}{u_B}$.
The normalized equation is not much simpler than the previous version, but now the variables are of the order of unity.
