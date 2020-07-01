# dispersion_MTSI

In this project we will try to write some Python code to solve the dispersion relation and study the MTSI mode in different situations.


## The MTSI Dispersion relation

From S. Janhunen et. al, Physics of Plasmas 25, 082308 (2018);  doi: 10.1063/1.5033896

We have to solve :
$$ 1 -  \frac{\omega_{pi}^2}{\omega^2} - \frac{\omega_{pe}^2 k_z^2}{ (\omega - k_y v_0)^2 k^2} - \frac{\omega_{pe}^2 k_y^2}{ ((\omega - k_y v_0)^2 - \Omega_{ce}^2 ) k^2} = 0 $$ 

with $\omega_{pi}$ and $\omega_{pi}$ the ion and electron plasma frequencies, respectively, $\omega = \omega_r + i \gamma$ the complex pulsation, $k$ the norm of the wave vector $\vec{k} = k_x \vec{e_x} + k_y \vec{e_y} + k_z \vec{e_z}$.
The direction $z$ is parallel to the magnetic field $B$, $x$ is in the direction of the electric field $E$, and $y$ is in the direction of the $E x B$ drift.
The drift velocity $\vec{v_0} = \frac{E x B}{B^2} = \frac{E}{B} \vec{e_y}$.
Lastly, $\Omega_{ce} = \frac{q B}{m_e}$ is the electron cyclotron frequency. 

## Approche

The usual approache is to 

   1. Normalise the relation, using the usual $\lambda_{De}$ Debye lengh and the Ions plasma pulsation
   2. Fixe the wave vector $\vec{k}$
   3. Solve for the complex pulsation  $\omega = \omega_r + i \gamma$
   4. Itterate step 2. and 3. over the parameter space needed.

## Normalised Relation


   