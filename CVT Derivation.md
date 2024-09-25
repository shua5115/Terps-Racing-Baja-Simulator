# CVT Dynamics Derivation

Main reference: Skinner, Sean Sebastian. Modeling and Tuning of CVT Systems for SAE® Baja Vehicles. West Virginia University, 2020.

## Variables

### Constants
Name|Unit|Description
---|:---:|---:
$r_p$     | m    | inner radius of primary for belt
$R_p$     | m    | outer radius of primary for belt
$r_s$     | m    | inner radius of secondary for belt
$R_s$     | m    | outer radius of secondary for belt
$\phi$    | rad  | angle between sheaves
$g_p$     | m    | max gap between primary sheaves
$g_s$     | m    | max gap between secondary sheaves
$L$       | m    | center to center distance between primary and secondary
$l_b$     | m    | belt length
$A_b$     | m<sup>2</sup> | belt cross sectional area
$m_b$     | kg   | total mass of belt
$E_b$     | Pa   | Young's modulus of belt
$G_b$     | Pa   | shear modulus of belt
$\mu_b$   | none | Coefficient of friction of belt with sheaves
$x_{0p}$ | m    | Primary spring initial displacement
$x_{0s}$ | m    | Secondary spring initial displacement
$r_{helix}$ | m  | Radius of secondary helix ramp
$m_{fly,arm}$ | kg | Mass of arm connected to flyweight
$R_g$     | none | Fixed gear ratio between the secondary and the wheels

### Derived Constants
Formula|Unit|Description
---|:---:|---:
$\mu_e=\mu_b/\sin(\frac{\phi}{2})$ | none | Effective friction coefficient between belt and sheaves


### CVT Tune Variables
Name|Unit|Description
---|:---:|---:
$ramp(x)$ | function | Function for ramp height vs. x, where x and height are 0 at the beginning of shift, all units in meters.
$k_p$     | N/m  | Linear spring constant for primary spring
$m_{fly}$ | kg   | Mass of flyweights, total
$k_s$     | N/m  | Linear spring constant for secondary spring
$\kappa_s$  | N-m/rad | Torsional spring constant for secondary spring
$\gamma_s$ | rad  | Angular pretension of secondary spring
$\theta_{helix}$ | rad | Angle of the secondary helix ramp

### Variables
Name|Unit|Description
---|:---:|---:
$\tau_s$ | N-m | torque load on the secondary
$\tau_p$ | N-m | torque load on the primary
$\omega_p$ | rad/s | angular velocity of primary
$\omega_s$ | rad/s | angular velocity of secondary
$r_p$ | m | pulley radius of primary
$r_s$ | m | pulley radius of secondary
$d_p$ | m | compression of primary spring during shift, range $[0,g_p]$, initially 0
$d_s$ | m | compression of secondary spring during shift, range $[0,g_s]$, initially 0


### Derived Variables
Formula|Unit|Description
---|:---:|---:
$\alpha = 2\arccos(\frac{r_s-r_p}{L})$ | rad | belt wrap angle around primary
$\beta = 2\arccos(\frac{r_p-r_s}{L})$  | rad | belt wrap angle around secondary
$T_0 = \frac{2\sin(\frac{\beta}{2})*[2F_{clamp}\tan(\frac{\phi}{2})] + \frac{m_b*r^2*\omega^2}{12}}{\cos(\frac{\beta-\pi}{2})*(e^{\mu_e\beta} + 1)}$ | N | Slack side tension for either pulley. $F_{clamp}$, $r$, and $\omega$ must be replaced with values for the primary or secondary.
$T_1 = T_0 e^{\mu_e\beta}$ | N | Taut side tension force derived from slack side tension, applicable for primary and secondary
$F_{clamp,s} = \frac{1}{\beta}[k_s*(x_{0s} + x_s) + \frac{1}{2\tan(\theta_{helix})}(\kappa_s*(\gamma_s + \theta_s) + \tau/r_{helix})]$ | N | Clamping force at secondary

## FBDs
(optional) Overall (engine, gearcase input shaft)

### Primary sheaves (belt replaced with tension forces)

### Secondary sheaves (vice versa)

### Primary internals (shaft, rollers, linkages, flyweights, spring, sheaves)

![image](figures/Primary%20Internals%20FBD.png)
![image](figures/Primary%20Flyweight%20FBD.png)
The following FBD mixes up symbols due to being from a different source.
The formulas below use the correct symbols.
![image](figures/Simple%20Belt%20FBD.png)


$F_{spring} = k_p*(d - d_p)$ where d_p is the initial displacement of the primary spring

$F_{sheave}$ depends on belt tension and sheave geometry:

$R = N\sin(\phi/2)$ where R is a vertical force due to tension, and $\phi$ is the V angle. 

$N = F_{sheave}/\cos(\phi/2)$

$R = F_{sheave}\frac{\sin(\phi/2)}{\cos(\phi)} = F_{sheave}\tan(\phi/2)$

$F_{sheave} = R/\tan(\phi/2)$

R results from a distributed load of belt tension. From the Bestorq belt theory reference:

$R = T_0 \alpha$, where $T_0$ is the slack-side tension and $\alpha$ is the wrap angle around the sheave ($\alpha$ for primary, $\beta$ for secondary)

$F_{sheave} = T_0 \alpha/\tan(\phi/2)$

Finally, Rx, the reaction force from the roller linkage, must be determined using kinematics. The system acts like a 2-bar linkage where the first link extends from the point "R" to the roller with length $L$ and angle $\theta_R$. The second link extends from the roller to the edge of the ramp with length $r_{roller}$ and angle $\theta_N$

![image](figures/two_bar_linkage.png)

$D=\frac{x_2^2+y_2^2-L_1^2-L_2^2}{2 L_1 L_2}$

$\theta_2=\arctan(s_{elbow}\sqrt{1-D^{2}}, D)$, where $s_{elbow}$ is either -1 (elbow up)or 1 (elbow down)

$\theta_1=\arctan(y_2,x_2)-\arctan(L_2\sin(\theta_2),L_1+L_2\cos(\theta_2))$

An important note is that the ramp is defined by an arbitrary function, ramp(d), which returns a height from the base of the ramp. It is assumed that any result of this function will have the two-bar linkage in an "elbow-up" state.

Solving this system is not trivial, since $ramp(d)$ depends on the unknown $d$, which prevents the system from having any givens. We have constraints which must hold, though:

Free variables (7): $x_1, y_1, x_2, y_2, \theta_1, \theta_2, d$
Knowns: $L_1, L_2, ramp(d), ramp'(d), x_0, r_{cage}$

$x_1 = L_1\cos(\theta_1)$

$y_1 = L_1\sin(\theta_1)$

$x_2 = x_0 - d$

$x_2 = x_1 + L_2\cos(\theta_1+\theta_2)$

$y_2 = r_{cage} - ramp(d)$

$y_2 = y_1 + L_2\sin(\theta_1+\theta_2)$

$\theta_1 + \theta_2 = \arctan(ramp'(d)) + \frac{\pi}{2}$ where $ramp'(x) = \frac{d}{dx}ramp(x)$

With these 7 equations, we can in theory solve for every variable, but try as I might, I can't reach an analytic solution. This is mostly due to ramp(x) not being a known function. So, this nonlinear system must be solved numerically. It can be simplified to this system of 3 equations where the variables are $d, \theta_R, \theta_N$:

$\theta_R = \theta_1$

$\theta_N = \theta_1 + \theta_2$

eq1: $L_1*\cos(\theta_R) + L_2*\cos(\theta_N) - x_0 + d$
eq2: $L_1*\sin(\theta_R) + L_2*\sin(\theta_N) - r_{cage} + ramp(d)$
eq3: $\arctan(ramp'(d)) + \pi/2 - \theta_N$

Where each equation is set equal to 0.

The penalty method can be used, where penalty is minimized:

$E = eq_1^2 + eq_2^2 + eq_3^2$

Bounds need to be placed on $\theta_R$ and $\theta_N$ such that they remain in the range $[-\pi,\pi]$, otherwise the numerical solver may not converge due to the symmetry of angles.

The next step is to solve for the forces of the system:

To find the reaction force $R_x$, we need to use the centripetal force due to rotation of the primary.

$F_c = mr\omega^2$ where $m = m_{fly}/4, r = r_{fly}, \omega=\omega_p$

The two-bar linkage is usually applied in robotics, where the middle linkage can move. However, in this case the second link is a result of the roller on the ramp surface, so the entire ramp and roller piece can be treated as a single rigid body.

Solving forces and moments, with the quasi-static assumption:

$\sum{F_x} = R_x - N\cos(\theta_N)$

$\sum{F_y} = R_y - N\sin(\theta_N) + F_c$

$\sum{M_z} = F_c L \cos(\theta_R) - N x_2 \sin(\theta_N) - N y_2 \cos(\theta_N)$
where $x_2 = L\cos(\theta_1) + r_{roller}\cos(\theta_2)$
and $y_2 = L\sin(\theta_1) + r_{roller}\sin(\theta_2)$

With the constraints $\sum{F_y} = 0, \sum{M_z} = 0$, since the roller cannot move in the y direction and the roller's rotation can be treated as quasi-static.

Solving for N in the moment equation:

$N = \frac{F_c L \cos(\theta_R)}{x_2 \sin(\theta_N) + y_2 \cos(\theta_N)}$

With some careful observation, we actually do not care about the $F_y$ equation at all. It means nothing to us, and never will. We can just plug N directly into the $F_x$ equation to get 



### Secondary internals (shaft, rollers, helix ramp, spring, sheaves)

Assumptions
- Constant temperature (no effects of temp on friction or size of objects)
- Engine torque scales linearly with throttle

Time-Dependent Variables
- Primary w
- Secondary w
- Primary sheave displacement
- Secondary sheave displacement
- Belt length?
- Torque at primary
- Torque at secondary

## Equations

Taut side tension equivalent between primary and secondary:

$T_{1p} = T_{1s}$



