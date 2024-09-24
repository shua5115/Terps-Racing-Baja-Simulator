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
$\tau$ | N-m | torque load on the secondary
$\omega_p$ | rad/s | angular velocity of primary
$\omega_s$ | rad/s | angular velocity of secondary
$r_p$ | m | pulley radius of primary
$r_s$ | m | pulley radius of secondary
$x_p$ | m | compression of primary spring during shift, range $[0,g_p]$, initially 0
$x_s$ | m | compression of secondary spring during shift, range $[0,g_s]$, initially 0


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

Primary sheaves (belt replaced with tension forces)

Secondary sheaves (vice versa)

Primary internals (shaft, rollers, linkages, flyweights, spring, sheaves)

Secondary internals (shaft, rollers, helix ramp, spring, sheaves)


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



