# CVT Dynamics Derivation

![Full FBD](figures/CVT%20System.svg)

# Assumptions
- Constant temperature (no thermal expansion, friction coefficient of belt is constant)
- Engine torque scales linearly with throttle
- Forces and moments related to the CVT ratio are assumed to be quasi-static
    - Justification: shifting occurs at a slower rate than the acceleration of the vehicle itself
    - Benefit: shift ratio is independent of time, only on position and velocity at the current instant
    - Drawback: response of shift ratio deviates from reality, which is not quasi-static
    - True reason: this system would be extremely difficult to model without this assumption
- Friction between certain CVT components are considered negligible
    - Rollers on ramps/helix
    - Friction of sheaves moving during shift
- Belt never slips on secondary sheave

# Constants
Name|Unit|Description
---|:---:|---:
**Global**
$r_{absmin,p}$ | m | Inner radius of primary, bottom of where sheaves touch
$r_{absmin,s}$ | m | Inner radius of secondary, bottom of where sheaves touch
$\phi$ | rad | Half of the angle between sheaves (also applies to belt V-shape)
$L$ | m | Center to center distance between primary and secondary
$L_{b0}$ | m | Unstretched belt length
$b$ | m | Average width of belt V-shaped section
$A_b$ | m^2 | Belt cross sectional area
$m_b$ | kg | Total mass of belt
$E_b$ | Pa | Young's modulus of belt
$G_b$ | Pa | Shear modulus of belt
$\mu_b$ | none | Coefficient of friction of belt with sheaves
$N_g$ | none | Fixed gear ratio between the secondary and the wheels
$N_{fly}$ | none | Number of flyweight linkages in primary
$I_e$ | kg-m^2 | Total moment of inertia of spinning engine components
$I_p$ | kg-m^2 | Total moment of inertia of primary components with constant inertia values
$I_s$ | kg-m^2 | Total moment of inertia of secondary components with constant inertia values
**Primary CVT Tune**
$ramp(x)$ | function(m) -> m | Function for ramp height vs. x, where x is 0 at the highest ramp point.<br/>This function must be continuous, once differentiable, and $ramp'() < 0$.<br/>Function must be defined over range $[0, d_{p,max}]$.
$k_p$ | N/m | Linear spring constant for primary spring
$m_{fly}$ | kg | Mass of flyweights, total
**Secondary CVT Tune**
$k_s$ | N/m | Linear spring constant for secondary spring
$\kappa_s$ | N-m/rad | Torsional spring constant for secondary spring
$\theta_{0s}$ | rad | Angular pretension of secondary spring
$\theta_{helix}$ | rad | Angle of the secondary helix ramp
**Primary Subsystem**
$d_{p,max}$ | m | Max linear gap between primary sheaves
$d_{0p}$ | m | Primary spring initial displacement
$r_{cage}$ | m | Radius from primary axis to outer edge of primary ramp
$r_{shldr}$ | m | Radius from primary axis to flyweight arm pivot
$L_{arm}$ | m | Length of flyweight arm
$r_{roller}$ | m | Radius of rollers in primary
$x_{ramp}$ | m | Offset from flyweight arm pivot to furthest edge of ramp
$m_{flyarm}$ | kg | Mass of arm connected to flyweight
**Secondary Subsystem**
$d_{s,max}$ | m | Max linear gap between secondary sheaves
$d_{0s}$ | m | Secondary spring initial linear displacement
$r_{helix}$ | m | Radius of secondary helix ramp


# Derived Constants
Formula|Unit|Description
---|:---:|---:
$\rho_b = \frac{m_b}{A_b L_{b0}}$ | kg/m^3 | Density of the belt
$\theta_{s,max} = \frac{d_{s,max}}{r_{helix}\tan(\theta_{helix})}$ | rad | Max angular displacement of secondary sheave and torsional spring


# Time-Dependent Variables
Name|Unit|Description
---|:---:|---:
**Global**
$\tau_s$ | N-m | Torque load on the secondary
$\tau_p$ | N-m | Torque load on the primary
$\omega_p$ | rad/s | Angular velocity of primary
$\omega_s$ | rad/s | Angular velocity of secondary
$T_0$ | N | Slack-side belt tension (quasi-static)
$T_1$ | N | Taut-side belt tension (quasi-static)
$L_b$ | m | Current belt length (quasi-static)
**Primary Subsystem**
$d_p$ | m | Linear displacement of primary sheave during shift, range $[0,d_{p,max}]$, initially 0 (quasi-static)
$d_r$ | m | Linear displacement of roller from innermost edge of ramp (quasi-static)
$\theta_1$ | rad | Angle between flyweight arm and primary axis (quasi-static)
$\theta_2$ | rad | Angle between primary ramp surface normal at roller contact point and primary axis (quasi-static)
**Secondary Subsystem**
$d_s$ | m | Linear displacement of secondary sheave during shift, range $[0,d_{s,max}]$, initially 0 (quasi-static)


# Derived Variables
Formula|Unit|Description
---|:---:|---:
$r_p = d_p/\tan(\phi) + r_{absmin,p}$ | m | Pulley radius of primary
$r_s = (d_{s,max} - d_s)/\tan(\phi) + r_{absmin,s}$ | m | Pulley radius of secondary
$\alpha = 2\arccos(\frac{r_s-r_p}{L})$ | rad | Belt wrap angle around primary
$\beta = 2\arccos(\frac{r_p-r_s}{L})$  | rad | Belt wrap angle around secondary
$\theta_s = \frac{d_s}{r_{helix}\tan(\theta_{helix})}$ | rad | Angular twist of secondary spring during shift, range $[0,\theta_{s,max}]$
$r_{fly} = r_{shldr} + L_{arm} \sin(\theta_1)$ | m | Distance from flyweight to primary axis
$\theta_1 =$ solve numerically | rad | Angle of flyweight arm from primary axis
$\theta_2 =$ solve numerically | rad | Angle of surface normal of primary ramp at roller contact point
$d_r =$ solve numerically | m | Displacement of roller from rightmost edge of ramp
$F_f = \tau_p/r_p + \tau_s/r_s$ | N | Friction force at primary, assuming no slip
$N_p = \frac{\alpha(F_f)}{\sin(\phi)\ln(F_f + 1)} - \frac{\alpha}{\sin(\phi)}(T_0 - 1)$ | N | Net normal force applied to belt surface at primary


# Forces and Moments
Formula | Description
---|---:
$F_f \le N_p \mu_b$ | No-slip condition
$T_0 = E_b A_b * ((r_p\alpha + r_s\beta + 2\sqrt{L^2 - (r_p - r_s)^2})/L_{b0} - 1) - \rho_b A_b (r_s\omega_s)^2$ | Slack side tension
$F_f = T_1 - T_0 = \tau_p/r_p + \tau_s/r_s$ | Relation between slack and taut tension under no-slip condition
$F_f = T_1 - T_0 = \mu_b N_p + \tau_s/r_s$ | Relation between slack and taut tension under slipping condition
**Primary Subsystem**
$F_{sp} = k_p (d_{0p} + d_p)$ | Force from linear primary spring
$F_{bp} = \frac{\alpha(F_f)}{\tan(\phi)\ln(F_f + 1)} - \frac{\alpha}{\tan(\phi)}(T_0 - 1)$ | Force from belt
$F_{flyarm} = \frac{0.25 m_{fly}(r_{shldr} + L_{arm}\sin(\theta_1))\omega_p^2 L_{arm} \cos(\theta_1) \cos(\theta_2)}{L_{arm}\sin(\theta_1 + \theta_2) + r_{roller}\sin(2\theta_2)}$ | Force from flyweights and ramp
**Secondary Subsystem**
$F_{ss} = k_s (d_{0s} + d_s)$ | Force from linear secondary spring
$F_{bs} = \frac{\beta(F_f)}{\tan(\phi)\ln(F_f + 1)} - \frac{\beta}{\tan(\phi)}(T_0 - 1)$ | Force from belt
$T_{ss} = \kappa_s (\theta_{0s} + \theta_s)$ | Torque from torsional secondary spring
**Belt Drive**
$\tau_{ps} = \tau_p \frac{r_s}{r_p}$ | Torque applied to secondary from primary
$M_s = \tau_{ps} - \tau_s$ | Net moment applied to secondary


# Derivation of $r_p, r_s$

![Sheave Diagram](figures/belt%20on%20sheave.svg)


# Derivation of $\theta_s$

Tangent = opposite/adjacent

$\tan(\theta_{helix}) = \frac{d_s}{\theta_s r_{helix}}$

Solve for $\theta_s$

$\theta_s r_{helix} \tan(\theta_{helix}) = d_s$

$\theta_s = \frac{d_s}{r_{helix} \tan(\theta_{helix})}$


# Derivation of $F_b$ vs Tension

![Belt Cross-Section FBD](figures/belt%20cross%20section%20fbd.svg)

![Belt Section FBD](figures/belt%20tension%20fbd.svg)

$F_b$ depends on belt tension and sheave geometry:

$N = F_b/\cos(\phi) = R/\sin(\phi)$

$F_b = R/\tan(\phi)$

$R = F_b\tan(\phi)$ where R is a radial force on the sheave due to belt tension, and $\phi$ is half the V angle.

$R = N\sin(\phi)$

$F_b = N\cos(\phi)$

Something to note is that $R$ is influenced by centripetal force, $F_c$.

$F_c = mr\omega^2$

In the case of the belt, centripetal force along each element is based on this:

$dF_c = \rho A r d\theta\omega_p^2$

If the sheave is not slipping, then the static friction condition holds:

$F_f \le F_{f,max} = \mu_b N$, and $F_f$ is like a reaction force.
This reaction force is applied at every point along the belt contact surface. If the total moment from these reaction forces is balanced with the applied torque on the sheave and belt tension, then we can solve for F_f and verify that the static friction condition holds.

Taking sum of forces along the tangential and radial directions for a small slice of the belt contact area:

$\sum{F_{tangential}} = 0 = (-T) + (T + dT) - F_f$

$dT = F_f$

Applying small angle approximations:

$\sum{F_{radial}} = 0 = T d\theta + dT d\theta/2 - dR$

$dR = T d\theta + dT/2 d\theta$

Assuming that dT/2 is very small, approximately:

$dR = T d\theta$

$dN\sin(\phi) = T d\theta$

$\frac{dN}{d\theta}\sin(\phi) = T$

Assuming that $T(\theta)$ is in the form of an exponential, based on the derivations from WVU and Bestorq:

$T(\theta) = e^{C\theta} + Q$ where $T(0) = T_0$ and $T(\alpha) = T_1$

$Q = T_0 - 1$

$C = \frac{\ln(T_1-T_0+1)}{\alpha}$

$T(\theta) = (T_1 - T_0 + 1)^{\theta/\alpha} - 1 + T_0$

$\sin(\phi) \frac{dN}{d\theta} = (T_1 - T_0 + 1)^{\theta/\alpha} - 1 + T_0$

Integrating both sides wrt $d\theta$ with bounds $[0, \alpha]$ for $\theta$

$N \sin(\phi) = \int_0^\alpha{((T_1 - T_0 + 1)^{\theta/\alpha} - 1 + T_0) d\theta}$

$= \int_0^\alpha{((T_1 - T_0 + 1)^{\theta/\alpha}) d\theta} - \alpha(T_0 - 1)$

$= \int_0^\alpha{(e^{\ln(T_1 - T_0 + 1)\theta/\alpha}) d\theta} - \alpha(T_0 - 1)$

$= \frac{\alpha}{\ln(T_1 - T_0 + 1)} [e^{\ln(T_1 - T_0 + 1)\theta/\alpha}]_0^\alpha - \alpha(T_0 - 1)$

$= \frac{\alpha}{\ln(T_1 - T_0 + 1)} (e^{\ln(T_1 - T_0 + 1)} - 1) - \alpha(T_0 - 1)$

$= \frac{\alpha(T_1 - T_0)}{\ln(T_1 - T_0 + 1)} - \alpha(T_0 - 1)$

Using the fact that all torques must cancel:

$T_1 - T_0 = F_f$

$N \sin(\phi) = \frac{\alpha(F_f)}{\ln(F_f + 1)} - \alpha(T_0 - 1)$

$N = \frac{\alpha(F_f)}{\sin(\phi)\ln(F_f + 1)} - \frac{\alpha}{\sin(\phi)}(T_0 - 1)$

Substituting for $F_b$:

$F_b = N\cos(\phi)$

$F_b = \frac{\alpha(F_f)}{\tan(\phi)\ln(F_f + 1)} - \frac{\alpha}{\tan(\phi)}(T_0 - 1)$

For the secondary, replace $\alpha$ with $\beta$ and $\omega_p$ with $\omega_s$.

This gives us a relationship between slack side tension force and the normal force.

Once we can determine $T_0$, the rest of the system falls into place.

To determine whether the system will be slipping, we need to calculate $F_f$ assuming no slip:

$F_f = \tau_p/r_p + \tau_s/r_s$

If the friction "reaction" force $F_f > \mu_b N$, then the sheave is slipping and torque is limited by kinetic friction:

$F_f = \mu_b N$


# Derivation of $T_0$

We are basing the slack tension of the belt on the spring-like behavior of materials.
When the belt stretches, there is a static tension through the whole belt, which is T0.

$T_0 = E_b A_b * (L_b - L_{b0}) / L_{b0}$

$L_b = L_p + L_s + 2L_{mid}$

$L_b = r_p\alpha + r_s\beta + 2\sqrt{L^2 - (r_p - r_s)^2}$

Substituting $L_b$:

$T_0 = E_b A_b * ((r_p\alpha + r_s\beta + 2\sqrt{L^2 - (r_p - r_s)^2})/L_{b0} - 1)$

In addition, there is an effect due to centripetal force. The force effectively acts in opposition to slack side tension.

$T_c = \rho_b A_b (r_s\omega_s)^2$ Secondary-based linear belt speed is used, because the secondary is always assumed to have the no-slip condition.

(tec-science.com)

$T_0 = E_b A_b * ((r_p\alpha + r_s\beta + 2\sqrt{L^2 - (r_p - r_s)^2})/L_{b0} - 1) - \rho_b A_b (r_s\omega_s)^2$


# Derivation of $\theta_1, \theta_2, d_p$

The reaction force from the flyweight arm must be determined using kinematics. The system acts like a 2-bar linkage where the first link extends from the point "R" to the roller with length $L$ and angle $\theta_R$. The second link extends from the roller to the edge of the ramp with length $r_{roller}$ and angle $\theta_N$

![image](figures/two_bar_linkage.png)

An important note is that the ramp is defined by an arbitrary function, ramp(d), which returns a height from the base of the ramp. It is assumed that any result of this function will have the two-bar linkage in an "elbow-up" state, as shown in the figure. This is true if the slope of the ramp() function is always negative.

Solving this system is not trivial, since the function $ramp(d)$ is arbitrary. It means that the equations involving it cannot isolate its arguments, since the function is not known to be invertable.

Free variables (7): $x_1, y_1, x_2, y_2, \theta_1, \theta_2, d_r$

Knowns: $L_1, L_2, x_{ramp}, d_p, r_{cage}, r_{shldr}$ (Note: $L_1 = L_{arm}$ and $L_2 = r_{roller}$)

Arbitrary: $ramp(), ramp'()$

$x_1 = L_1\cos(\theta_1)$

$y_1 = L_1\sin(\theta_1)$

$x_1 = x_{ramp} + d_p - d_r$

$x_2 = x_1 + L_2\cos(\theta_2)$

$y_2 = r_{cage} - r_{shldr} - ramp(d_r - L_2\cos(\theta_2))$

$y_2 = y_1 + L_2\sin(\theta_2)$

$\theta_2 = \arctan(ramp'(d_r - L_2\cos(\theta_2))) + \pi/2$ where $ramp'(x) = \frac{d}{dx}ramp(x)$

With these 7 equations, we could in theory solve for every variable, if $ramp(x)$ wasn't arbitrary.
This nonlinear system must be solved numerically.
It can be simplified to this system of 3 equations where the variables are $d_r, \theta_1, \theta_2$:

eq1: $0 = L_1*\cos(\theta_1) - x_{ramp} + d_r$

eq2: $0 = L_1\sin(\theta_1) + L_2\sin(\theta_2) - r_{cage} + r_{shldr} + ramp(d_r - L_2\cos(\theta_2))$

eq3: $0 = \arctan(ramp'(d_r - L_2\cos(\theta_2))) + \pi/2 - \theta_2$

The penalty method can be used, where penalty (E) is minimized. Penalty is the sum of the right side of each equation squared.
The square is to make the value unsigned and to make larger deviations from zero have much more penalty.

$E = eq_1^2 + eq_2^2 + eq_3^2$

Bounds need to be placed on $\theta_1$ and $\theta_2$ such that they remain in the range $[-\pi,\pi]$, otherwise the numerical solver may not converge due to the symmetry of angles.
An alternative way to address the problem is to evaluate all $\theta$ as $\theta \mod 2\pi$.


# Derivation of $F_{flyarm}$

To find the reaction force $F_{flyarm}$, we need to use the centripetal force due to rotation of the primary.

$F_c = mr\omega^2$ where $m = m_{fly}/4, r = r_{fly}, \omega=\omega_p$

The two-bar linkage is usually applied in robotics, where the middle linkage can move. However, in this case the second link is a result of the roller on the ramp surface, so the entire ramp and roller piece can be treated as a single rigid body.

When solving for the forces and moments, forces in the x and y direction must be cancelled.
This is because the linkage arm itself does not move relative to its pivot point, it only rotates.
Because shifting occurs over a relatively long period of time, and the components of the flyweight arm assembly are so light compared to the rest of the assembly,
the quasi-static assumption can be used on this subsystem. The only information we would gain by not using it is the slight effect the momentum of the flyweight arm on shifting.
Most of the influence of momentum is captured by the flyweight's centripetal force, and not solving this system as quasi-static would be extremely difficult.

$\sum{F_x} = 0 = F_{flyarm} - N\cos(\theta_2)$

$\sum{F_y} = 0 = R_y - N\sin(\theta_2) + F_c$

$\sum{M_z} = 0 = F_c L_{arm} \cos(\theta_1) - N x_2 \sin(\theta_2) - N y_2 \cos(\theta_2)$

where $x_2 = L_{arm}\cos(\theta_1) + r_{roller}\cos(\theta_2)$

and $y_2 = L_{arm}\sin(\theta_1) + r_{roller}\sin(\theta_2)$

Solving for N in the moment equation:

$N = \frac{F_c L_{arm} \cos(\theta_1)}{x_2 \sin(\theta_2) + y_2 \cos(\theta_2)}$

The $F_y$ equation does not reveal any relevant information to the rest of the problem.

Plugging in $N$ into the $F_x$ equation:

$F_{flyarm} = \frac{F_c L_{arm} \cos(\theta_1) \cos(\theta_2)}{x_2 \sin(\theta_2) + y_2 \cos(\theta_2)}$

Expanding all intermediate values:

$F_{flyarm} = \frac{(0.25 m_{fly}r_{fly}\omega_p^2) L_{arm} \cos(\theta_1) \cos(\theta_2)}{\sin(\theta_2) (L_{arm}\cos(\theta_1) + r_{roller}\cos(\theta_2)) + \cos(\theta_2) (L_{arm}\sin(\theta_1) + r_{roller}\sin(\theta_2))}$

$F_{flyarm} = \frac{(0.25 m_{fly}r_{fly}\omega_p^2) L_{arm} \cos(\theta_1) \cos(\theta_2)}{(L_{arm}\cos(\theta_1)\sin(\theta_2) + r_{roller}\cos(\theta_2)\sin(\theta_2)) + (L_{arm}\sin(\theta_1)\cos(\theta_2) + r_{roller}\sin(\theta_2)\cos(\theta_2))}$

Reducing number of trig evaluations with identities:

$\sin(\theta_2)\cos(\theta_2) = \sin(2\theta_2)/2$

$\cos(\theta_1)\sin(\theta_2) + \cos(\theta_2)\sin(\theta_1) = \sin(\theta_1 + \theta_2)$

$F_{flyarm} = \frac{(0.25 m_{fly}r_{fly}\omega_p^2) L_{arm} \cos(\theta_1) \cos(\theta_2)}{L_{arm}\cos(\theta_1)\sin(\theta_2) + r_{roller}\sin(2\theta_2)/2 + L_{arm}\sin(\theta_1)\cos(\theta_2) + r_{roller}\sin(2\theta_2)/2}$

$F_{flyarm} = \frac{(0.25 m_{fly}r_{fly}\omega_p^2) L_{arm} \cos(\theta_1) \cos(\theta_2)}{L_{arm}\cos(\theta_1)\sin(\theta_2) + L_{arm}\sin(\theta_1)\cos(\theta_2) + r_{roller}\sin(2\theta_2)}$

$F_{flyarm} = \frac{0.25 m_{fly}r_{fly}\omega_p^2 L_{arm} \cos(\theta_1) \cos(\theta_2)}{L_{arm}\sin(\theta_1 + \theta_2) + r_{roller}\sin(2\theta_2)}$

$F_{flyarm} = \frac{0.25 m_{fly}(r_{shldr} + L_{arm}\sin(\theta_1))\omega_p^2 L_{arm} \cos(\theta_1) \cos(\theta_2)}{L_{arm}\sin(\theta_1 + \theta_2) + r_{roller}\sin(2\theta_2)}$


# Derivation of $d_p$

$\sum{F_x} = 0 = F_{sp} + F_{bp} - N_fly F_{flyarm}$

Expand values which depend on $d_p$:

$0 = k_p (d_{0p} + d_p) + F_{bp} - N_fly F_{flyarm}$

$k_p (d_{0p} + d_p) = 4F_{flyarm} - F_{bp}$

$d_p = (N_fly F_{flyarm} - F_{bp})/k_p - d_{0p}$


# Derivation of $d_s$

$\sum{F_x} = 0 = F_{ss} - F_{bs} + N_{helix} \cos(\theta_{helix})$

$\sum{M_z} = 0 = T_{ss} - r_{helix} N_{helix} \sin(\theta_{helix})$

Expand values:

$0 = k_s (d_{0s} + d_s) - F_{bs} + N_{helix} \cos(\theta_{helix})$

$0 = \kappa_s (\theta_{0s} + \frac{d_s}{r_{helix}\tan(\theta_{helix})}) - r_{helix} N_{helix} \sin(\theta_{helix})$

Solve moment equation for $N_{helix}$:

$r_{helix} N_{helix} \sin(\theta_{helix}) = \kappa_s (\theta_{0s} + \frac{d_s}{r_{helix}\tan(\theta_{helix})})$

$N_{helix}  = \kappa_s (\theta_{0s} + \frac{d_s}{r_{helix}\tan(\theta_{helix})})/(r_{helix} \sin(\theta_{helix}))$

Solve force equation for $d_s$:

$0 = k_s (d_{0s} + d_s) - F_{bs} + \cos(\theta_{helix}) \kappa_s (\theta_{0s} + \frac{d_s}{r_{helix}\tan(\theta_{helix})})/(r_{helix} \sin(\theta_{helix}))$

$0 = k_s (d_{0s} + d_s) - F_{bs} + \kappa_s (\theta_{0s} + \frac{d_s}{r_{helix}\tan(\theta_{helix})})/(r_{helix} \tan(\theta_{helix}))$

$0 = k_s d_{0s} + k_s d_s - F_{bs} + \frac{\theta_{0s} \kappa_s}{r_{helix}\tan(\theta_{helix})} + \frac{d_s \kappa_s}{(r_{helix}\tan(\theta_{helix}))^2}$

$k_s d_s + \frac{d_s \kappa_s}{(r_{helix}\tan(\theta_{helix}))^2} = -k_s d_{0s} + F_{bs} - \frac{\theta_{0s} \kappa_s}{r_{helix}\tan(\theta_{helix})}$

$d_s (k_s + \frac{\kappa_s}{(r_{helix}\tan(\theta_{helix}))^2}) = -k_s d_{0s} + F_{bs} - \frac{\theta_{0s} \kappa_s}{r_{helix}\tan(\theta_{helix})}$

$d_s = (-k_s d_{0s} + F_{bs} - \frac{\theta_{0s} \kappa_s}{r_{helix}\tan(\theta_{helix})})/(k_s + \frac{\kappa_s}{(r_{helix}\tan(\theta_{helix}))^2})$ Note how all the units work out to (N)/(N/m)


# Misc Notes

- Evaluate fitness of shift curve by optimizing SSE
- Tune accuracy of system by optimizing the sum of SSE optimizations for all test cases


# References

Skinner, Sean Sebastian. "Modeling and Tuning of CVT Systems for SAE® Baja Vehicles". West Virginia University, 2020.

BESTORQ. "Belt Theory". https://www.bestorq.com/techinfo.asp, Accessed Oct. 2024.

tec-science.com. "Centrifugal forces in the belt of a belt drive". https://www.tec-science.com/mechanical-power-transmission/belt-drive/centrifugal-forces/, Accessed Oct. 2024.