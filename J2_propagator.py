import numpy as np
from constants_1U import MU, J2, v_w_ei_i, R_EARTH, B_COEFF


def f1(pos, vel):
    '''Differential Equation - Function:
          d(x1) = x2
          -----           ,    x1 = pos,   x2 = v
           dt     
        Parameters:
        -----------
        pos: array
                This initializes the position state variable.
        vel: array
                This initializes the velocity state variable.
        Output:
        -------
        Returns instantaneous derivative of position.
    '''

    pos = vel
    return pos


def f2(pos, vel):
    '''Differential Equation - Function:
          d(x2) =  -u . x1       P             P 
          -----    -------   +     oblate  +     drag      , x1 = r, x2 = v
           dt       (x1)^3   
        Parameters:
        -----------
        pos: array
                This initializes the position state variable.
        vel: array
                This initializes the velocity state variable.
        Output:
        -------
        Returns instantaneous derivative of velocity.
    '''

    if DRAG_BOOLEAN == True:
        v_pertOblate = oblate_pertubations(pos)
        v_pertDrag = drag(pos, vel)
        v_pertTotal = (((-MU / ((norm(pos)) ** 3))) * pos) + v_pertOblate + v_pertDrag
        return v_pertTotal
    else:
        v_pertOblate = oblate_pertubations(pos)
        v_pertTotal = (((-MU / ((norm(pos)) ** 3))) * pos) + v_pertOblate
        return v_pertTotal


def oblate_pertubations(pos):
    ''' Calculates the Pertubations from the Effects of 
        Oblateness of Earth through J2 Model.
        J2 Perturbations: Pg 664 HD Curtis
        Parameters:
        -----------
        pos: array
                This initializes the position state variable.
        Output:
        -------
        Returns the perturbating accleration due 
        to Oblateness of the Earth.      
    '''

    scalar_pos = norm(pos)
    CONST_J2 = (1.5 * J2 * MU * R_EARTH ** 2) / (scalar_pos ** 4)
    x, y, z = pos
    k = ((5 * ((z / scalar_pos) ** 2) - 1) / scalar_pos)

    v_pertOblate = np.array([x * k, y * k, z * ((5 * ((z / scalar_pos) ** 2) - 3) / scalar_pos)])
    v_pertOblate = v_pertOblate * CONST_J2
    return v_pertOblate


def drag(pos, vel):
    '''Calculates the Pertubations from the Effects of Atmospheric Drag,
       assuming velocity of atmosphere at a particular point,
       is appx. equal to cross product of angular velocity of Earth 
       and position vector of that point.
       ~~~~~Drag Equation~~~~~

        Parameters:
        -----------
        pos: array
                This initializes the position state variable.
        vel: array
                This initializes the velocity state variable.
        Output:
        -------
        Returns the perturbating accleration due 
        to drag from Earth's Atmosphere.
    '''

    vel = np.array(vel, dtype=np.float64) * 1000  # Converting to m/s
    pos = np.array(pos, dtype=np.float64) * 1000  # Converting to m

    v_velAtm = np.cross(v_w_ei_i, pos)
    v_velRel = vel - v_velAtm  # Relative velocity being calculated.

    denst = density(pos / 1000)  # Requires argument in Km
    v_pertDrag = (-0.5) * denst * norm(v_velRel) * B_COEFF * v_velRel  # Accleration in m/s
    return (v_pertDrag / 1000)  # Returning Accleration in km/s


def density(pos):
    '''Calulates the atomosphereic density from a model based on 
       exponential decay, where the input is the height
       of the satellite from surface of the Earth,
       and the output is density of the atmosphere in kg/m^-3
       NOTE: The model is made specifically for LEO, i.e, 
       h ranging from 200km to 1000km.
       ~~~~~~Density Equation~~~~~~~
       Parameters:
        -----------
        pos: array
                This initializes the position state variable.
        Output:
        -------
        Returns the density of the atmosphere at that altitude.
    '''
    height = norm(pos) - R_EARTH  # Appx height of satellite from surface of Earth.

    I, alpha1, alpha2, alpha3, beta = (
    -55.80854359317351, 17771.64895643925, -3718462.067004107, 291861748.7626916, 0.008582907557446885)

    rho = np.exp(I + beta * height + alpha1 / height + alpha2 / (height ** 2) + alpha3 / (height ** 3))
    return rho  # In kg/m^3


def propagate(pos, vel, time, h_step_size=1, drag=False):
    """Propogate the State Variables.
       The Propogator here uses the initial State Vectors 
       and RK - 4 to approximate the subsequent system 
       of State Variables.
       RK - 4: X(n+1) = X(n) + h      (a + 2*b + 2*c + d)
                              ---  *  
                               6

       Parameters:
       -----------
       r: list of length (3)
            This initializes the position of the satellite.

       vel: list of length (3)
            This initializes the velocity of the satellite.

       time: floating-point number
            This initializes the time interval over which the propagator
            runs, in seconds.

       h_step_size: floating-point number, optional
            This sets the value of the step size, for RK-4.
            Default = 1

       drag: boolean, optional
            Takes Drag into consideration as a perturbation as well.  
            Default = False      
        Output:
        -------
        Returns the perturbating accleration due 
        to Oblateness of the Earth.
    """
    v_x, v_y = np.array(pos, np.float64), np.array(vel, np.float64)  # x is position, y is velocity
    steps = int(time / h_step_size)
    h = h_step_size
    global DRAG_BOOLEAN
    DRAG_BOOLEAN = drag

    for i in range(steps):
        v_ax, v_ay = RK4_a(v_x, v_y, h)
        v_bx, v_by = RK4_b(v_x, v_y, v_ax, v_ay, h)
        v_cx, v_cy = RK4_c(v_x, v_y, v_bx, v_by, h)
        v_dx, v_dy = RK4_d(v_x, v_y, v_cx, v_cy, h)
        v_x = v_x + ((h) * (v_ax + 2 * (v_bx + v_cx) + v_dx)) / 6
        v_y = v_y + ((h) * (v_ay + 2 * (v_by + v_cy) + v_dy)) / 6

    return v_x, v_y


def RK4_a(pos, vel, h):
    """ To Calculate Value of (a) of RK-4:
        v_a = f(v_x).
    """
    v_ax, v_ay = np.zeros(3), np.zeros(3)
    v_ax = f1(pos, vel)
    v_ay = f2(pos, vel)
    return v_ax, v_ay


def RK4_b(pos, vel, v_ax, v_ay, h):
    """ To Calculate Value of (b) of RK-4:
        v_b = f(v_x + (h/2)v_a)
    """
    v_l, v_m = np.zeros(3), np.zeros(3)
    v_l = pos + ((h / 2) * v_ax)
    v_m = vel + ((h / 2) * v_ay)
    v_bx, v_by = f1(v_l, v_m), f2(v_l, v_m)
    return v_bx, v_by


def RK4_c(pos, vel, v_bx, v_by, h):
    """ To Calculate Value of (c) of RK-4:
        v_c = f(v_x + (h/2)*v_b)
    """
    v_l, v_m = np.zeros(3), np.zeros(3)
    v_l = pos + ((h / 2) * v_bx)
    v_m = vel + ((h / 2) * v_by)
    v_cx, v_cy = f1(v_l, v_m), f2(v_l, v_m)
    return v_cx, v_cy


def RK4_d(pos, vel, v_cx, v_cy, h):
    """ To Calculate Value of (d) of RK-4:
        v_d = f(v_x + (h)*v_c)
    """
    v_l, v_m = np.zeros(3), np.zeros(3)
    v_l = pos + ((h) * v_cx)
    v_m = vel + ((h) * v_cy)
    v_dx, v_dy = f1(v_l, v_m), f2(v_l, v_m)
    return v_dx, v_dy


def norm(v_arr):
    ''' Calculates 2-norm of a vector.
        Parameters:
        -----------
        v_arr: array
                This initializes the vector variable.
        Output:
        -------
        Returns the 2-norm of the vector.  
    '''
    return (v_arr[0] ** 2 + v_arr[1] ** 2 + v_arr[2] ** 2) ** 0.5