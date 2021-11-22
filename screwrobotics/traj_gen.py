from typing import Callable, List, Tuple
import numpy as np
from . import representations as rep



def quintic_time_scaling(t: float, T_f: float = 1) -> float:
    """
    Quintic time scaling function.

    :param t: current time
    :param T_f: final time
    :return: time scaling factor
    """
    t = t/T_f
    return 10 * t**3 - 15 * t**4 + 6 * t**5

def screw_trajectory(start: np.array, end: np.array, N: int,
                     scaling: Callable[[float, float], float] = quintic_time_scaling, 
                     T_f: float = 1) -> List[Tuple[float, np.array]]:
    """
    Generates a trajectory from a start transformation to an end transformation using 
    the "screw path", made from interpolation of the transfroms on the lie algebra.

    :param start: start transformation
    :param end: end transformation
    :param N: number of discretisation points
    :param scaling: time scaling function. Default: quintic scaling
    :param T_f: final time
    :return: trajectory
    """
    ts = scaling(np.linspace(0, T_f, N))
    trajectory = []
    for t in ts:
        lie_diff = np.linalg.pinv(start) @ end
        transform = start @ (rep.Transformation.from_SE3(lie_diff) * t).SE3
        trajectory.append((t, transform))
    return trajectory


