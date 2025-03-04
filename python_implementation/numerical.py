import numpy as np
import scipy.optimize
from typing import Callable, Union, Optional

class NumericalError(Exception):
    """Custom exception for numerical method errors."""
    pass

def bisection(f: Callable[[float], float], 
              x0: float, 
              x1: float, 
              max_iter: int = 20, 
              tol: float = 1e-6) -> float:
    """
    Bisection method for root finding.
    
    Args:
        f (callable): Function to find root of
        x0 (float): Lower bound of interval
        x1 (float): Upper bound of interval
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        float: Approximate root of the function
    
    Raises:
        NumericalError: If method fails to converge or input is invalid
    """
    def sign(x):
        return np.sign(x)
    
    # Input validation
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    if max_iter <= 0:
        max_iter = 1
    if x0 == x1:
        raise NumericalError("ERROR: Two guesses needed.")
    
    # Check that function changes sign in the interval
    if sign(f(x0)) * sign(f(x1)) >= 0:
        raise NumericalError("ERROR: Two guesses must give opposite signs.")
    
    # Actual bisection method
    for _ in range(max_iter):
        mid = (x0 + x1) / 2
        
        # Check if mid is close enough to a root
        if abs(f(mid)) <= tol:
            return mid
        
        # Update interval
        if sign(f(mid)) == sign(f(x0)):
            x0 = mid
        else:
            x1 = mid
    
    # Return midpoint if max iterations reached
    return (x0 + x1) / 2

def fixed_point(g: Callable[[float], float], 
                x0: float, 
                max_iter: int = 20, 
                tol: float = 1e-6) -> float:
    """
    Fixed point iteration method.
    
    Args:
        g (callable): Fixed point iteration function
        x0 (float): Initial guess
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        float: Approximate fixed point
    
    Raises:
        NumericalError: If method fails to converge
    """
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    if max_iter <= 0:
        max_iter = 1
    
    for _ in range(max_iter):
        x_new = g(x0)
        
        if not np.isfinite(x_new):
            raise NumericalError("ERROR: Unable to evaluate function.")
        
        if abs(x_new - x0) <= tol:
            return x_new
        
        x0 = x_new
    
    raise NumericalError(f"ERROR: Unable to converge within {max_iter} iteration(s).")

def newton_root(f: Callable[[float], float], 
                df: Optional[Callable[[float], float]] = None, 
                x0: float = 1.0, 
                max_iter: int = 20, 
                tol: float = 1e-6) -> float:
    """
    Newton's method for root finding.
    
    Args:
        f (callable): Function to find root of
        df (callable, optional): Derivative of function. If None, numerical derivative used.
        x0 (float): Initial guess
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        float: Approximate root
    
    Raises:
        NumericalError: If method fails to converge
    """
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    if max_iter <= 0:
        max_iter = 1
    
    # If no derivative provided, use secant method
    if df is None:
        return scipy.optimize.root_scalar(
            f, 
            method='secant', 
            x0=x0, 
            x1=x0*1.0001, 
            options={'maxiter': max_iter, 'xtol': tol}
        ).root
    
    # Newton's method with analytical derivative
    for _ in range(max_iter):
        fx = f(x0)
        dfx = df(x0)
        
        if dfx == 0:
            raise NumericalError("ERROR: Zero derivative.")
        
        dx = -fx / dfx
        
        if abs(dx) <= tol:
            return x0
        
        x0 += dx
        
        if not np.isfinite(x0):
            raise NumericalError("ERROR: Unable to evaluate function.")
    
    raise NumericalError(f"ERROR: Unable to converge within {max_iter} iteration(s).")

def secant_root(f: Callable[[float], float], 
                x0: float, 
                x1: float, 
                max_iter: int = 20, 
                tol: float = 1e-6) -> float:
    """
    Secant method for root finding.
    
    Args:
        f (callable): Function to find root of
        x0 (float): First initial guess
        x1 (float): Second initial guess
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        float: Approximate root
    
    Raises:
        NumericalError: If method fails to converge
    """
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    if max_iter <= 0:
        max_iter = 1
    if x0 == x1:
        raise NumericalError("ERROR: Two guesses needed.")
    
    for _ in range(max_iter):
        fx0, fx1 = f(x0), f(x1)
        dfx = fx1 - fx0
        
        if dfx == 0:
            raise NumericalError("ERROR: Zero denominator.")
        
        dx = -fx1 * (x1 - x0) / dfx
        
        if abs(dx) <= tol:
            return x0
        
        x0 = x1
        x1 += dx
        
        if not np.isfinite(x1):
            raise NumericalError("ERROR: Unable to evaluate function.")
    
    raise NumericalError(f"ERROR: Unable to converge within {max_iter} iteration(s).")

def vector_fixed_point(g: Callable[[np.ndarray], np.ndarray], 
                       x: np.ndarray, 
                       max_iter: int = 20, 
                       tol: float = 1e-6) -> np.ndarray:
    """
    Fixed point iteration for vector-valued functions.
    
    Args:
        g (callable): Vector-valued fixed point iteration function
        x (ndarray): Initial guess vector
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        ndarray: Approximate fixed point
    
    Raises:
        NumericalError: If method fails to converge
    """
    if max_iter < 1:
        max_iter = 1
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    
    for _ in range(max_iter):
        y = g(x)
        
        if not np.all(np.isfinite(y)):
            raise NumericalError("ERROR: Unable to evaluate function.")
        
        # Check relative error
        if np.all(np.abs((y - x) / np.clip(np.abs(x), 1e-100, None)) <= tol):
            return y
        
        x = y
    
    raise NumericalError(f"ERROR: Unable to converge within {max_iter} iteration(s).")

def vector_newton(f: Callable[[np.ndarray], np.ndarray], 
                  df: Optional[Callable[[np.ndarray], np.ndarray]] = None,
                  x: Optional[np.ndarray] = None, 
                  max_iter: int = 20, 
                  tol: float = 1e-6) -> np.ndarray:
    """
    Newton's method for vector-valued nonlinear systems.
    
    Args:
        f (callable): Vector-valued function to find root of
        df (callable, optional): Jacobian matrix function. If None, uses numerical approximation.
        x (ndarray, optional): Initial guess vector
        max_iter (int): Maximum number of iterations
        tol (float): Convergence tolerance
    
    Returns:
        ndarray: Approximate root vector
    
    Raises:
        NumericalError: If method fails to converge
    """
    if max_iter < 1:
        max_iter = 1
    if tol <= 0:
        raise NumericalError("ERROR: Invalid tolerance.")
    
    if x is None:
        # Create initial guess of zeros if not provided
        x = np.zeros(f(np.zeros(1)).size)
    
    for _ in range(max_iter):
        try:
            # Compute Jacobian (numerically if not provided)
            if df is None:
                jac = scipy.optimize.approx_fprime(x, f)
            else:
                jac = df(x)
            
            # Compute function value
            fx = f(x)
            
            if not np.all(np.isfinite(fx)):
                raise NumericalError("ERROR: Unable to evaluate function.")
            
            # Solve for correction step
            dx = np.linalg.solve(jac, -fx)
            
            # Update solution
            x += dx
            
            # Check convergence
            if np.linalg.norm(dx) <= tol * (1 + np.linalg.norm(x)):
                return x
        
        except (np.linalg.LinAlgError, ValueError) as e:
            raise NumericalError(f"ERROR in numerical solution: {str(e)}")
    
    raise NumericalError(f"ERROR: Unable to converge within {max_iter} iteration(s).")

# Optional: Numerical derivative if needed
def numerical_derivative(f: Callable[[float], float], 
                         x: float, 
                         h: float = 1e-8) -> float:
    """
    Compute numerical derivative using central difference method.
    
    Args:
        f (callable): Function to differentiate
        x (float): Point at which to compute derivative
        h (float): Step size for difference approximation
    
    Returns:
        float: Approximate derivative
    """
    return (f(x + h) - f(x - h)) / (2 * h)

# Export selected functions to match C++ interface
def deriv(f: Callable[[float], float], x: float) -> float:
    """
    Wrapper for numerical derivative to match C++ function signature.
    
    Args:
        f (callable): Function to differentiate
        x (float): Point at which to compute derivative
    
    Returns:
        float: Approximate derivative
    """
    return numerical_derivative(f, x)
