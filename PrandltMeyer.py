import numpy as np

def PMF(M,Gamma):
    Gp=Gamma+1 # γ-1
    Gm=Gamma-1 # γ+1
    if Gamma>1:
        nu = np.sqrt(Gp/Gm)*np.arctan(np.sqrt(Gm*(np.power(M,2)-1)/Gp))-np.arctan(np.sqrt(np.power(M,2)-1)) # Pranldt-Meyer angle in rad
        nu = np.rad2deg(nu) # Pranldt-Meyer angle in degrees
    elif Gamma==1:
        nu = np.arctan(1/(np.sqrt(np.power(M,2)-1)))+np.deg2rad(np.sqrt(np.power(M,2)-1))-np.deg2rad(45) # Pranldt-Meyer angle in rad
        nu = np.rad2deg(nu) # Pranldt-Meyer angle in degrees
    else:
        nu = complex(0,-1)/np.sqrt(complex(Gm/Gp,0))*np.arctanh(-np.sqrt(complex(Gm/Gp,0))/complex(0,1)*np.sqrt(np.power(M,2)-1))-np.arctan(np.sqrt(np.power(M,2)-1)) # Pranldt-Meyer angle in rad
        assert(np.imag(nu.any())==0),"Prandlt-Mayer nu angle calculation error : The imaginary part is not zero" # Stop excecution if Prandlt-Meyer angle is not real number
        nu = np.rad2deg(np.real(nu)) # Pranldt-Meyer angle in degrees
    mu = np.arcsin(1/M) # Mach angle in rad
    mu = np.rad2deg(mu) # Mach angle in degrees
    return M,nu,mu

def goal_seek(target,_threshold,Gamma): # Gold seek algorithm for finding the proper Mach number for a know Prandlt-Meyer angle
    Mmax = 4; Mmin = 1
    numax = PMF(Mmax,Gamma)[1]; numin = PMF(Mmin,Gamma)[1]
    fmax = target - numax; fmin = target - numin
    Mprev = Mmin; fprev = fmin; 
    L_ = fmin; R_ = fmax; Range = abs(Mmax-Mmin)
    while Range >= _threshold:
        Mnext = (Mmin*R_-Mmax*L_)/(R_-L_)
        nunext = PMF(Mnext,Gamma)[1]; fnext = target - nunext
        if fnext == 0:
            break
        elif fmin*fnext<=0:
            Mmax = Mnext; fmax = fnext; R_ = fnext
            if fprev*fnext>0:
                L_ = L_/2
        else:
            Mmin = Mnext; fmin = fnext; L_ = fnext
            if fprev*fnext>0:
                R_ = R_/2
        Mprev = Mnext; fprev = fnext; Range = (abs(Mmax-Mmin))
#         print(f'Range is: {Mmin}  ----  {Mnext}   ----   {Mmax}, ----> {fnext} --------------------------------',end="\r", flush=True)
    return Mnext

def secant(x1,x2,f,tol, *args):
    err = 1
    stepNum = 0
    while err > tol:
        stepNum=stepNum+1
        A = np.array([[1,1],[f(x1), f(x2)]])
        B = np.array([[1], [0]])
        p = np.matmul(np.linalg.inv(A), B)
        x_n = float(sum([k1*k2 for k1,k2 in zip(p, [x1,x2])]))
        if 'show' in args:
            print(x_n, x1, x2)
        if x_n < 0 and 'only_possitive' in args:
            x_n = float(np.random.rand(1)*(x1 + x2)/2)
        err = abs(x2 - x_n)
        x1,x2 = x2, x_n
        # print(x1, x2)
    
    return x_n, f(x_n), stepNum


# Inputs
def MNM(M,nu,mu,Gamma): # Mach-Nu-Mu Function: One input only must be a number, the other values must be zero
    # ex: M = 2, mu = 0, nu = 0, Gamma = 1.4 or M = 0,mu = 30, nu = 0, Gamma 1.4 etc.
    M = np.array(M)
    mu = np.array(mu)
    nu = np.array(nu)

    Error = True
    if M.any()!=0 and nu.any()==0 and mu.any()==0:
        Error = False
        [M,nu,mu] = PMF(M,Gamma)
    elif mu.any()!=0 and M.any()==0 and nu.any()==0:
        Error = False
        M = 1/np.sin(np.deg2rad(mu))
        [M,nu,mu] = PMF(M,Gamma)
    elif nu.any()!=0 and M.any()==0 and mu.any()==0:
        Error = False
        nu = [nu] if np.size(nu) == 1 else nu
        M = np.array([secant(1, 1.5, lambda x: PMF(x,Gamma)[1] - y, 1e-5, 'only_possitive')[0] for y in nu])
        mu = np.rad2deg(np.arcsin(1/M))    
        if np.size(M) == 1:
            M, nu, mu = M[0], nu[0], mu[0]

        # if np.size(nu)==1:
        #     M = goal_seek(nu,1e-6,Gamma)
        #     mu = np.arcsin(1/M)
        #     mu = np.rad2deg(mu)
        # else:
        #     M = np.ones(len(nu))
        #     mu = np.ones(len(nu))
        #     for i in range(0,len(nu)):
        #         M[i] = goal_seek(nu[i],1e-6,Gamma)
        #         mu[i] = np.arcsin(1/M[i])
        #         mu[i] = np.rad2deg(mu[i])
    else:
        if Error == True:
            raise ValueError("Error! Only one input must be not zero - the other 2 must be equal with zero")
    return M,nu,mu

