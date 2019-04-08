import numpy as np
import scipy
import scipy.linalg
import scipy.optimize
import fb
import expm

def FBLLH(theta, gamma, O, A, B, n, method='hg', withvol=False):
    #p = (len(parameters) + 1) / 2
    #theta = parameters[:p-1] # length p-1
    #theta = np.concatenate([0], theta)
    #gamma = parameters[p-1:] # length p

    p = len(theta)

    # get the parameters sorted before we evaluate anything
    #sorting_indices = np.argsort(theta)
    #sorted_theta = theta[sorting_indices]
    #sorted_gamma = gamma[sorting_indices]

    alpha = np.concatenate((theta, gamma))

    if method == 'hg':
        c_res = fb.hgm_FB_2(alpha, withvol=withvol)
    elif method == 'hg_linear':
        c_res = fb.hgm_FB(alpha, withvol=withvol)
    elif method == 'MC':
        c_res = fb.G_FB(alpha, method='MC', withvol=withvol)
    elif method == 'SPA':
        #c_res = fb.saddleapprox_FB_revised(np.array(sorted(alpha[:p]))+1.0, M=alpha[p:p*2]/2.0, dub=1, order=3) * np.exp(1)
        c_res = fb.SPA(alpha, withvol=withvol)
    elif method == 'lpinv':
        raise ValueError('Not sure what lpinv is supposed to do.')
        c_res = lpinv(alpha)
    else:
        raise ValueError('Not a valid method.')

    c = c_res[0]

    #l = -n*np.log(c) - np.sum(A.dot(O.T).dot(np.diag(theta)).dot(O) + gamma.dot(B.T).dot(O))
    l = -n*np.log(c) - np.trace(O.dot(B).dot(gamma[:,None].T) + A.dot(O.T).dot(np.diag(theta)).dot(O))

    #grad = -n*c_res[1:] / c + np.concatenate((-np.diag(O.dot(A).dot(O.T)),O.dot(B)))
    theta_grad = np.diag(O.dot(A).dot(A.T))
    gamma_grad = np.squeeze(B.T.dot(O.T))
    param_grad = np.concatenate((theta_grad, gamma_grad))
    grad = -n*c_res[1:] / c + param_grad

    return l, grad

# I think this is useless in the context of this problem...
def orthogonal_matrix(x):
    a = x[:,None]
    q = np.identity(len(x)) - 2.0*a.dot(a.T)/float(a.T.dot(a))
    return q

def grad_mle_update(theta, gamma, O, A, B, n, method='hg', withvol=False, min_step=1e-3):
    p = len(theta)

    # get the initial llh and gradient
    l, grad = FBLLH(theta, gamma, O, A, B, n, method=method, withvol=withvol)

    # one of the thetas must be zero
    # assuming it's theta[0], let's fix that component of the gradient to zero
    grad[0] = 0

    # try stepping once
    delta = 1
    new_theta = theta + delta*grad[:p]
    new_theta[0] = 0
    new_gamma = gamma + delta*grad[p:]
    new_l, new_grad = FBLLH(new_theta, new_gamma, O, A, B, n, method=method, withvol=withvol)

    # try some smaller steps
    while new_l < l and delta > min_step:
        delta /= 2.0
        new_theta = theta + delta*grad[:p]
        new_theta[0] = 0
        new_gamma = gamma + delta*grad[p,:]
        new_l, new_grad = FBLLH(new_theta, new_gamma, O, A, B, n, method=method, withvol=withvol)

    return new_theta, new_gamma, new_l, new_grad, delta

def grad_mle_orth_update(theta, gamma, O, A, B, n, method='hg', withvol=False, tol=1e-3):
    m_theta = np.diag(theta)
    #AA = m_theta.dot(O).dot(A).dot(O.T).dot(m_theta) + gamma.dot(B.T).dot(O.T)
    AA = m_theta.dot(O).dot(A).dot(O.T) - O.dot(A).dot(O.T).dot(m_theta) + gamma.dot(np.squeeze(B.T)) * O.T
    vhat = AA - AA.T
    #cc = -np.sum(np.diag(A.dot(O.T).dot(m_theta).dot(O) + gamma.dot(B.T).dot(O.T)))
    l, grad = FBLLH(theta, gamma, O, A, B, n, method=method, withvol=withvol)
    delta = np.sign(np.sum(AA.dot(vhat)))
    new_O = expm.expm_Higham08(vhat*delta).dot(O)
    new_l, new_grad = FBLLH(theta, gamma, new_O, A, B, n, method=method, withvol=withvol)

    f0 = lambda t0: -FBLLH(theta, gamma, expm.expm_Higham08(vhat*t0).dot(O), A, B, n, method=method, withvol=withvol)[0]

    res = scipy.optimize.minimize_scalar(f0, bounds=[-100.,100.])

    new_O_1 = expm.expm_Higham08(vhat*res.fun).dot(O)

    new_l_1, new_grad_1 = FBLLH(theta, gamma, new_O_1, A, B, n, method=method, withvol=withvol)

    if new_l > new_l_1:
        new_O = new_O_1
    else:
        print "vhat is too small?"

    return new_O, new_l_1, new_grad_1, AA

def grad_mle_update_optim_sqrt(theta, gamma, O, A, B, n, method='hg', withvol=False, dmax=1.0):
    p = len(theta)
    l, grad = FBLLH(theta, gamma, O, A, B, n, method=method, withvol=withvol)
    x = np.sqrt(theta)
    y = np.sqrt(gamma)
    delta = 10.0
    new_x = x*(1.0+delta * grad[:p])
    new_y = y*(1.0+delta * grad[p:p*2])

    # one of the thetas must be zero
    # assuming it's theta[0], let's fix that component of the gradient to zero
    new_x[0] = 0
    grad[0] = 0

    def f(dl):
        x = np.sqrt(theta)
        y = np.sqrt(gamma)
        new_x = x*(1.0+dl * grad[:p])
        new_y = y*(1.0+dl * grad[p:p*2])
        return -FBLLH(new_x**2.0, new_x**2.0, O, A, B, n, method=method, withvol=withvol)[0]

    res = scipy.optimize.minimize_scalar(f, bounds=[0.0, dmax])

    delta = res.fun
    new_x = x*(1.0+delta * grad[:p])
    new_y = y*(1.0+delta * grad[p:p*2])
    new_theta = new_x**2.0
    new_gamma = new_y**2.0

    new_l, new_grad = FBLLH(new_theta, new_gamma, O, A, B, n, method=method, withvol=withvol)

    return new_theta, new_gamma, new_l, new_grad, delta

def MLE(A, B, n, theta=None, gamma=None, O=None, method='hg', withvol=False, max_iter=200, tol=1e-3):

    p = A.shape[0]

    # get the initial estimates for theta, gamma, O if they are not provided
    if theta is None:
        theta = np.arange(p)
    theta = theta - np.amin(theta)
    assert(np.all(theta >= 0))
    if gamma is None:
        gamma = theta + 1
    if O is None:
        O = np.identity(p)

    # try a first update of gamma and theta with n=1 for some reason?
    theta, gamma, l, grad, delta = grad_mle_update(theta, gamma, O, A, B, 1, withvol=withvol)
    # try a second update with n=n now?
    new_theta, new_gamma, new_l, new_grad, new_delta = grad_mle_update(theta, gamma, O, A, B, n, method=method, withvol=withvol)

    a = scipy.linalg.norm(theta - new_theta, ord=1)
    b = scipy.linalg.norm(gamma - new_gamma, ord=1)

    i = 1
    while a+b > tol and i < max_iter:
        a = scipy.linalg.norm(theta - new_theta, ord=1)
        b = scipy.linalg.norm(gamma - new_gamma, ord=1)
        theta, gamma, l, grad, delta = new_theta, new_gamma, new_l, new_grad, new_delta

        O, _, _, _ = grad_mle_orth_update(theta, gamma, O, A, B, n, method=method)
        new_theta, new_gamma, new_l, new_grad, new_delta = grad_mle_update_optim_sqrt(theta, gamma, O, A, B, n, method=method, withvol=withvol)

        i += 1

    return {'theta':new_theta, 'gamma':new_gamma, 'O':O, 'l':new_l, 'grad':new_grad}

def test():
    # astronomy data test
    A = np.zeros((3,3))
    A[0,0]=0.3119
    A[0,1]=0.0292/2.0; A[1,0]=A[0,1]
    A[0,2]=0.0707/2.0;A[2,0]=A[0,2]
    A[1,1]=0.3605
    A[1,2]=0.0462/2.0; A[2,1]=A[1,2]
    A[2,2]=0.3276
    A = 2.0*A - np.diag(np.diag(A))

    B = np.array([-0.0063,-0.0054,-0.0762])[:,None]
    B = -B

    res = MLE(A, B, n=1, theta=None, gamma=None, O=None, method='hg', withvol=True, max_iter=200)

    print res

    res = MLE(A, B, n=1, theta=None, gamma=None, O=None, method='MC', withvol=True, max_iter=200)

    print res

test()

