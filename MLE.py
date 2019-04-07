import numpy as np
import scipy
import scipy.linalg
import fb

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

    alpha = np.concatenate(theta, gamma)

    if method == 'hg':
        c_res = fb.hgm_FB_2(alpha, withvol=withvol)
    elif method == 'hg_linear':
        c_res = fb.hgm_FB(alpha, withvol=withvol)
    elif method == 'MC':
        c_res = fb.G_FB(alpha, method='MC', withvol=withvol)
    elif method == 'SPA':
        c_res = fb.saddleapprox_FB_revised(sorted(alpha[:p])+1.0, M=alpha[p:p*2]/2.0, dub=1, order=3) * np.exp(1)
    elif method == 'lpinv':
        raise ValueError('Not sure what lpinv is supposed to do.')
        c_res = lpinv(alpha)
    else:
        raise ValueError('Not a valid method.')

    c = c_res[0]

    l = -n*np.log(c) - np.sum(A.dot(O.T).dot(np.diag(theta)).dot(O) + gamma.dot(B.T).dot(O))

    grad = -n*c_res[1:] / c + np.concatenate((-np.diag(O.dot(A).dot(O.T)),O.dot(B)))

    return l, grad

# I think this is useless in the context of this problem...
def orthogonal_matrix(x):
    a = x[:,None]
    q = np.identity(len(x)) - 2.0*a.dot(a.T)/float(a.T.dot(a))
    return q

def grad_mle_update(theta, gamma, O, A, B, n, withvol=False, min_step=1e-3):
    p = len(theta)

    # get the initial llh and gradient
    l, grad = FBLLH(theta, gamma, O, A, B, n, method=method, withvol=withvol)

    # one of the thetas must be zero
    # assuming it's theta[0], let's fix that component of the gradient to zero
    grad[0] = 0

    # try stepping once
    delta = 1
    new_theta = theta + delta*grad[:,p]
    new_theta[0] = 0
    new_gamma = gamma + delta*grad[p,:]
    new_l, new_grad = FBLLH(new_theta, new_gamma, O, A, B, n, method=method, withvol=withvol)

    # try some smaller steps
    while new_l < l and delta > min_step:
        delta /= 2.0
        new_theta = theta + delta*grad[:p]
        new_theta[0] = 0
        new_gamma = gamma + delta*grad[p,:]
        new_l, new_grad = FBLLH(new_theta, new_gamma, O, A, B, n, method=method, withvol=withvol)

    return new_theta, new_gamma, new_l, new_grad, delta

def grad_mle_orth_update():
    pass

def grad_mle_update_optim_sqrt():
    pass

def MLE(A, B, n, theta=None, gamma=None, O=None, method='hg', withvol=False, max_iter=200):

    p = A.shape[0]

    # get the initial estimates for theta, gamma, O if they are not provided
    if theta is None:
        theta = np.range(p)
    theta = theta - np.amin(theta)
    assert(np.all(theta >= 0))
    if gamma is None:
        gamma = theta + 1
    if O is None:
        O = np.identity(p)

    # try a first update of gamma and theta with n=1 for some reason?
    theta, gamma, l, grad, delta = grad_mle_update(theta, gamma, O, A, B, 1, withvol=withvol)
    # try a second update with n=n now?
    new_theta, new_gamma, new_l, new_grad, new_delta = grad_mle_update(theta, gamma, O, A, B, n, withvol=withvol)

    a = scipy.linalg.norm(theta - new_theta, ord=1)
    b = scipy.linalg.norm(gamma - new_gamma, ord=1)

    i = 1
    while a+b > tol and i < max_iter:
        a = scipy.linalg.norm(theta - new_theta, ord=1)
        b = scipy.linalg.norm(gamma - new_gamma, ord=1)
        theta, gamma, l, grad, delta = new_theta, new_gamma, new_l, new_grad, new_delta

        O, _ = grad_mle_orth_update(theta, gamma, O, A, B, n)
        new_theta, new_gamma, new_l, new_grad, new_delta = grad_mle_update_optim_sqrt(theta, gamma, O, A, B, n, withvol=withvol)

        i += 1

    return {'theta':new_theta, 'gamma':new_gamma, 'O':O, 'l':new_l, 'grad':new_grad}
