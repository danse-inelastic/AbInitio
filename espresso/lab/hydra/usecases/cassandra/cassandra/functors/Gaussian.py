

from pyre.components.Component import Component


class Gaussian(Component):
    """
    Component that implements the normal distribution with mean mu and variance sigma 

        g(x; mu,sigma) = \frac{1}{\sqrt{2pi} sigma} e^{-\frac{|x-mu|^2}{2sigma^2}}

    mu, and sigma are implemented as component properties so that Gaussian can conform to the
    Functor interface. See quadrature.interfaces.Functor for more details.
    """

    class Inventory(Component.Inventory):
        # public state
        
        import pyre.inventory
        
        mu = pyre.inventory.array("mu")
        mu.doc = "the position of the mean of the Gaussian distribution"
        mu.default = [0.0]
        
        sigma = pyre.inventory.float("sigma")
        sigma.doc = "the variance of the Gaussian distribution"
        sigma.default = 1.0


    def eval(self, points):
        """
        Compute the value of the gaussian
        """

        # access the math symbols
        from math import exp, sqrt
        from math import pi as pi
        # cache the inventory items
        mu = self.inventory.mu
        sigma = self.inventory.sigma
        # precompute the normalization factor
        normalization = 1 / sqrt(2*pi) / sigma
        # loop over points and yield the computed value
        for x in points:
            # compute |x - mu|^2
            # this works as long as x and mu have the same length
            r2 = sum((x_i - mu_i)**2 for x_i, mu_i in zip(x, mu))
            # yield the value at the xurrent x
            yield normalization * exp(- r2 / 2 / sigma**2)

        return
