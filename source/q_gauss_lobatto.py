class QGaussLobatto:
    # Class to calculate and store points and weights for lobatto quadrature
    def __init__(self, n_points):
        # Initialises the class with the provided number of points
        assert n_points >= 2, "Gauss-Lobatto points for totalling less than 2 doesn't make sense"
        assert n_points <= 6, "Currently more than 6 Gauss-Lobatto points are not supported"
        self.n_points = n_points
        if n_points == 2:
            self.q_points = [-1.0, 1.0]
            self.q_weights = [1.0, 1.0]
        elif n_points == 3:
            self.q_points = [-1.0, 0.0, 1.0]
            temp = 1.0/3
            self.q_weights = [temp, 4*temp, temp]
        elif n_points == 4:
            temp = 0.2**0.5
            self.q_points = [-1.0, -temp, temp, 1.0]
            temp = 1.0/6
            self.q_weights = [temp, 5*temp, 5*temp, temp]
        elif n_points == 5:
            temp = (3.0/7)**0.5
            self.q_points = [-1.0, -temp, 0.0, temp, 1.0]
            temp = 49.0/90
            self.q_weights = [0.1, temp, 32.0/45, temp, 0.1]
        else:
            temp = 2*(7.0**0.5)/21
            temp1 = (1.0/3 - temp)**0.5
            temp2 = (1.0/3 + temp)**0.5
            self.q_points = [-1.0, -temp2, -temp1, temp1, temp2, 1.0]
            temp = (14.0 + 7**0.5)/30
            temp1 = (14.0 - 7**0.5)/30
            temp2 = 1.0/15
            self.q_weights = [temp2, temp1, temp, temp, temp1, temp2]