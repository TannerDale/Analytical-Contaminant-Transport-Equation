# Calculating contaminant concentration at a given point
# using the Analytical Contaminant Transport Equation
# in a 2-dimensional environment.

# Velocity_time signifies the distance a contaminant has travelled (meters)
# after an unknown length of time (days).


from scipy import special
import matplotlib.pyplot as plt

# Set-up #


class Source:
    def __init__(self, width, concentration):
        self.concentration = concentration  # micro-gram per liter
        self.width = width


class Material:
    def __init__(self, dispersivity_x, dispersivity_y):
        self.ax = dispersivity_x  # meters
        self.ay = dispersivity_y  # meters


class Well:
    def __init__(self, distance_x, distance_y):
        self.dx = distance_x  # meters
        self.dy = distance_y  # meters


# Unit Conversions #


def ft_to_m(ft):
    """Unit conversion for feet to meters."""
    meters = ft / 0.3048

    return meters


# Calculating Directional Spread #


def spread_x_positive(material, well, vel_time):
    """Finds x_value of contaminant at requested distance in x direction to be converted."""
    dist_x = well.dx
    vt = vel_time
    ax = material.ax

    x_value = special.erfc((dist_x - vt) / (2 * ((ax * vt) ** (1 / 2))))

    return x_value


def spread_y_positive(source, material, well):
    """Finds y_value of contaminant at requested distance in positive y direction to be converted."""
    s_width = source.width
    disperse_y = material.ay
    dist_x = well.dx
    dist_y = well.dy

    pos_y_value = special.erf((dist_y + (s_width / 2)) / (2 * ((disperse_y * dist_x) ** (1 / 2))))

    return pos_y_value


def spread_y_negative(source, material, well):
    """Finds y_value of contaminant at requested distance in negative y direction to be converted."""
    s_width = source.width
    disperse_y = material.ay
    dist_x = well.dx
    dist_y = well.dy

    neg_y_value = special.erf((dist_y - (s_width / 2)) / (2 * ((disperse_y * dist_x) ** (1 / 2))))

    return neg_y_value


# Converting Spread to Concentration by Volume #


def convert_to_volume(source, value):
    """Applies error function value to source concentration, returning resulting concentration."""
    s = source.concentration
    concentration = (s / 4) * value

    if concentration <= 0.01:  # volume under 0.01 micro-grams per liter is assumed 0
        return 0
    else:
        return format(concentration, '0.2f')


def calc_concentration(source, material, well, vel_time):
    """Calculates contaminant concentration at given well (point)."""
    spread_x = spread_x_positive(material, well, vel_time)
    spread_y_pos = spread_y_positive(source, material, well)
    spread_y_neg = spread_y_negative(source, material, well)
    to_convert = spread_x * (spread_y_pos - spread_y_neg)
    concentration = convert_to_volume(source, to_convert)

    return float(concentration)


# Calculating Concentrations #


def get_data(source, material, vt):
    """Calculates data in an x,y range and sorts by concentration."""
    over_900 = []
    over_500 = []
    over_300 = []
    over_100 = []

    for dist_x in range(2, 450, 2):
        for dist_y in range(0, 201):
            well_1 = Well(ft_to_m(dist_x), ft_to_m(dist_y))
            well_con = calc_concentration(source, material, well_1, vt)

            if well_con >= 900.0:
                lst = [dist_x, dist_y, well_con]
                over_900.append(lst)
            elif well_con >= 500.0:
                lst = [dist_x, dist_y, well_con]
                over_500.append(lst)
            elif well_con >= 300.0:
                lst = [dist_x, dist_y, well_con]
                over_300.append(lst)
            elif well_con >= 100.0:
                lst = [dist_x, dist_y, well_con]
                over_100.append(lst)

    data = [over_900, over_500, over_300, over_100]

    return data


def reflect_over_y_axis(d):
    """Reflects coordinates and concentrations over y axis to get full data set."""
    use = []
    for coord in d:
        new = [coord[0], -coord[1], coord[2]]
        use.append(new)

    return d + use


def get_all_coord(data):
    """Takes concentration lists needing to be reflected.
    Positive and negative y spread are assumed equivalent."""
    coord = [reflect_over_y_axis(lst) for lst in data]

    return coord


# Plot Creation #


def create_plot(coordinates):
    """Creates plot of concentration amounts."""
    all_x = []
    all_y = []
    for amounts in coordinates:
        x, y, con = zip(*amounts)
        all_x.append(x)
        all_y.append(y)

    colors = ['red', 'orange', 'yellow', 'green']
    labels = ['901+', '501-900', '301-500', '100-300']

    fig, plot = plt.subplots()
    for i in range(4):
        plot.scatter(all_x[i], all_y[i], color=colors[i], label=labels[i])

    fig.set_size_inches(9, 7)
    plot.set_xlabel('Horizontal Distance (ft)', size='large')
    plot.set_ylabel('Lateral Distance (ft)', size='large')
    plt.legend(bbox_to_anchor=(0.934, 0.925), loc='center', fontsize='small')
    plt.title('Contaminant Concentration (\u03BCg/L)', size='xx-large')
    plt.show()


# Find User Concentration #


def get_specific_concentration(coord_x, coord_y, data):
    """Returns concentration at requested x,y position."""
    use_set = [coord_x, coord_y]
    for d in data:
        for lst in d:
            if set(use_set).issubset(lst):
                return lst[2]


# Main #


def main():
    """Establish variables and run program."""
    source_width = 200  # given in feet
    source_concentration = 1000
    material_dispersivity_x = 10  # given in feet
    material_dispersivity_y = 6  # given in feet
    velocity_time = ft_to_m(350)  # given in feet

    source_1 = Source(ft_to_m(source_width), source_concentration)
    material_1 = Material(ft_to_m(material_dispersivity_x), ft_to_m(material_dispersivity_y))

    positive_y_data = get_data(source_1, material_1, velocity_time)
    all_data = get_all_coord(positive_y_data)

    specific = input('Would you like the concentration for a specific point? (Y/N) ').upper()
    if 'Y' in specific:
        user_x = int(input('What is the x-value for the point you would like between 2 and 400? '
                           'Please select an even number. '))
        user_y = int(input('What is the y-value for the point you would like between -200 and 200? '))
        user_con = get_specific_concentration(user_x, user_y, all_data)
        if user_con is None:
            user_con = 0
        print('The concentration at point ({}, {}) is {} \u03BCg/L.'.format(user_x, user_y, user_con))

    print('\nThe concentration plot has been displayed.')
    print('The program will terminate upon closure of the plot window.')
    create_plot(all_data)


if __name__ == '__main__':
    main()
