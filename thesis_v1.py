# Calculating contaminant concentration at a given point
# using the Analytical Contaminant Transport Equation
# in a 2-dimensional environment.

# Velocity_time signifies the distance a contaminant has travelled (meters)
# after an unknown length of time (days).


import numpy as np
import pandas as pd
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

    x_data = []
    y_data = []
    con_data = []
    tri_data = []

    for dist_x in range(2, 450, 2):
        for dist_y in range(0, 201):
            well_1 = Well(ft_to_m(dist_x), ft_to_m(dist_y))
            well_con = calc_concentration(source, material, well_1, vt)
            x_data.append(dist_x)
            y_data.append(dist_y)
            con_data.append(well_con)
            tri_data.append((dist_x, dist_y, well_con))

    data = (x_data, y_data, con_data), tri_data

    return data


def reflect_over_y_axis(d):
    """Reflects coordinates and concentrations over y axis to get full data set."""
    x_data, y_data, z_data = d
    inverted_y = []
    for value in y_data:
        inverted_y.append(-value)

    return x_data + x_data, y_data + inverted_y, z_data + z_data


def get_all_coord(data):
    """Takes concentration lists needing to be reflected.
    Positive and negative y spread are assumed equivalent."""
    coord = reflect_over_y_axis(data)

    return coord


# Plot Creation #


def create_plot(coordinates, specific_point):
    """Creates plot of concentration amounts."""
    all_x, all_y, all_z = coordinates
    data = []
    for i, x in enumerate(all_x):
        data.append((x, all_y[i], all_z[i]))

    contour_data = pd.DataFrame.from_records(data, columns=['x', 'y', 'z'])
    Z = contour_data.pivot_table(index='x', columns='y', values='z').T.values

    X_unique = np.sort(contour_data.x.unique())
    Y_unique = np.sort(contour_data.y.unique())
    X, Y = np.meshgrid(X_unique, Y_unique)

    fig1, plot = plt.subplots()

    plt.title('Contaminant Concentration (\u03BCg/L)', size='x-large')

    colors = ('olivedrab', 'goldenrod', 'chocolate', 'crimson')
    con_plot = plot.contour(X, Y, Z, levels=(100, 300, 500, 900), colors=colors)
    manual_locations = [(50, 50), (100, 100), (160, 125), (200, 250)]
    plt.clabel(con_plot, inline=1, fontsize=10, manual=manual_locations)

    plot.set_xlim(right=500)
    plot.set_xlabel('Horizontal Distance (ft)')
    plot.set_ylabel('Lateral Distance (ft)')

    if specific_point is not None:
        plt.scatter(specific_point[0], specific_point[1], color='black', alpha=0.7, marker='+')

    fig1.set_size_inches(8, 6)
    plt.show()


# Finding Specific Concentration #


def get_specific_concentration(data):
    """Returns concentration at requested x,y position."""
    while True:
        user_x = int(input('What is the x-value for the point you would like between 2 and 400? '
                           'Please select an even number. '))
        if user_x % 2 == 0 and 2 <= user_x <= 400:
            break
        else:
            print('That is not a valid location.')

    while True:
        user_y = int(input('What is the y-value for the point you would like between -200 and 200? '))
        if -200 <= user_y <= 200:
            break
        else:
            print('That is not a valid location.')

    user_con = 0
    for pair in data:
        if user_x == pair[0] and abs(user_y) == pair[1]:
            user_con = pair[2]
    print('The concentration at point {}, {} is {} \u03BCg/L.'.format(user_x, user_y, user_con))

    return user_x, user_y


# Main #


def main():
    """Establish variables and run program."""
    source_width = 200            # given in feet
    source_concentration = 1000   # given in micro-gram/Liter
    material_dispersivity_x = 10  # given in feet
    material_dispersivity_y = 6   # given in feet
    velocity_time = ft_to_m(350)  # given in feet

    source_1 = Source(ft_to_m(source_width), source_concentration)
    material_1 = Material(ft_to_m(material_dispersivity_x), ft_to_m(material_dispersivity_y))

    positive_y_data = get_data(source_1, material_1, velocity_time)
    all_data = get_all_coord(positive_y_data[0])

    specific = input('Would you like the concentration for a specific point? (Y/N) ').upper()
    if 'Y' in specific:
        specific_point = get_specific_concentration(positive_y_data[1])
    else:
        specific_point = None

    print('\nThe concentration plot has been displayed.')
    print('The program will terminate upon closure of the plot window.')
    create_plot(all_data, specific_point)


if __name__ == '__main__':
    main()
