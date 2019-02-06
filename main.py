def fit_linear(filename):
    # this is the main function of this program.
    # first it checks if the file is valid and can be used.
    # in order to check if that's true, it will use a designated function.
    # that function accepts the file, and turns it into usable data, under specific [x,dx,y,dy] values.
    # if that function returns an error with the file, it will be printed as a string
    if type(input_file_function(filename)) == type("string"):
        return print(input_file_function(filename))
    # if the file is valid, the program will now calculate all the values needed.
    (a, da, b, db, chi2, chi2_reduced) = (get_a(filename), get_da(filename), get_b(filename),
                                          get_db(filename), get_chi2(filename), get_chi2_reduced(filename))
    # after acquiring those values, the main program uses them to build a graph
    # and displays them in this format:
    return print("a =", a, "+-", da, '\n'
                 "b =", b, "+-", db, '\n'
                 "chi2 =", chi2,     '\n'
                 "chi2_reduced =", chi2_reduced), plot_function(filename)


def input_file_function(file):
    # this function accepts the file and reads it
    # it separates it according to the file sorting of data (rows/columns)
    # this function will not process the data if it has:
    # (a) (rows / columns) does not have the same length.
    # (b) the uncertainty (rows / columns) are less than zero.
    # (c) there are duplicate value names e.g [x,dx,dx,dy]
    # --------------------------------------------------------
    # the function opens the file and receives the information
    file_pointer = open(file, "r")
    data_raw = file_pointer.readlines()
    # it creates two lists, one for the data itself, one for the axis names.
    working_file = []
    axis_names = []
    # in order to fill the new lists, it determines if the data is presented in rows or columns
    # the first row of the data is examined for more than 1 value names in it.
    # if that is true, the data is sorted in columns, otherwise, in rows.
    row_or_column_variable = data_raw[0].lower().strip("\n").split()
    is_the_first_place_a_string = "x" in row_or_column_variable[0] or "y" in row_or_column_variable[0]
    is_the_second_place_a_string = "x" in row_or_column_variable[1] or "y" in row_or_column_variable[1]
    # if the data received is sorted in columns, the function turns it into lists
    # at this point it also extracts the axis names and puts them in their respective list.
    if is_the_first_place_a_string and is_the_second_place_a_string:
        unsorted_rows = []
        for row in data_raw:
            line = row.lower().strip("\n").split()
            if len(line) == 0:  # this makes sure there are no blank data lines
                continue
            if "axis:" in line:
                axis_names.append(line)
                continue
            unsorted_rows.append(line)
        # to check condition (a)
        # the length of the first column is compared to the lengths of the other columns.
        number_of_columns_not_in_the_same_length = 0
        for j in range(len(unsorted_rows)):
            if len(unsorted_rows[0]) != len(unsorted_rows[j]):
                number_of_columns_not_in_the_same_length += 1
        if number_of_columns_not_in_the_same_length > 0:
            # to which it returns an error
            return "Input file error: Data lists are not the same length."
        # the columns are sorted in lists under each respective variable name [x,dx,y,dy]
        for column_index in range(len(unsorted_rows[0])):
            columns = []
            for row_index in range(len(unsorted_rows)):
                columns.append(unsorted_rows[row_index][column_index])
            working_file.append(columns)
    else:  # the data is sorted in rows.
        for row in data_raw:
            line = row.lower().strip("\n").split()
            if len(line) == 0:  # this makes sure there are no blank data lines
                continue
            if "axis:" in line:
                axis_names.append(line)
            if line[0] == "x" and line not in axis_names:
                working_file.append(line)
            if line[0] == "dx" and line not in axis_names:
                working_file.append(line)
            if line[0] == "y" and line not in axis_names:
                working_file.append(line)
            if line[0] == "dy" and line not in axis_names:
                working_file.append(line)
    # check condition (c)
    x_counter, dx_counter, y_counter, dy_counter = (0, 0, 0, 0)
    for number in range(len(working_file)):
        if working_file[number][0] == "x":
            x_counter += 1
        if working_file[number][0] == "dx":
            dx_counter += 1
        if working_file[number][0] == "y":
            y_counter += 1
        if working_file[number][0] == "dy":
            dy_counter += 1
    if x_counter > 1 or y_counter > 1 or dx_counter > 1 or dy_counter > 1:
        return "Input file error: Duplicate value names "
    # to make the data usable - the function sorts for [x,dx,y,dy]
    # and it convert numbers inside the file into floating point numbers.
    data = []
    for line in working_file:
        data_temp = []
        if line[0] == "x":
            line.pop(0)
            for number in line:
                data_temp.append(float(number))
            data.append(data_temp)
    # to uphold condition (b)
    # check that the uncertainties are bigger than zero.
    for line in working_file:
        data_temp = []
        if line[0] == "dx":
            line.pop(0)
            for number in line:
                if float(number) < 0:
                    return "Input file error: Not all uncertainties are positive."
                data_temp.append(float(number))
            data.append(data_temp)
    for line in working_file:
        data_temp = []
        if line[0] == "y":
            line.pop(0)
            for number in line:
                data_temp.append(float(number))
            data.append(data_temp)
    # to uphold condition (b)
    # check that the uncertainties are bigger than zero.
    for line in working_file:
        data_temp = []
        if line[0] == "dy":
            line.pop(0)
            for number in line:
                if float(number) < 0:
                    return "Input file error: Not all uncertainties are positive."
                data_temp.append(float(number))
            data.append(data_temp)
    # it ends with making sure that after handling the data
    # the rows of data are in the same length.
    # hence condition (a) is reviewed again
    if len(working_file[0]) != len(working_file[1]) \
            or len(working_file[0]) != len(working_file[2]) \
            or len(working_file[0]) != len(working_file[3]):
        return "Input file error: Data lists are not the same length."
    return data, axis_names


def Z_GaG(variable, filename):
    # in order to create the equations for the values we look for,
    # a new mathematical operation has to be determined
    # according to formula (6)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    numerator = 0
    for i in range(len(dy)):
        numerator += (variable[i]/((dy[i])**2))
    denominator = 0
    for j in range(len(dy)):
        denominator += 1/((dy[j])**2)
    z_gag = numerator/denominator
    return z_gag

# These 6 following functions get the values [a,da,b,db,chi2,chi2_reduced] respectively.


def get_a(filename):
    # according to formula (4)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    x_times_y = []
    x_squared = []
    for i in range(len(x)):
        x_times_y.append(x[i]*y[i])
        x_squared.append((x[i])**2)
    numerator_a = (Z_GaG(x_times_y, filename))-(Z_GaG(x, filename)*Z_GaG(y, filename))
    denominator_a = (Z_GaG(x_squared, filename)) - (Z_GaG(x, filename))**2
    a = numerator_a / denominator_a
    return a


def get_da(filename):
    # according to formula (4)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    dy_squared = []
    x_squared = []
    N = len(x)  # number of samples
    for i in range(len(dy)):
        dy_squared.append((dy[i])**2)
        x_squared.append((x[i])**2)
    da_squared = (Z_GaG(dy_squared, filename))/(N*(Z_GaG(x_squared, filename)-(Z_GaG(x, filename))**2))
    da = da_squared**0.5
    return da


def get_b(filename):
    # according to formula (5)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    b = Z_GaG(y, filename) - (get_a(filename))*Z_GaG(x, filename)
    return b


def get_db(filename):
    # according to formula (5)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    dy_squared = []
    x_squared = []
    N = len(x)  # number of samples
    for i in range(len(dy)):
        dy_squared.append((dy[i])**2)
        x_squared.append((x[i])**2)
    db_squared = (Z_GaG(dy_squared, filename)*Z_GaG(x_squared, filename)) /\
                 (N*(Z_GaG(x_squared, filename)-(Z_GaG(x, filename))**2))
    db = db_squared**0.5
    return db


def get_chi2(filename):
    # according to formula (3)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    chi_squared = 0
    a = get_a(filename)
    b = get_b(filename)
    for i in range(len(x)):
        chi_squared += ((y[i]-(a*x[i]+b))/dy[i])**2
    return chi_squared


def get_chi2_reduced(filename):
    # according to formula (7)
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    chi2 = get_chi2(filename)
    N = len(x)  # number of samples
    chi2_reduced = chi2/(N-2)
    return chi2_reduced


def plot_function(filename):
    # this function creates the plot for the data, using "matplotlib"
    # to handle the data, we assign each row to its respective variable
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    (a, da, b, db, chi2, chi2_reduced) = (get_a(filename), get_da(filename), get_b(filename),
                                          get_db(filename), get_chi2(filename), get_chi2_reduced(filename))
    # to create the trendline, "y_axis" values for the fit will be calculated and put into a list.
    linear_fit = []
    for i in range(len(x)):
        linear_fit.append(a*x[i]+b)
    # axis names are extracted from their list, in order to become the titles of the plot.
    x_axis = ""
    y_axis = ""
    for k in range(len(axis_names)):
        if "x" in axis_names[k]:
            x_axis = axis_names[k][2] + axis_names[k][3]
        elif "y" in axis_names[k]:
            y_axis = axis_names[k][2] + axis_names[k][3]
    # Module imported to create the plot:
    import matplotlib.pyplot as plt
    plt.plot(x, linear_fit, "r-")
    plt.errorbar(x, y, xerr=dx, yerr=dy, fmt='b+')
    plt.xlabel(x_axis.title())
    plt.ylabel(y_axis.title())
    plt.savefig("linear_fit.svg", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format="SVG",
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    return plt.show()

# -------------------------------------------------- bonus section -------------------------------------------------- #


def search_best_parameter(filename):
    # this is the main program for this part.
    # it uses the first part of the program to obtain the data and axis names
    # already sorted in lists.
    # it then proceeds to acquire the guesses for the fit parameters "a" and "b"
    # using a designated function to do so.
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    (a, b, da, db) = a_and_b_extractor(filename)
    # it later calculates the chi squared values for all parameter guesses.
    # for all those values it picks the minimal chi value, and presents the value alongside
    # the "a" and "b" parameters who suit that value.
    (best_a, best_b, chi_2_min) = chi_squared_for_non_linear_functions(filename)
    # the function will return values "a", "b", "chi2", and a plot of the fit.
    # it will also plot the distribution function of chi squared values to each "a"
    # parameter holding the "b" parameter constant to its optimal value.
    return print("a =", best_a, "+-", da, '\n'
                 "b =", best_b, "+-", db, '\n'
                 "chi2 =", chi_2_min), \
                 plot_that_bonus_graph(filename), \
                 plot_that_bonus_graph_but_the_interesting_one(filename)


def a_and_b_extractor(filename):
    # this is the designated function for getting the values of "a" and "b"
    # the text file will have initial anf final guesses for those values, including a step size.
    # this function receives this data, and creates a list for that range.
    file_pointer = open(filename, "r")
    data_raw = file_pointer.readlines()
    a_and_b_dict = {}
    # it will find the rows that have "a" and "b" in them
    for row in data_raw:
        line = row.lower().strip("\n").split()
        if "a" in line or "b" in line:
            a_and_b_dict[line[0]] = line[1:4]
    # to access the dictionary values more easily, it will create variables for each key.
    a = a_and_b_dict["a"]
    b = a_and_b_dict["b"]
    for i in range(len(a)):
        a[i] = float(a[i])
        b[i] = float(b[i])
    # these lists include the range from initial guess to the final guess with the step size.
    a_list = []
    b_list = []
    for i in range(int(abs(a[0] - a[1]) / abs(a[2])) + 1):
        a_list.append(a[0] + i * a[2])
    for n in range(int(abs(b[0] - b[1]) / abs(b[2])) + 1):
        b_list.append(b[0] + n * b[2])
    # it returns these lists, accompanied by the step sizes,
    # which will later become the uncertainties for the parameters.
    return a_list, b_list, a[2], b[2]


def chi_squared_for_non_linear_functions(filename):
    # this function retrieves the data gathered by th other functions
    # and calculates the minimal chi squared value for a specific pair of "a" and "b"
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    (a, b, da, db) = a_and_b_extractor(filename)
    # according to formula (1):
    chi_2 = {}
    for n in range(len(b)):
        for k in range(len(a)):
            # while holding "a" and "b" constant, the loop goes over all the [x,y] data set.
            # that means there is a chi squared value for each set of parameters.
            chi_2_i = 0
            for i in range(len(x)):
                numerator = y[i] - (a[k]*x[i] + b[n])
                denominator_squared = (dy[i])**2 + (a[k]*(x[i]+dx[i])+b[n]-(a[k]*(x[i]-dx[i])+b[n]))**2
                denominator = denominator_squared**0.5
                chi_2_i += (numerator/denominator)**2
            # the chi squared is now the key for both parameters.
            chi_2[chi_2_i] = [a[k], b[n]]
    chi_2_min = min(chi_2.keys())
    (best_a, best_b) = chi_2[chi_2_min]
    return best_a, best_b, chi_2_min


def plot_that_bonus_graph(filename):
    # this function plots the linear fit for the data, in a similar way to the last part of the project.
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    (best_a, best_b, chi_2_min) = chi_squared_for_non_linear_functions(filename)
    # the linear fit is constructed from the "a" and "b" who produce the minimal chi value.
    linear_fit = []
    for i in range(len(x)):
        linear_fit.append(best_a * x[i] + best_b)
    # axis names are extracted from their list, in order to become titles in the plot
    x_axis = ""
    y_axis = ""
    for k in range(len(axis_names)):
        if "x" in axis_names[k]:
            x_axis = axis_names[k][2] + axis_names[k][3]
        elif "y" in axis_names[k]:
            y_axis = axis_names[k][2] + axis_names[k][3]
    import matplotlib.pyplot as plt
    plt.plot(x, linear_fit, "r-")
    plt.errorbar(x, y, xerr=dx, yerr=dy, fmt='b+')
    plt.xlabel(x_axis.title())
    plt.ylabel(y_axis.title())
    plt.savefig("linear_fit.svg", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format="SVG",
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    return plt.show()


def plot_that_bonus_graph_but_the_interesting_one(filename):
    # this function plots the chi values vs. the "a" values
    # while holding "b" to its optimal value.
    (data, axis_names) = input_file_function(filename)
    (x, dx, y, dy) = data
    (a, b, da, db) = a_and_b_extractor(filename)
    (best_a, best_b, chi_2_min) = chi_squared_for_non_linear_functions(filename)
    # it creates a list for all chi values, in the order of the "a" values that match them.
    chi_2 = []
    for k in range(len(a)):
        chi_2_i = 0
        for i in range(len(x)):
            numerator = y[i] - (a[k] * x[i] + best_b)
            denominator_squared = (dy[i]) ** 2 + (a[k] * (x[i] + dx[i]) + best_b - (a[k] * (x[i] - dx[i]) + best_b)) ** 2
            denominator = denominator_squared ** 0.5
            chi_2_i += (numerator / denominator) ** 2
        chi_2.append(chi_2_i)
    import matplotlib.pyplot as plt
    plt.plot(a, chi_2, "b-")
    plt.xlabel("a")
    plt.ylabel("chi2(a,b=" + str(best_b)+")")
    plt.savefig("numeric_sampling.svg", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format="SVG",
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    return plt.show()

# -------------------------------------------------- end project -------------------------------------------------- #


filename_1 = 'input_bonus.txt'
# search_best_parameter(filename_1)
filename_2 = 'input1.txt'
fit_linear(filename_2)

