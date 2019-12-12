#!/usr/bin/env python3

import os, re, sys

first = True
init = False
num_semicol_header = 0

# analyze parameters at the end of the file
def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:
       
        splitpars = line.split(";")

        if len(splitpars) > 1:

            pars[splitpars[0]] = splitpars[1]

    return(pars)

def analyze_data(lines):

    data = [ [ float(celli) ] for celli in lines[0].split(";")[0:-1] ]

    # loop through the lines to collect the data
    for line in lines[1:]:
        splitline = line.split(";")[0:-1]

        for i in range(0,len(splitline)):
            data[i].append(float(splitline[i]))

    # now take averages

    avgs = []
    for i in range(0,len(data)):
        avgs.append(sum(data[i])/len(data[i]))

    return(avgs)

# processes the first line headers
# when making line headers for initial values
def process_first_line(line):

    # get the column names and split them into a list
    line_cols = line.strip().split(";")

    new_cols = ""

    for colname in line_cols:
        if not colname:
            continue

        new_cols += colname + "_t_0;" 

    return(new_cols)

def analyze_file(filename):

    # whether this is the first file we parse
    # (in which case we need to print headers)
    global first;

    # how many semicolons to print the header?
    global num_semicol_header;

    global init;

    # indicator variable whether we are at first line
    # of the file to be read
    firstline = False

    flhead = ""

    lc = 0

    # indicator variable whether we are at the part
    # involving parameters
    parameter_part = False

    parameter_lines = []

    # the line where the parameter output
    # starts
    parline = -1

    # store the last line of data
    last_data_line = ""

    # store the first line of data
    # in case we need initial values too
    first_data_line =""

    # the header of the resulting data file
    flhead = ""

    # open the file and read the stuff
    with open(filename) as infile:
        for lineno, line in enumerate(infile):

            if re.search("^(generation|time|timestep);",line.lower()) is not None:
                firstline = True

            # see whether we also have to store the initial values
            if init:
                if firstline:

                    flhead += process_first_line(line)

            # update line count
            lc += 1

            # get the first line of data
            if lc == 2:
                first_data_line = line.strip()

            # if this is the first line store
            # the header
            if firstline:
                the_line = line.strip()
                assert(len(the_line.split(";")) > 0)
                flhead += the_line
                firstline = False

            # if this is any other line starting
            # with a numerical value store the line
            # as it might be potentially the last one
            elif re.match("^\d",line):
                last_data_line = line

            # hold this until we have the parameter file
            if not parameter_part:
                if lineno > 3 and re.match("^\n",line) is not None:
                    parline = lc
                    parameter_part = True
                    parameter_lines += [line.strip()]
            elif parameter_part:
                parameter_lines += [line.strip()]

    if parline < 1:
        return

    parameters = analyze_parameters(parameter_lines)

    # prepare the data to be printed
    data = ""

    if init:

        # do error checking in terms of the csv values
        count_semicol = len(re.findall(";", first_data_line))
        count_semicol_data = len(re.findall(";", last_data_line))

        assert(count_semicol == count_semicol_data)

        data += first_data_line.strip()

    data += last_data_line.strip()

    if first:
        #check whether header is present
        assert(len(flhead) > 0)

        header_line = ";".join(parameters.keys()) + ";" + flhead.strip() + "file"

        # count number of occurrences of semicolon for error checking
        num_semicol_header = len(re.findall(";",header_line))

        print(header_line)

        first = False

    data_line = ";".join(parameters.values()) + ";" + data +  filename

    print(data_line)

    # count number of occurrences of semicolon for error checking
    num_semicol_data = len(re.findall(";",data_line))

    assert(num_semicol_header == num_semicol_data)


if len(sys.argv) > 2:

    # get initial values from the files as well
    init = True


# run the function on the indicated dir
for root, dir, files in os.walk(sys.argv[1]):

    for file in files:
        if re.search("(sim|iter).*\d$",file) is not None:
            analyze_file(os.path.join(root, file))
