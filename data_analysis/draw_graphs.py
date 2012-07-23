'''
This script is part of the mRF experiment.
It takes a file which concatenates dot scripts as an input
and outputs a PNG for each dot script (probably a mRF graph).
(you must set input_filename to your input filename in the main function)

It can also convert the generated PNGs into a video thought the function
generate_movie() and adds the main neurons inputs and outputs as well as
their connections to the inputs/outputs in every chip with add_input_output_edges ()

The main goal of this script is to either: 
* analyze an entire generation,
* analyze the evolution of the best individual,
* see how the same individual evolved over generations.

WARNING: 
Some functions call external binaries.
Therefore, convert.exe (ImageMagick), dot.exe and ffmpeg.exe must be in $PATH!
* convert.exe : http://www.imagemagick.org/script/binary-releases.php#windows
* dot.exe : http://www.graphviz.org/Download..php
* ffmpeg.exe : http://ffmpeg.arrozcru.org/autobuilds/
If you're on Linux, sudo apt-get install them, it should work.
'''

#!/usr/bin/python
# -*- coding: latin-1 -*-

import os
import re

def cut_dot_files(input_filename, output_filename_prefix):
    ''' 
    This function cuts a file containing several graphs in dot format
    into several dot file, one for each graph.
    '''
    file = open(input_filename, 'r')
    graph_count = 0
    file_output = open(output_filename_prefix + str(graph_count).zfill(4)
                       + '.dot', 'w')

    for cur_line in file:
        # If we are at the beginning of a graph, clear the characters 
        # before the string 'digraph G {'
        if (cur_line.find('digraph G {') >= 0):
            cur_line = 'digraph G {'
        file_output.write(cur_line)
        # If the line contains the character '}', it means we reach the 
        # end of a dot script
        if (cur_line.find('}') >= 0):
            file_output.close()
            graph_count += 1
            file_output = open(output_filename_prefix
                          + str(graph_count).zfill(4) + '.dot', 'w')

    file.close()
    return graph_count

def add_input_output_edges(nb_dot_files, input_filename_prefix, output_filename_prefix):
    ''' 
    This function adds the main neurons inputs and outputs as well as their
    connections to the inputs/outputs in every chip.
    '''
    # for each dot file
    for num_dot_file in range(nb_dot_files):
        file = open(output_filename_prefix + str(num_dot_file).zfill(4)
                    + '.dot', 'r')
        file_output = open(output_filename_prefix + "io" 
                           + str(num_dot_file).zfill(4) + '.dot', 'w')
        # Find the number of inputs/outputs/chips
        max_i = 0
        max_o = 0
        max_nb_chips = 4
        max_nb_neurons = 0
        for cur_line in file:
            # Copy each line of the dot script into the new dot script
            # except for the close bracket '}' which means we reach the
            # end of the script, i.e. we need to add lines there 
            if (cur_line.find('}') < 0):
                file_output.write(cur_line)
            
            line = cur_line.split(" ")
            # We don't care about blank line or the first script line
            # containing 'digraph G {'
            if (len(line) == 0 or cur_line.find('digraph G {') >= 0 
                or line[0].strip() == ''):
                continue
            # Search max_i
            if line[0].find("i") >= 0:
                line[0] = line[0].replace("i", "")
                if int(line[0]) > max_i:
                    max_i = int(line[0])
            # Search max_o
            #TODO: Use startswith instead
            if line[0].find("o") >= 0:
                line[0] = line[0].replace("o", "")
                if int(line[0]) > max_o:
                    max_o = int(line[0])
            # Search max_nb_neurons
            if (line[0].find("i") < 0 and line[0].find("o") < 0 
                and cur_line.find('}') < 0 and int(line[0]) > max_nb_neurons):
                max_nb_neurons = int(line[0])
            
            # Search max_nb_chips
            if (cur_line.find('chip=') >= 0):
                regexp_chip = re.compile(r'chip=[0-9]*', re.IGNORECASE)
                results = regexp_chip.findall(cur_line)
                chip_num = int(results[0].replace('chip=', ''))
                if (chip_num > max_nb_chips):
                    max_nb_chips = chip_num
                
            # If the line contains the character '}', it means we reach the 
            # end of a dot script. We can now add the neurons and edges.
            if (cur_line.find('}') >= 0):
                #print("max_i=" + str(max_i) + "max_o=" + str(max_o))
                #print("max_nb_neurons=" + str(max_nb_neurons) + "max_nb_chips=" + str(max_nb_chips))
                
                # Deal with input neurons
                num_inputs = (max_i+1)/max_nb_chips
                for input_num in range(num_inputs):
                    max_nb_neurons += 1
                    # Add main input neuron
                    file_output.write(str(max_nb_neurons) 
                                + ' [label="INPUT #' + str(input_num) + '" shape=Mdiamond]\n')
                    # Add connections between the main input neuron and each chip
                    for chip_num in range(max_nb_chips):
                        file_output.write(str(max_nb_neurons) 
                                          + ' -> i' + str((chip_num * num_inputs) 
                                                          + input_num) 
                                          + '[label=" "]\n')
                        
                # Deal with output neurons
                num_outputs = (max_o+1)/max_nb_chips
                for ouput_num in range(num_outputs):
                    max_nb_neurons += 1
                    # Add main output neuron
                    file_output.write(str(max_nb_neurons) 
                                      + ' [label="OUTPUT #' + str(ouput_num) + '"shape=Msquare]\n')
                    # Add connections between the main output neuron and each chip
                    for chip_num in range(max_nb_chips):
                        file_output.write('o' 
                                          + str((chip_num * num_outputs) + ouput_num) 
                                          + ' -> ' + str(max_nb_neurons) 
                                          + '[label=" "]\n')
                # End the dot script
                file_output.write('}')
        file_output.close()
        file.close()


def divide_in_subgraphs(filename, output_filename):
    ''' 
    This function divide the graph into several subgraphs, so as to make
    mRF's clusters easily visible. 
    WARNING: the dot file must have the following structure:
    digraph G { 
    Node declarations. Chip number must be in ascending order.
    Edge declarations.
    '''
    # open input and output files
    file = open(filename, 'r')
    file_output = open(output_filename, 'w')
    
    # Find the number of inputs/outputs/chips
    current_chip_number = 0
    node_declaration = False
    for cur_line in file:
        # If this is the first line
        if (cur_line.find('digraph G {') >= 0):
            file_output.write('digraph G {\n\n')
            file_output.write('subgraph cluster0 {\n')
            file_output.write('label="cluster0"\n')
            file_output.write('color=blue\n')
            file_output.write('bgcolor=azure\n')
            file_output.write(cur_line.replace('digraph G {', ''))
            node_declaration = True
            continue
        # If we're in the node declaration part, check out the chip number
        # and create subgraph if we encounter a new chip number
        if (cur_line.find('chip=') >= 0):
            # get current chip number
            regexp_chip = re.compile(r'chip=[0-9]*', re.IGNORECASE)
            results = regexp_chip.findall(cur_line)
            line_chip_number = int(results[0].replace('chip=', ''))
            # chip number can only increase
            assert(current_chip_number <= line_chip_number)
            if (line_chip_number > current_chip_number):
                current_chip_number = line_chip_number
                file_output.write('}\n\n')
                file_output.write('subgraph cluster' + str(current_chip_number) +' {\n')
                file_output.write('label="cluster' + str(current_chip_number) +'"\n')
                file_output.write('color=blue\n')
                file_output.write('bgcolor=azure\n')
                file_output.write(cur_line)
                continue
        # If we reach the end of the node declaration part, then end the last
        # subgraph with "}" , set flag node_declaration to false.
        if ((cur_line.find('chip=') < 0) and node_declaration):
            node_declaration = False
            file_output.write('}\n\n')
            file_output.write(cur_line)
            continue

        file_output.write(cur_line)
            
    # close input and output files
    file_output.close()
    file.close()
        
        
def generate_png(nb_dot_files, input_filename_prefix, output_filename_prefix):
    ''' 
    This function converts a list of dot files into graphs (1 file = 1 graph).
    WARNING: dot.exe must be in $PATH!
    '''
    #for num_dot_file in range(nb_dot_files):
    for num_dot_file in range(3, nb_dot_files):
        print("Generating graph#" + str(num_dot_file) + " out of "
              + str(nb_dot_files - 1))
        os.system("dot -Tpng -o" + output_filename_prefix
                  + str(num_dot_file).zfill(4) + ".png -Kdot "
                  + input_filename_prefix + str(num_dot_file).zfill(4) + ".dot")
        

def generate_movie(nb_dot_files, input_filename_prefix, output_filename_prefix):
    ''' 
    This function converts a list of PNG into a mp4 video
    WARNING: ffmpeg.exe and convert.exe (ImageMagick) must be in $PATH!
    '''
    for num_dot_file in range(nb_dot_files):
        # adding the character '!' to the size makes convert ignore the 
        # aspect ratio and distort the image so it always generates an
        # image exactly the size specified.
        os.system("convert out" + str(num_dot_file).zfill(4)
                  + ".png -resize 1800x900! out" + str(num_dot_file).zfill(4)
                  + "_reduce.png")
        print("Converting graph#" + str(num_dot_file) + " out of "
              + str(nb_dot_files - 1))
    # ffmpeg options:
    # -b : video bit rate
    # -i : input files, %04d says that we have four numbers in the
    #      filename where the number is filled with zeros left of it.
    # -qscale 5 : define fixed video quantizer scale (VBR) where 1 is 
    #             the best and 31 the worst. Since mpeg/jpeg has problems
    #             to compress line graphics it's a good idea to set this 
    #             variable close to 1. You get a big movie file, but 
    #             otherwise the movie doesn't look, well, that good.
    # -r : set the frame rate
    # -y : overwrite output files.
    os.system("ffmpeg -qscale 3 -y -r 1 -b 9600 -i out%04d_reduce.png movie.mp4")


def main():
    '''
    This is the main function
    '''
    input_filename = 'test10.dot'
    output_filename_prefix = 'out'
    nb_dot_files = cut_dot_files(input_filename, output_filename_prefix)
    add_input_output_edges(nb_dot_files, output_filename_prefix, output_filename_prefix)
    for num_dot_file in range(nb_dot_files):
        print("Clustering graph#" + str(num_dot_file) + " out of "
              + str(nb_dot_files - 1))
        filename = output_filename_prefix + "io" + str(num_dot_file).zfill(4) + '.dot'
        output_filename = output_filename_prefix + "ios" + str(num_dot_file).zfill(4) + '.dot'
        divide_in_subgraphs(filename, output_filename)
    generate_png(nb_dot_files, output_filename_prefix + "ios", output_filename_prefix)  
    #generate_movie(nb_dot_files, output_filename_prefix, output_filename_prefix)
    print("ok")


if __name__ == "__main__":
    main()
