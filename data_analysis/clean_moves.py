'''
This script is called by draw_moves_stats.m
It cleans all_phenotypes_moves.dat the log file.
'''


import sys
def main():
    '''
    This is the main function
    '''
    print("hi")
    
    # Parse arguments
    if (len(sys.argv) != 2 ):
        print 'ERROR: you must specify 1 argument'
        return
    input_filename = str(sys.argv[1]) # Eg. "exp02"
    
    # Read file and clean it
    file = open(input_filename + '.dat', 'r')
    file_output = open(input_filename + "clean"  + '.dat', 'w')
    
    for cur_line in file:
        if (cur_line.find('NewIndividual') >= 0):
            continue
        line = cur_line.split(" ")
        if (cur_line.split(",") < 12):
            continue
        file_output.write(cur_line)
    
    file.close()
    file_output.close()

if __name__ == "__main__":
    main()