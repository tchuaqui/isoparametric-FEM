import sys, os
import matplotlib.pyplot as plt


def run(cpp_file, exe_file):
    os.system("echo Compiling " + cpp_file)
    os.system('g++ ' + cpp_file + ' -o' + exe_file)
    os.system("echo Running " + exe_file)
    os.system('./' + exe_file)

if __name__=='__main__':
    cpp_file = 'main.cpp'
    exe_file = 'out'
    run(cpp_file, exe_file)
    with open('initial_coordinates.txt', 'r') as f:
        initial = [[float(num) for num in line.split()] for line in f]
    with open('final_coordinates.txt', 'r') as f:
        final = [[float(num) for num in line.split()] for line in f]
    for i in range(len(final)):
        plt.scatter(initial[i][0], initial[i][1], marker='o', s=1, color='blue')
        plt.scatter(final[i][0], final[i][1], marker='o', s=1, color='red')
    plt.legend(('Undeformed','Deformed'))        
    plt.show()



