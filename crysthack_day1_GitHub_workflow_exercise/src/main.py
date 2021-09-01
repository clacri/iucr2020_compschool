import glob
import random


def read_names():
    """
    Read the names in ../names
    :return: names_list
    """
    names_list = []
    for name_file in glob.glob("./names/*"):
        with open(name_file, 'r') as f:
            first_line = f.readline()
            first_name = first_line.split()[0]
            names_list.append(first_name)
    return names_list


def greet(names):
    greeting_list = ["Hello", "Hi", "Howdy"]
    print("%s, %s!" % (random.choice(greeting_list),
                       random.choice(names)))


if __name__ == '__main__':
    greet(read_names())
