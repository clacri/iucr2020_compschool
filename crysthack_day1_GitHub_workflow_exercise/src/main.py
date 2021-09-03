from glob import glob
import random
from nameparser import HumanName
from pathlib import Path


def read_names(path):
    """
    Reads the first token from the first line of all files in the given path

    :param path: path to the directory you want to read the files from
    :return names: a list of the names read from the files
    """
    first_names = []
    for name_file in path.glob("*.txt"):
        with open(name_file, 'r') as f:
            first_line = f.readline()
            name = HumanName(first_line)
            first_names.append(name.first)
    return first_names


GREETINGS = ["Hello", "Hi", "Howdy", "Hola", "Ciao", "Ahoj", "Buna", "G'day", "What's up", "Howzgarn"]


def greet_random(names):
    print(f"{random.choice(GREETINGS)}, {random.choice(names)}!")


def greet_all(names):
    for name in names:
        print(f"{random.choice(GREETINGS)}, {name}!")


def greet_each_other(names):
    name1, name2 = random.sample(names, 2)
    print(f"{random.choice(GREETINGS)}, my name is {name1}, what's your name?")
    print(f"{random.choice(GREETINGS)}, my name is {name2}, How are you enjoying this course?")
    print("I am loving it, but I wish the conference was in-person and we could all be there...")


def main():
    # TODO: build a CLI to allow the path to be set when running the program
    p = Path("../names/")

    names = read_names(p)
    greet_random(names)
    greet_all(names)
    greet_each_other(names)


if __name__ == '__main__':
    main()
