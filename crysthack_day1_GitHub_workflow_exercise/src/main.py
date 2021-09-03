from glob import glob
import random
from nameparser import HumanName
from pathlib import Path


def read_names():
    """
    Read the names in ../names

    :return: list of first names
    """
    first_names = []
    p = Path("./names")
    for name_file in p.glob("*"):
        with open(name_file, 'r') as f:
            first_line = f.readline()
            name = HumanName(first_line)
            first_names.append(name.first)

    return first_names


def greet(names):
    greetings = ["Hello", "Hi", "Howdy", "Hola", "Ciao", "Ahoj", "Buna"]
    greeting_choice = random.choice(greetings)
    name_choice = random.choice(names)
    print(f"{greeting_choice}, {name_choice}!")


if __name__ == '__main__':
    greet(read_names())
